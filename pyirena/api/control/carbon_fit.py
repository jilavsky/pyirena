"""pyirena.api.control.carbon_fit — agent-drivable Carbon model fitting.

Full-range SAXS+WAXS analysis of disordered carbonaceous materials, as the
sum of grain Porod scattering, micropore scattering and turbostratic
diffraction (Saurel et al., *Energy Storage Materials* **21** (2019) 162–173).

The surface is deliberately narrow.  The Carbon model has four sections, a
variable peak list and three optional geometry links, and exposing one setter
per control would be forty tools; instead every fittable quantity is addressed
by the same dotted key the core uses — ``background.S_macro``,
``saxs.pore_radius``, ``waxs.delta_z2``, ``peak.002.Q0`` — through one
:func:`set_carbon_parameter` / :func:`set_carbon_parameter_fit` /
:func:`set_carbon_parameter_bounds` trio, and the *shape* of the model is set
through :func:`configure_carbon_model`.  :func:`list_carbon_parameters` is how
an agent discovers which keys currently exist, which depends on the shape.

These functions operate on the *same*
:class:`~pyirena.api.control.session.Session` objects as the other control
surfaces, so the shared session-lifecycle and Q-range tools are reused as-is.
The model state lives in ``session.model`` (a
:class:`~pyirena.core.carbon_fit.CarbonFitModel`) with
``session.model_name == "carbon_fit"``.

There is no per-region fit.  The contrast that scales the background and the
micropore term is computed from the fitted WAXS peak positions and the
porosity, so the three regions are coupled through the physics and
:func:`run_carbon_fit` refines all of them together.

Typical workflow
----------------
>>> import pyirena.api.control as ctrl                      # doctest: +SKIP
>>> sid = ctrl.open_dataset("/data/hard_carbon.h5")["session_id"]
>>> ctrl.select_carbon_model(sid)
>>> ctrl.configure_carbon_model(sid, saxs_mode="teubner_strey",
...                             use_roughness=True)
>>> ctrl.list_carbon_peaks(sid)
>>> ctrl.set_carbon_parameter(sid, "peak.002.Q0", 1.80)
>>> ctrl.run_carbon_fit(sid)
>>> ctrl.get_carbon_results(sid)        # derived quantities, not coefficients
>>> ctrl.get_carbon_fit_image(sid)
>>> ctrl.save_carbon_fit(sid)

All functions return plain dicts.  Errors are
``{"error": ..., "code": ..., "suggestion": ...}`` dicts rather than exceptions.
"""
from __future__ import annotations

from typing import Optional

import numpy as np

from pyirena.api._paths import PathSecurityError, resolve_safe
from pyirena.api.control._images import render_png
from pyirena.api.control.errors import make_error, no_fit, no_session
from pyirena.api.control.session import fit_mask, get_session
from pyirena.api.control.unified_fit import _quality_scalars

__all__ = [
    "list_carbon_options",
    "select_carbon_model",
    "get_carbon_config",
    "configure_carbon_model",
    "list_carbon_parameters",
    "set_carbon_parameter",
    "set_carbon_parameter_fit",
    "set_carbon_parameter_bounds",
    "list_carbon_peaks",
    "add_carbon_peak",
    "remove_carbon_peak",
    "set_carbon_material",
    "run_carbon_fit",
    "get_carbon_results",
    "get_carbon_fit_image",
    "save_carbon_fit",
]


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _require_carbon(session_id: str):
    """Return (session, None) when a Carbon model is ready, else (None, error)."""
    s = get_session(session_id)
    if s is None:
        return None, no_session(session_id)
    if s.model is None or s.model_name != "carbon_fit":
        return None, make_error(
            f"Session '{session_id}' has no Carbon model.",
            suggestion="Call select_carbon_model() first.",
            code="NO_CARBON_MODEL",
        )
    return s, None


def _find_ref(model, key: str):
    """The parameter reference for a dotted key, or None if it is not active."""
    for ref in model.parameter_refs(active_only=True):
        if ref.key == key:
            return ref
    return None


def _bad_key(model, key: str) -> dict:
    active = [r.key for r in model.parameter_refs(active_only=True)]
    return make_error(
        f"'{key}' is not an active parameter of this model.",
        suggestion=(
            "Call list_carbon_parameters() for the current keys. Which "
            "parameters exist depends on the model's shape — the SAXS-region "
            "mode, whether roughness and the crumpled envelope are on, the "
            "geometry links, and the peak list. "
            f"Active now: {', '.join(active[:12])}"
            + (" …" if len(active) > 12 else "")
        ),
        code="BAD_PARAM",
    )


def _peak_summary(model) -> list:
    """One row per peak, with its live d-spacing."""
    out = []
    for index, peak in enumerate(model.peaks):
        d = 2.0 * np.pi / peak.Q0 if peak.Q0 > 0 else float("nan")
        out.append({
            "index": index,
            "label": peak.label,
            "enabled": bool(peak.enabled),
            "Q0": float(peak.Q0),
            "d_spacing": float(d),
            "K": float(peak.K),
            "FWHM_G": float(peak.FWHM_G),
            "FWHM_L": float(peak.FWHM_L),
        })
    return out


def _finite(value) -> Optional[float]:
    """A JSON-safe float, or None for NaN/inf — agents choke on NaN."""
    try:
        v = float(value)
    except (TypeError, ValueError):
        return None
    return v if np.isfinite(v) else None


# ---------------------------------------------------------------------------
# Discovery and model creation
# ---------------------------------------------------------------------------

def list_carbon_options() -> dict:
    """The choices that define a Carbon model's shape.

    Returns:
        dict with ``saxs_modes`` and ``waxs_envelopes`` (each key → a
        human-readable label), the switchable options, and the meaning of the
        three geometry links.
    """
    from pyirena.core.carbon_fit import CARBON_SAXS_MODES, CARBON_WAXS_ENVELOPES

    return {
        "ok": True,
        "saxs_modes": dict(CARBON_SAXS_MODES),
        "waxs_envelopes": dict(CARBON_WAXS_ENVELOPES),
        "switches": {
            "use_roughness": "Add a nanoscale surface-roughness Porod term "
                             "(eq. 3) to the grain background.",
            "use_fractal": "Aggregate the micropores with a Teixeira mass "
                           "fractal (eq. 5). Only for saxs_mode='fractal'.",
            "use_orientation_factor": "Apply the 1/Q² powder average to the "
                                      "diffraction term. Leave on.",
            "link_R_to_pore": "Use the SAXS pore radius as the WAXS layer "
                              "transition radius instead of fitting both.",
            "link_D_to_saxs": "Share the fractal dimension between the SAXS "
                              "and WAXS regions.",
            "link_sigma_to_saxs": "Share the fractal cutoff Σ between the two "
                                  "regions.",
        },
        "note": "All three regions are fitted together; the contrast that "
                "scales them is computed from the fitted peak positions.",
    }


def select_carbon_model(session_id: str, formula: str = "C") -> dict:
    """Create a Carbon model on this session, with the usual carbon defaults.

    Args:
        session_id: An open session.
        formula: Chemical formula of the solid, e.g. ``"C"`` or
            ``"C0.95N0.05"`` for a doped carbon.

    Returns:
        The new model's configuration, as :func:`get_carbon_config` reports it.
    """
    from pyirena.core.carbon_fit import CarbonFitModel

    s = get_session(session_id)
    if s is None:
        return no_session(session_id)

    model = CarbonFitModel()
    model.material.formula = str(formula or "C")
    s.model = model
    s.model_name = "carbon_fit"
    s.last_fit_result = None
    return get_carbon_config(session_id)


def get_carbon_config(session_id: str) -> dict:
    """The model's current shape, peaks, material settings and contrast chain."""
    s, err = _require_carbon(session_id)
    if err:
        return err
    m = s.model
    material = m.resolve_material()
    return {
        "ok": True,
        "session_id": session_id,
        "background_enabled": bool(m.background.enabled),
        "use_roughness": bool(m.background.use_roughness),
        "saxs_enabled": bool(m.saxs.enabled),
        "saxs_mode": m.saxs.mode,
        "use_fractal": bool(m.saxs.use_fractal),
        "waxs_enabled": bool(m.waxs.enabled),
        "waxs_envelope": m.waxs.envelope,
        "use_orientation_factor": bool(m.waxs.use_orientation_factor),
        "link_R_to_pore": bool(m.waxs.link_R_to_pore),
        "link_D_to_saxs": bool(m.waxs.link_D_to_saxs),
        "link_sigma_to_saxs": bool(m.waxs.link_sigma_to_saxs),
        "formula": m.material.formula,
        "rho_struc_mode": m.material.rho_struc_mode,
        "porosity_mode": m.material.porosity_mode,
        "contrast_mode": m.material.contrast_mode,
        "peaks": _peak_summary(m),
        "material": {k: _finite(v) for k, v in material.items()},
        "n_free_parameters": len(m.fitted_refs()),
    }


def configure_carbon_model(
    session_id: str,
    saxs_mode: Optional[str] = None,
    waxs_envelope: Optional[str] = None,
    background_enabled: Optional[bool] = None,
    saxs_enabled: Optional[bool] = None,
    waxs_enabled: Optional[bool] = None,
    use_roughness: Optional[bool] = None,
    use_fractal: Optional[bool] = None,
    use_orientation_factor: Optional[bool] = None,
    link_R_to_pore: Optional[bool] = None,
    link_D_to_saxs: Optional[bool] = None,
    link_sigma_to_saxs: Optional[bool] = None,
) -> dict:
    """Set the model's shape — which terms exist and which are linked.

    Every argument is optional; omitted ones are left alone.  Changing the
    shape changes which parameter keys are active, so follow this with
    :func:`list_carbon_parameters` rather than assuming.

    Returns:
        The updated configuration, or an error naming the valid choices.
    """
    from pyirena.core.carbon_fit import CARBON_SAXS_MODES, CARBON_WAXS_ENVELOPES

    s, err = _require_carbon(session_id)
    if err:
        return err
    m = s.model

    if saxs_mode is not None:
        if saxs_mode not in CARBON_SAXS_MODES:
            return make_error(
                f"Unknown SAXS-region mode '{saxs_mode}'.",
                suggestion=f"Choose one of {sorted(CARBON_SAXS_MODES)}.",
                code="BAD_SAXS_MODE",
            )
        m.saxs.mode = saxs_mode
    if waxs_envelope is not None:
        if waxs_envelope not in CARBON_WAXS_ENVELOPES:
            return make_error(
                f"Unknown WAXS envelope '{waxs_envelope}'.",
                suggestion=f"Choose one of {sorted(CARBON_WAXS_ENVELOPES)}.",
                code="BAD_WAXS_ENVELOPE",
            )
        m.waxs.envelope = waxs_envelope

    for value, owner, attr in (
        (background_enabled, m.background, "enabled"),
        (use_roughness, m.background, "use_roughness"),
        (saxs_enabled, m.saxs, "enabled"),
        (use_fractal, m.saxs, "use_fractal"),
        (waxs_enabled, m.waxs, "enabled"),
        (use_orientation_factor, m.waxs, "use_orientation_factor"),
        (link_R_to_pore, m.waxs, "link_R_to_pore"),
        (link_D_to_saxs, m.waxs, "link_D_to_saxs"),
        (link_sigma_to_saxs, m.waxs, "link_sigma_to_saxs"),
    ):
        if value is not None:
            setattr(owner, attr, bool(value))
    return get_carbon_config(session_id)


# ---------------------------------------------------------------------------
# Parameters
# ---------------------------------------------------------------------------

def list_carbon_parameters(session_id: str, active_only: bool = True) -> dict:
    """Every parameter the model currently has, with value, bounds and Fit? flag.

    Args:
        session_id: An open session with a Carbon model.
        active_only: Only the parameters the current shape actually uses.
            Pass False to see the greyed-out ones too — useful for working out
            what a mode switch would expose.

    Returns:
        dict with ``parameters`` (a list of ``{key, label, units, value, fit,
        lo, hi}``) and ``n_free``.
    """
    s, err = _require_carbon(session_id)
    if err:
        return err
    rows = []
    for ref in s.model.parameter_refs(active_only=bool(active_only)):
        lo, hi = ref.limits
        rows.append({
            "key": ref.key, "label": ref.label, "units": ref.unit,
            "value": _finite(ref.value), "fit": bool(ref.fit),
            "lo": _finite(lo), "hi": _finite(hi),
        })
    return {"ok": True, "session_id": session_id, "parameters": rows,
            "n_free": len(s.model.fitted_refs())}


def set_carbon_parameter(session_id: str, key: str, value: float) -> dict:
    """Set one parameter's value by its dotted key."""
    s, err = _require_carbon(session_id)
    if err:
        return err
    ref = _find_ref(s.model, key)
    if ref is None:
        return _bad_key(s.model, key)
    try:
        ref.value = float(value)
    except (TypeError, ValueError):
        return make_error(f"'{value}' is not a number.",
                          suggestion="Pass a float.", code="BAD_VALUE")
    lo, hi = ref.limits
    return {"ok": True, "key": key, "value": _finite(ref.value),
            "lo": _finite(lo), "hi": _finite(hi),
            "within_bounds": bool(lo <= ref.value <= hi)}


def set_carbon_parameter_fit(session_id: str, key: str, fit: bool) -> dict:
    """Mark one parameter as refined or held fixed."""
    s, err = _require_carbon(session_id)
    if err:
        return err
    ref = _find_ref(s.model, key)
    if ref is None:
        return _bad_key(s.model, key)
    setattr(ref.owner, f"fit_{ref.attr}", bool(fit))
    return {"ok": True, "key": key, "fit": bool(fit),
            "n_free": len(s.model.fitted_refs())}


def set_carbon_parameter_bounds(session_id: str, key: str,
                                lo: Optional[float] = None,
                                hi: Optional[float] = None) -> dict:
    """Set one parameter's fitting bounds; omitted sides are left alone."""
    s, err = _require_carbon(session_id)
    if err:
        return err
    ref = _find_ref(s.model, key)
    if ref is None:
        return _bad_key(s.model, key)
    cur_lo, cur_hi = ref.limits
    new_lo = float(lo) if lo is not None else cur_lo
    new_hi = float(hi) if hi is not None else cur_hi
    if new_lo >= new_hi:
        return make_error(
            f"Lower bound {new_lo} is not below upper bound {new_hi}.",
            suggestion="Pass lo < hi.", code="BAD_BOUNDS")
    setattr(ref.owner, f"{ref.attr}_limits", (new_lo, new_hi))
    return {"ok": True, "key": key, "lo": new_lo, "hi": new_hi}


# ---------------------------------------------------------------------------
# Peaks
# ---------------------------------------------------------------------------

def list_carbon_peaks(session_id: str) -> dict:
    """Every diffraction peak with its position, d-spacing, widths and amplitude."""
    s, err = _require_carbon(session_id)
    if err:
        return err
    return {"ok": True, "session_id": session_id, "peaks": _peak_summary(s.model)}


def add_carbon_peak(session_id: str, label: str, Q0: float,
                    K: float = 1.0, FWHM_G: float = 0.3,
                    FWHM_L: float = 0.15) -> dict:
    """Add a reflection.

    Args:
        label: Miller index, e.g. ``"004"``.  The Material section looks peaks
            up by label, so ``"002"`` and ``"100"`` are the two that feed the
            density calculation.
        Q0: Peak centre [Å⁻¹]; d = 2π/Q₀.
        K: Amplitude, before the 1/Q² and Debye-Waller factors.
        FWHM_G: Gaussian width [Å⁻¹] — crystallite-size broadening.
        FWHM_L: Lorentzian width [Å⁻¹] — layer-curvature broadening.
    """
    s, err = _require_carbon(session_id)
    if err:
        return err
    if any(p.label == label for p in s.model.peaks):
        return make_error(
            f"A peak labelled '{label}' already exists.",
            suggestion="Use a different label, or remove the existing peak.",
            code="DUPLICATE_PEAK")
    peak = s.model.add_peak(label=str(label), Q0=float(Q0))
    peak.K, peak.FWHM_G, peak.FWHM_L = float(K), float(FWHM_G), float(FWHM_L)
    return {"ok": True, "peaks": _peak_summary(s.model)}


def remove_carbon_peak(session_id: str, index: int) -> dict:
    """Remove the peak at ``index`` (see :func:`list_carbon_peaks`)."""
    s, err = _require_carbon(session_id)
    if err:
        return err
    if not isinstance(index, int) or not 0 <= index < len(s.model.peaks):
        return make_error(
            f"No peak at index {index} (there are {len(s.model.peaks)}).",
            suggestion="Call list_carbon_peaks() to see the current indices.",
            code="BAD_PEAK")
    s.model.remove_peak(index)
    return {"ok": True, "peaks": _peak_summary(s.model)}


# ---------------------------------------------------------------------------
# Material
# ---------------------------------------------------------------------------

def set_carbon_material(
    session_id: str,
    formula: Optional[str] = None,
    rho_struc_mode: Optional[str] = None,
    rho_struc: Optional[float] = None,
    porosity_mode: Optional[str] = None,
    porosity: Optional[float] = None,
    contrast_mode: Optional[str] = None,
    contrast_porod: Optional[float] = None,
    contrast_micropore: Optional[float] = None,
) -> dict:
    """Set the composition → density → contrast chain, or override a stage of it.

    Each stage is normally computed: the fitted (002) and (100) positions give
    the structural density, the SAXS porosity gives the sample density, and
    each density gives an SLD and a contrast.  Override a stage when the data
    cannot supply it — no usable (100) peak, or the Teubner-Strey branch, whose
    porosity is derived from the very contrast it would feed.

    Args:
        formula: Chemical formula of the solid.
        rho_struc_mode: ``"from_peaks"`` or ``"manual"``.
        rho_struc: Structural density [g/cm³] for the manual mode.
        porosity_mode: ``"auto"`` (from the SAXS region) or ``"manual"``.
        porosity: Micropore volume fraction for the manual mode.
        contrast_mode: ``"auto"`` (from the SLDs) or ``"manual"``.
        contrast_porod: Grain-vs-vacuum (Δρ)² [10²⁰ cm⁻⁴], manual mode.
        contrast_micropore: Pore-vs-matrix (Δρ)² [10²⁰ cm⁻⁴], manual mode.

    Returns:
        The resolved chain, so the effect of the change is visible immediately.
    """
    s, err = _require_carbon(session_id)
    if err:
        return err
    mat = s.model.material

    for value, allowed, attr in (
        (rho_struc_mode, ("from_peaks", "manual"), "rho_struc_mode"),
        (porosity_mode, ("auto", "manual"), "porosity_mode"),
        (contrast_mode, ("auto", "manual"), "contrast_mode"),
    ):
        if value is None:
            continue
        if value not in allowed:
            return make_error(
                f"'{value}' is not valid for {attr}.",
                suggestion=f"Choose one of {list(allowed)}.", code="BAD_MODE")
        setattr(mat, attr, value)

    for value, attr in ((formula, "formula"),
                        (rho_struc, "rho_struc_manual"),
                        (porosity, "porosity_manual"),
                        (contrast_porod, "contrast_porod_manual"),
                        (contrast_micropore, "contrast_micropore_manual")):
        if value is None:
            continue
        setattr(mat, attr, str(value) if attr == "formula" else float(value))

    return {"ok": True,
            "material": {k: _finite(v)
                         for k, v in s.model.resolve_material().items()}}


# ---------------------------------------------------------------------------
# Fitting and results
# ---------------------------------------------------------------------------

def run_carbon_fit(session_id: str, weighting: str = "auto",
                   n_mc_runs: int = 0) -> dict:
    """Refine every ticked parameter across the whole Q range at once.

    Args:
        weighting: ``"auto"`` uses the measured uncertainties when they exist
            and relative weighting otherwise; ``"sigma"``, ``"relative"`` and
            ``"log"`` force a choice.  Relative weighting is usually right for
            this tool even when errors exist: the fit spans five decades in Q
            and ten or more in intensity, and absolute weighting lets the
            low-Q Porod region set every parameter.
        n_mc_runs: Monte-Carlo passes for the uncertainty estimate; 0 uses the
            covariance estimate alone.

    Returns:
        dict with ``success``, χ², the fitted parameters with 1-σ
        uncertainties, and the derived quantities.
    """
    s, err = _require_carbon(session_id)
    if err:
        return err

    m = s.model
    if weighting not in ("auto", "sigma", "relative", "log"):
        return make_error(
            f"Unknown weighting '{weighting}'.",
            suggestion="Choose auto, sigma, relative or log.",
            code="BAD_WEIGHTING")
    m.weighting = weighting
    m.n_mc_runs = int(n_mc_runs or 0)
    m.q_min = float(s.fit_q_min or 0.0)
    m.q_max = float(s.fit_q_max or 0.0)

    mask = fit_mask(s)
    if not np.any(mask):
        return make_error(
            "No data points in the current fit Q range.",
            suggestion="Call reset_fit_q_range() or widen set_fit_q_range().",
            code="EMPTY_RANGE")

    try:
        result = m.fit(s.q, s.intensity, s.error)
    except ValueError as exc:
        return make_error(str(exc),
                          suggestion="Tick at least one Fit? box with "
                                     "set_carbon_parameter_fit().",
                          code="NOTHING_TO_FIT")
    except Exception as exc:                            # pragma: no cover
        return make_error(f"Fit failed with exception: {exc}",
                          suggestion="Check the starting values and the Q range.",
                          code="FIT_EXCEPTION")
    finally:
        m.n_mc_runs = 0

    s.last_fit_result = result
    return get_carbon_results(session_id)


def get_carbon_results(session_id: str) -> dict:
    """The last fit's parameters, uncertainties and derived quantities.

    The derived block is the answer to most carbon questions — BET-comparable
    specific surface areas, pore and wall widths, stack height, layer count,
    lattice spacings and densities — rather than the fit coefficients that
    produced them.
    """
    s, err = _require_carbon(session_id)
    if err:
        return err
    if s.last_fit_result is None:
        return no_fit(session_id)

    result = s.last_fit_result
    params = {}
    for key, value in (result.params or {}).items():
        params[key] = {"value": _finite(value),
                       "std": _finite((result.errors or {}).get(key))}
    return {
        "ok": True,
        "session_id": session_id,
        "success": bool(result.success),
        "message": result.message,
        "chi_squared": _finite(result.chi_squared),
        "reduced_chi_squared": _finite(result.reduced_chi_squared),
        "n_points": int(result.n_points),
        "n_params": int(result.n_params),
        "saxs_mode": s.model.saxs.mode,
        "waxs_envelope": s.model.waxs.envelope,
        "parameters": params,
        "derived": {k: _finite(v) for k, v in (result.derived or {}).items()},
        # A pinned parameter is the one failure that looks like success: the
        # fit "ran", the value is reported, and it came from the bound.
        "warnings": list(result.warnings),
        "pinned_parameters": [
            {"key": key, "bound": side}
            for key, side in s.model.pinned_fitted_parameters()
        ],
        # Scalars only: fit_quality_metrics also returns per-point
        # arrays, and this layer must stay JSON-serialisable.
        "fit_quality": _quality_scalars(result.quality),
    }


def get_carbon_fit_image(session_id: str, width: int = 1100,
                         height: int = 850, dpi: int = 120) -> dict:
    """Render the fit as a PNG: data, total model, the three components, residuals.

    Seeing which component owns which decade is how you spot the failure this
    model actually has — the grain Porod term creeping up under the micropore
    region, or a diffraction peak absorbing the high-Q background.
    """
    s, err = _require_carbon(session_id)
    if err:
        return err
    if s.last_fit_result is None:
        return no_fit(session_id)

    import matplotlib  # noqa: PLC0415
    matplotlib.use("Agg", force=True)
    import matplotlib.pyplot as plt  # noqa: PLC0415

    result = s.last_fit_result
    q = np.asarray(result.q, dtype=float)
    fig, (ax1, ax2) = plt.subplots(
        2, 1, figsize=(width / dpi, height / dpi), sharex=True,
        gridspec_kw={"height_ratios": [3, 1]})

    ax1.loglog(q, np.asarray(result.I_data, dtype=float), "o", markersize=2.5,
               alpha=0.5, color="0.3", label="Data")
    ax1.loglog(q, np.asarray(result.I_model, dtype=float), "-", linewidth=2,
               color="red", label="Total fit")
    for arr, label, colour in ((result.I_porod, "Grain Porod", "tab:blue"),
                               (result.I_mp, "Micropores", "tab:orange"),
                               (result.I_waxs, "Diffraction", "tab:green")):
        curve = np.asarray(arr, dtype=float)
        good = np.isfinite(curve) & (curve > 0)
        if good.sum() > 1:
            ax1.loglog(q[good], curve[good], "--", linewidth=1.1,
                       color=colour, label=label)
    ax1.set_ylabel("I  (cm$^{-1}$)")
    ax1.legend(fontsize="small")
    ax1.grid(True, which="both", alpha=0.25)
    ax1.set_title(f"Carbon model — reduced χ² = "
                  f"{result.reduced_chi_squared:.4g}")

    ax2.semilogx(q, np.asarray(result.residuals, dtype=float), "o",
                 markersize=2.5, color="0.3")
    ax2.axhline(0.0, color="k", linewidth=1, linestyle="--")
    ax2.set_xlabel("Q  (Å$^{-1}$)")
    ax2.set_ylabel("Residuals")
    ax2.grid(True, which="both", alpha=0.25)
    fig.tight_layout()

    b64, path = render_png(fig, f"carbon_fit_{session_id}", dpi)
    return {"ok": True, "image_base64": b64, "image_path": path}


def save_carbon_fit(session_id: str, output_path: Optional[str] = None) -> dict:
    """Save the fit to NXcanSAS HDF5 under ``entry/carbon_fit_results``.

    The model's own settings go in as the embedded ``_pyirena_config``, so the
    saved file reopens in the GUI panel with every control where the agent left
    it.

    Args:
        session_id: An open session with a completed fit.
        output_path: Where to save.  Defaults to the original file, in place;
            a different path is created from the source with the previous
            results stripped.
    """
    s, err = _require_carbon(session_id)
    if err:
        return err
    if s.last_fit_result is None:
        return no_fit(session_id)

    from pyirena.io.nxcansas_carbon_fit import (  # noqa: PLC0415
        save_carbon_fit_results,
    )

    try:
        src = resolve_safe(s.file_path, must_exist=False)
        target = resolve_safe(output_path, must_exist=False) if output_path else src
    except PathSecurityError as exc:
        return make_error(str(exc),
                          suggestion="Save to a path inside PYIRENA_DATA_ROOT.",
                          code="PATH_NOT_ALLOWED")

    if target != src and not target.exists():
        from pyirena.io._nxcansas_common import (  # noqa: PLC0415
            copy_and_strip_results,
        )
        try:
            copy_and_strip_results(src, target)
        except Exception as exc:
            return make_error(
                f"Could not create output file '{target}' from source: {exc}",
                suggestion="Check the source exists and the target is writable.",
                code="SAVE_ERROR")

    try:
        save_carbon_fit_results(target, s.last_fit_result, s.model)
    except Exception as exc:
        return make_error(f"Could not save results: {exc}",
                          suggestion="Check the file is writable.",
                          code="SAVE_ERROR")
    return {"ok": True, "saved_to": str(target),
            "group": "entry/carbon_fit_results"}
