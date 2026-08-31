#!/usr/bin/env python3
"""Fit every synthetic validation dataset with pyIrena and tabulate the results.

    python3 validationData/run_validation_report.py

Reads ``ground_truth.json`` (written by ``generate_validation_data.py``), runs
the appropriate pyIrena tool on each ``.h5`` file, and writes

    VALIDATION_RESULTS.md    human-readable tables, one section per tool
    VALIDATION_RESULTS.csv   the same rows, machine readable

Both carry an **Irena** column, filled from ``irena_values.csv`` — the values
obtained by analysing the same files in the Igor Pro Irena package, keyed by
dataset and quantity.  Quantities with no entry there are left blank.  Because
those numbers live in the repository as data rather than as text inside a
document, regenerating the report never loses them; to add more, add rows to the
CSV (or edit a copy of the table and run ``fill_irena_deviations.py
--export-csv``) and run this script again.

``irena_notes.md``, if present, is appended verbatim as a commentary section.

Starting values for every fit are deliberately offset from the truth (see
``START_OFFSET``) so that convergence to the correct answer is demonstrated,
not assumed.
"""

from __future__ import annotations

import json
import sys
import traceback
from datetime import date
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(HERE))

from pyirena.io.hdf5 import readGenericNXcanSAS  # noqa: E402

#: Multiplicative offset applied to every fit starting value.  0.60 means each
#: fit starts 40 % below the true value (2.0 for a couple of shape parameters
#: where a low start is degenerate).
START_OFFSET = 0.60

ROWS: list[dict] = []
NOTES: list[str] = []

IRENA_CSV = HERE / "irena_values.csv"
IRENA_NOTES = HERE / "irena_notes.md"

#: {(dataset, quantity): (value or None, status, note)} from irena_values.csv.
#: A None value with status "not comparable" marks a quantity Irena defines
#: differently, rendered as ``*`` rather than left blank.
IRENA: dict = {}


def load_irena() -> None:
    """Read the hand-entered Irena results, if the CSV is present."""
    import csv
    if not IRENA_CSV.exists():
        return
    with IRENA_CSV.open(newline="") as fh:
        for r in csv.DictReader(fh):
            raw = (r.get("irena_value") or "").strip()
            try:
                val = float(raw) if raw else None
            except ValueError:
                NOTES.append(f"irena_values.csv: could not parse "
                             f"{r['dataset']}/{r['quantity']} = {raw!r}")
                continue
            IRENA[(r["dataset"], r["quantity"])] = (
                val, (r.get("status") or "").strip(),
                (r.get("note") or "").strip(), raw)


def irena_for(row: dict):
    """Return (display string, deviation from truth or nan) for one row.

    The value is rendered exactly as it was entered in the CSV: the significant
    figures quoted are the experimenter's claim about the precision of that
    result, and re-formatting them would overstate or lose it.  The deviation,
    by contrast, is computed against the full-precision true value rather than
    the rounded one printed in the table.
    """
    hit = IRENA.get((row["dataset"], row["quantity"]))
    if hit is None:
        return "", float("nan")
    val, status, _note, raw = hit
    if val is None:
        return ("*" if status == "not comparable" else ""), float("nan")
    truth = row["truth"]
    if isinstance(truth, (int, float)) and truth != 0:
        return raw, 100.0 * (val - truth) / truth
    return raw, float("nan")


def row(dataset, tool, quantity, unit, truth, fitted, tol_pct, settings="",
        comment=""):
    """Record one truth-vs-recovered comparison."""
    dev = float("nan")
    ok = ""
    if isinstance(truth, (int, float)) and isinstance(fitted, (int, float)):
        if truth != 0:
            dev = 100.0 * (fitted - truth) / truth
            if tol_pct is not None and np.isfinite(dev):
                ok = "PASS" if abs(dev) <= tol_pct else "FAIL"
        elif tol_pct is not None:
            # Quantities whose target is zero (a deviation, a residual offset)
            # are scored on the value itself, not on a relative deviation.
            ok = "PASS" if abs(fitted) <= tol_pct else "FAIL"
    ROWS.append(dict(dataset=dataset, tool=tool, quantity=quantity, unit=unit,
                     truth=truth, pyirena=fitted, dev_pct=dev,
                     tol_pct=tol_pct, status=ok, settings=settings,
                     comment=comment))


def curve_rows(meta, I_model, settings, tol_pct=2.0):
    """Compare the fitted curve with the exact generating curve.

    Tool-independent and easy to reproduce in Irena: export the fitted model
    intensity on the data Q grid and compare it with ``<name>_ideal.dat``.
    """
    I_ideal = ideal_of(meta)
    if I_model is None or len(I_model) != len(I_ideal):
        return
    d = 100.0 * np.abs(np.asarray(I_model, float) - I_ideal) / I_ideal
    row(meta["name"], meta["tool"], "model vs exact curve: median abs. deviation", "%",
        0.0, float(np.median(d)), tol_pct, settings,
        "fitted model compared point-by-point with <name>_ideal.dat")
    row(meta["name"], meta["tool"], "model vs exact curve: max abs. deviation", "%",
        0.0, float(np.max(d)), max(5.0 * tol_pct, 12.0), settings)


def truth_of(meta, name):
    for p in list(meta["parameters"]) + list(meta["derived"]):
        if p["name"] == name:
            return p["value"]
    raise KeyError(f"{name} not in ground truth of {meta['name']}")


def load(meta):
    """Return (q, I, dI) from the dataset's NXcanSAS file."""
    path = HERE / meta["files"]["nxcansas"]
    d = readGenericNXcanSAS(str(path.parent), path.name)
    return (np.asarray(d["Q"], float), np.asarray(d["Intensity"], float),
            np.asarray(d["Error"], float))


def ideal_of(meta):
    """The exact, noise-free curve stored with the file."""
    import h5py
    with h5py.File(HERE / meta["files"]["nxcansas"], "r") as f:
        return np.asarray(f["entry/ground_truth/I_ideal"][()], float)


# ===========================================================================
# Unified Fit
# ===========================================================================

def fit_unified(meta):
    from pyirena.core.unified import UnifiedFitModel

    q, I, dI = load(meta)
    n_levels = 2 if meta["name"] == "unified_two_level" else 1
    m = UnifiedFitModel(n_levels)
    m.background = truth_of(meta, "background") * START_OFFSET
    m.fit_background = True

    if n_levels == 1:
        specs = [("G", "Rg", "P", "B", None)]
    else:
        specs = [("G_1", "Rg_1", "P_1", "B_1", None),
                 ("G_2", "Rg_2", "P_2", "B_2", "RgCO_2")]

    for lv, (kG, kRg, kP, kB, kCO) in zip(m.levels, specs):
        lv.G = truth_of(meta, kG) * START_OFFSET
        lv.Rg = truth_of(meta, kRg) * START_OFFSET
        lv.P = truth_of(meta, kP) * 0.90
        lv.B = truth_of(meta, kB) * START_OFFSET
        lv.fit_G = lv.fit_Rg = lv.fit_B = lv.fit_P = True
        if kCO:
            lv.RgCO = truth_of(meta, kCO)
            lv.link_RGCO = True
            lv.fit_RgCO = False

    if meta.get("slit_length"):
        m.use_slit_smearing = True
        m.slit_length = float(meta["slit_length"])

    res = m.fit(q, I, dI)
    settings = (f"{n_levels} level(s), fit G/Rg/B/P + flat background, "
                f"TRF, start = {START_OFFSET:.2f}x truth"
                + (f", slit smearing dQl={meta['slit_length']:g}"
                   if meta.get("slit_length") else ""))

    tols = {"G": 3.0, "Rg": 2.0, "P": 2.0, "B": 25.0, "background": 15.0}
    for lv, (kG, kRg, kP, kB, _kCO) in zip(m.levels, specs):
        row(meta["name"], meta["tool"], kG, "1/cm", truth_of(meta, kG), lv.G,
            tols["G"], settings)
        row(meta["name"], meta["tool"], kRg, "A", truth_of(meta, kRg), lv.Rg,
            tols["Rg"], settings)
        row(meta["name"], meta["tool"], kP, "-", truth_of(meta, kP), lv.P,
            tols["P"], settings)
        row(meta["name"], meta["tool"], kB, "cm^-1 A^-P", truth_of(meta, kB),
            lv.B, tols["B"], settings,
            "B and P are strongly correlated (B carries units cm^-1 A^-P), so "
            "B has a wide tolerance whenever P is a free parameter")
    row(meta["name"], meta["tool"], "background", "1/cm",
        truth_of(meta, "background"), m.background, tols["background"], settings)
    row(meta["name"], meta["tool"], "reduced chi^2", "-", 1.0,
        float(res.get("reduced_chi_squared", float("nan"))), None, settings,
        "1.0 expected for a correct model and correct error bars")
    curve_rows(meta, m.calculate_intensity_smeared(q) if m.use_slit_smearing
               else m.calculate_intensity(q), settings)


# ===========================================================================
# Simple Fits
# ===========================================================================

SIMPLE_MAP = {
    "simple_guinier":       ("Guinier", ["I0", "Rg"], {}),
    "simple_guinier_rod":   ("Guinier Rod", ["I0", "Rc"], {}),
    "simple_guinier_sheet": ("Guinier Sheet", ["I0", "Rg"], {}),
    "simple_porod":         ("Porod", ["Kp", "Background"], {}),
    "simple_power_law":     ("Power Law", ["Prefactor", "Exponent", "Background"], {}),
    "simple_sphere":        ("Sphere", ["Scale", "R"], {"complex_bg": ["BG_flat"]}),
    "simple_debye_chain":   ("Debye Polymer Chain", ["Scale", "Rg"], {}),
    "simple_debye_bueche":  ("Debye-Bueche", ["Prefactor", "Eta", "CorrLength"], {}),
}


def fit_simple(meta):
    from pyirena.core.simple_fits import SimpleFitModel

    model_name, keys, opts = SIMPLE_MAP[meta["name"]]
    q, I, dI = load(meta)
    m = SimpleFitModel()
    m.select_model(model_name) if hasattr(m, "select_model") else None
    m.model = model_name
    m.reset_params() if hasattr(m, "reset_params") else None

    for k in keys:
        t = truth_of(meta, k)
        # An exponent must not start too far from truth or the fit is unstable
        start = t * (0.90 if k in ("Exponent",) else START_OFFSET)
        if k in ("Background", "BG_flat") and t == 0:
            start = 0.0
        m.params[k] = start

    fixed = None
    if "complex_bg" in opts:
        m.use_complex_bg = True
        for k in opts["complex_bg"]:
            m.params[k] = truth_of(meta, k) * START_OFFSET
        # Only the flat term of the complex background is present in the truth
        fixed = {"BG_B": 0.0, "BG_P": 4.0}

    res = m.fit(q, I, dI, fixed_params=fixed)
    settings = (f"model '{model_name}', all listed parameters free, "
                f"start = {START_OFFSET:.2f}x truth"
                + (", complex background (flat term only)" if "complex_bg" in opts else ""))

    if not res.get("success"):
        NOTES.append(f"{meta['name']}: Simple Fit failed - {res.get('error')}")
        return
    got = res["params"]
    tol = {"Exponent": 1.0}
    # Debye-Bueche is degenerate: only the product Prefactor*Eta^2 is determined
    degenerate = {"simple_debye_bueche": ("Prefactor", "Eta")}.get(meta["name"], ())
    all_keys = list(keys) + list(opts.get("complex_bg", []))
    for k in all_keys:
        note = ("degenerate with the other prefactor - only the product is "
                "determined by the data" if k in degenerate else "")
        row(meta["name"], meta["tool"], k,
            next(p["unit"] for p in meta["parameters"] if p["name"] == k),
            truth_of(meta, k), float(got[k]),
            None if k in degenerate else tol.get(k, 3.0), settings, note)
    if degenerate:
        t = truth_of(meta, "Prefactor") * truth_of(meta, "Eta") ** 2
        f = float(got["Prefactor"]) * float(got["Eta"]) ** 2
        row(meta["name"], meta["tool"], "Prefactor * Eta^2", "cm^-1 A^-3",
            t, f, 3.0, settings, "the combination the data actually determine")
    row(meta["name"], meta["tool"], "reduced chi^2", "-", 1.0,
        float(res["chi2"]) / max(res["dof"], 1), None, settings)
    curve_rows(meta, res.get("I_model"), settings)


# ===========================================================================
# Size Distribution
# ===========================================================================

#: Explicit tool settings for every Size Distribution dataset.  Use exactly
#: these in Irena so the two packages solve the same inverse problem.
SIZES_SETUP = {
    "sizes_sphere_lognormal": dict(
        r_min=20.0, r_max=320.0, n_bins=80, log_spacing=False,
        shape="sphere", contrast=100.0, method="maxent"),
    "sizes_sphere_bimodal": dict(
        r_min=10.0, r_max=600.0, n_bins=120, log_spacing=True,
        shape="sphere", contrast=100.0, method="maxent"),
    "sizes_sphere_flat_background": dict(
        r_min=20.0, r_max=400.0, n_bins=80, log_spacing=False,
        shape="sphere", contrast=100.0, method="maxent",
        fit_flat_bg=(0.30, 0.35)),
    "sizes_spheroid_lognormal": dict(
        r_min=30.0, r_max=450.0, n_bins=80, log_spacing=False,
        shape="spheroid", shape_params={"aspect_ratio": 0.30},
        contrast=50.0, method="maxent"),
}


def _dist_moments(r, p):
    """Volume fraction, modal radius, volume-weighted mean and RMS radius."""
    vf = float(np.trapezoid(p, r))
    mode = float(r[int(np.argmax(p))])
    mean = float(np.trapezoid(r * p, r) / vf) if vf > 0 else float("nan")
    rms = float(np.sqrt(np.trapezoid(r ** 2 * p, r) / vf)) if vf > 0 else float("nan")
    return vf, mode, mean, rms


def fit_sizes(meta):
    from pyirena.core.sizes import SizesDistribution

    cfg = SIZES_SETUP[meta["name"]]
    q, I, dI = load(meta)
    s = SizesDistribution()
    s.r_min, s.r_max = cfg["r_min"], cfg["r_max"]
    s.n_bins, s.log_spacing = cfg["n_bins"], cfg["log_spacing"]
    s.shape, s.contrast = cfg["shape"], cfg["contrast"]
    s.shape_params = dict(cfg.get("shape_params", {}))
    s.method = cfg["method"]

    bg_note = ""
    if "fit_flat_bg" in cfg:
        lo, hi = cfg["fit_flat_bg"]
        s.power_law_B = 0.0
        r_bg = s.fit_background_term(q, I, lo, hi)
        bg_note = f", flat background fitted over Q = {lo}-{hi} 1/A"
        row(meta["name"], meta["tool"], "background_flat", "1/cm",
            truth_of(meta, "background_flat"), float(s.background), 15.0,
            f"averaged over Q = {lo}-{hi} 1/A, the same recipe Irena uses; the "
            "residual particle signal there biases it slightly high")
        if not r_bg.get("success"):
            NOTES.append(f"{meta['name']}: background fit failed")

    res = s.fit(q, I, dI)
    if not res.get("success"):
        NOTES.append(f"{meta['name']}: Size Distribution fit failed - {res.get('message')}")
        return

    settings = (f"{cfg['method']}, {cfg['shape']}, r = {cfg['r_min']:g}-{cfg['r_max']:g} A "
                f"in {cfg['n_bins']} {'log' if cfg['log_spacing'] else 'linear'} bins, "
                f"contrast = {cfg['contrast']:g}e20 cm^-4" + bg_note)

    r, p = np.asarray(res["r_grid"]), np.asarray(res["distribution"])

    if meta["name"] == "sizes_sphere_bimodal":
        split = 100.0
        for tag, sel, kmed, kphi in (
            ("mode 1", r < split, "R_mode_1", "volume_fraction_1"),
            ("mode 2", r >= split, "R_mode_2", "volume_fraction_2"),
        ):
            vf, mode, mean, rms = _dist_moments(r[sel], p[sel])
            row(meta["name"], meta["tool"], f"{tag}: R at peak", "A",
                truth_of(meta, kmed), mode, 10.0, settings,
                f"population separated at r = {split:g} A")
            row(meta["name"], meta["tool"], f"{tag}: volume fraction", "-",
                truth_of(meta, kphi), vf, 10.0, settings)
        vf, mode, mean, rms = _dist_moments(r, p)
        row(meta["name"], meta["tool"], "total volume fraction", "-",
            truth_of(meta, "volume_fraction_total"), vf, 6.0, settings)
    else:
        vf, mode, mean, rms = _dist_moments(r, p)
        row(meta["name"], meta["tool"], "volume fraction", "-",
            truth_of(meta, "volume_fraction"), vf, 6.0, settings)
        row(meta["name"], meta["tool"], "R at distribution peak", "A",
            truth_of(meta, "R_mode"), mode, 10.0, settings,
            "regularised inversions broaden the distribution, so the modal "
            "radius has a wider tolerance than the integral moments")
        if any(p_["name"] == "R_mean_vol" for p_ in meta["derived"]):
            row(meta["name"], meta["tool"], "volume-weighted mean R", "A",
                truth_of(meta, "R_mean_vol"), mean, 6.0, settings)
        row(meta["name"], meta["tool"], "RMS radius (pyIrena 'Rg')", "A",
            truth_of(meta, "R_rms"), rms, 6.0, settings)

    row(meta["name"], meta["tool"], "reduced chi^2", "-", 1.0,
        float(res["chi_squared"]) / len(q), None, settings)
    I_model = np.asarray(res["model_intensity"]) + float(s.background)
    curve_rows(meta, I_model, settings, tol_pct=3.0)


# ===========================================================================
# WAXS Peak Fit
# ===========================================================================

WAXS_SETUP = {
    "waxs_three_gauss": dict(bg_shape="Linear", shape="Gauss", n_peaks=3),
    "waxs_pseudovoigt": dict(bg_shape="Cubic", shape="Pseudo-Voigt", n_peaks=2),
}


def fit_waxs(meta):
    from pyirena.core import waxs_peakfit as W

    cfg = WAXS_SETUP[meta["name"]]
    q, I, dI = load(meta)

    bg = W.default_bg_params(cfg["bg_shape"])
    n_bg = len(bg)
    for i, name in enumerate(W.bg_param_names(cfg["bg_shape"])):
        bg[name]["value"] = truth_of(meta, f"bg{i}") * START_OFFSET
        bg[name]["fit"] = True

    peaks = []
    for i in range(1, cfg["n_peaks"] + 1):
        pk = W.default_peak_params(
            cfg["shape"],
            Q0=truth_of(meta, f"peak{i}_Q0") * 1.01,       # 1 % off, not 40 %
            A=truth_of(meta, f"peak{i}_A") * START_OFFSET,
            FWHM=truth_of(meta, f"peak{i}_FWHM") * 1.30)
        if cfg["shape"] == "Pseudo-Voigt":
            pk["eta"]["value"] = 0.50
        peaks.append(pk)

    m = W.WAXSPeakFitModel(cfg["bg_shape"], peaks)
    res = m.fit(q, I, dI, bg, peaks)
    settings = (f"{cfg['n_peaks']} x {cfg['shape']} peaks, {cfg['bg_shape']} background, "
                f"all parameters free; A and background start at "
                f"{START_OFFSET:.2f}x truth, Q0 within 1 %, FWHM 1.30x truth")
    if not res.get("success"):
        NOTES.append(f"{meta['name']}: WAXS fit failed - {res.get('message')}")
        return

    bg_names = W.bg_param_names(cfg["bg_shape"])
    poly_note = ("polynomial background coefficients are strongly correlated and "
                 "individually ill-determined; compare the background CURVE below")
    for i, name in enumerate(bg_names):
        row(meta["name"], meta["tool"], f"bg{i}", "-",
            truth_of(meta, f"bg{i}"), float(res["bg_params"][name]["value"]),
            10.0 if cfg["bg_shape"] == "Linear" else None, settings, poly_note)
    coef_true = [truth_of(meta, f"bg{i}") for i in range(len(bg_names))]
    coef_fit = [float(res["bg_params"][n]["value"]) for n in bg_names]
    for q_probe in np.linspace(q[0], q[-1], 5):
        bt = sum(c * q_probe ** k for k, c in enumerate(coef_true))
        bf = sum(c * q_probe ** k for k, c in enumerate(coef_fit))
        row(meta["name"], meta["tool"], f"background at Q = {q_probe:.3g}", "1/cm",
            bt, bf, 5.0, settings, "background curve, the quantity that is determined")
    for i, pk in enumerate(res["peaks"], 1):
        for key, tol in (("A", 3.0), ("Q0", 0.5), ("FWHM", 3.0), ("eta", 15.0)):
            if key not in pk:
                continue
            row(meta["name"], meta["tool"], f"peak{i}_{key}",
                {"A": "1/cm", "Q0": "1/A", "FWHM": "1/A", "eta": "-"}[key],
                truth_of(meta, f"peak{i}_{key}"), float(pk[key]["value"]),
                tol, settings)
        area = W.peak_area(pk["shape"], pk)
        row(meta["name"], meta["tool"], f"peak{i}_area", "cm^-1 A^-1",
            truth_of(meta, f"peak{i}_area"), float(area), 4.0, settings,
            "analytic integral of the fitted peak")
    row(meta["name"], meta["tool"], "reduced chi^2", "-", 1.0,
        float(res["reduced_chi2"]), None, settings)
    curve_rows(meta, res.get("I_model"), settings)
    _ = n_bg


# ===========================================================================
# Modeling
# ===========================================================================

def _sd_pop(meta, mean_key, sd_key, scale_key, contrast_key, n_bins=200):
    from pyirena.core.modeling import SizeDistPopulation
    p = SizeDistPopulation()
    p.dist_type = "lognormal"
    p.dist_params = {"min_size": 0.0,
                     "mean_size": truth_of(meta, mean_key) * START_OFFSET,
                     "sdeviation": truth_of(meta, sd_key) * START_OFFSET}
    p.dist_params_fit = {"min_size": False, "mean_size": True, "sdeviation": True}
    p.dist_params_limits = {"min_size": (0.0, 1e6), "mean_size": (1.0, 1e6),
                            "sdeviation": (0.01, 3.0)}
    p.form_factor = "sphere"
    p.contrast = truth_of(meta, contrast_key)
    p.fit_contrast = False
    p.scale = truth_of(meta, scale_key) * START_OFFSET
    p.fit_scale = True
    p.n_bins = n_bins
    return p


def _build_modeling(meta):
    """Return (config, [(quantity, unit, truth_key, getter, tol)]) for a dataset."""
    from pyirena.core import modeling as MD

    name = meta["name"]
    cfg = MD.ModelingConfig()
    cfg.q_min, cfg.q_max = meta["q_min"] * 0.99, meta["q_max"] * 1.01
    cfg.background = truth_of(meta, "background") * START_OFFSET
    cfg.fit_background = True
    checks = []

    def sd_checks(i, mean_key, sd_key, scale_key):
        return [
            ("pop1 mean_size (log-normal median)", "A", mean_key,
             lambda c: c.populations[i].dist_params["mean_size"], 4.0),
            ("pop1 sdeviation", "-", sd_key,
             lambda c: c.populations[i].dist_params["sdeviation"], 6.0),
            ("pop1 scale", "-", scale_key,
             lambda c: c.populations[i].scale, 6.0),
        ]

    if name == "modeling_sizedist_plus_unified":
        p1 = _sd_pop(meta, "pop1_mean_size", "pop1_sdeviation",
                     "pop1_scale", "pop1_contrast")
        p2 = MD.UnifiedLevelPopulation(
            G=truth_of(meta, "pop2_G") * START_OFFSET,
            Rg=truth_of(meta, "pop2_Rg") * START_OFFSET,
            P=truth_of(meta, "pop2_P") * 0.9,
            B=truth_of(meta, "pop2_B") * START_OFFSET)
        p2.fit_G = p2.fit_Rg = p2.fit_B = p2.fit_P = True
        cfg.populations = [p1, p2]
        checks = sd_checks(0, "pop1_mean_size", "pop1_sdeviation", "pop1_scale") + [
            ("pop2 G", "1/cm", "pop2_G", lambda c: c.populations[1].G, 4.0),
            ("pop2 Rg", "A", "pop2_Rg", lambda c: c.populations[1].Rg, 4.0),
            ("pop2 P", "-", "pop2_P", lambda c: c.populations[1].P, 3.0),
            ("pop2 B", "cm^-1 A^-P", "pop2_B", lambda c: c.populations[1].B, 25.0),
        ]

    elif name == "modeling_sizedist_plus_peak":
        p1 = _sd_pop(meta, "pop1_mean_size", "pop1_sdeviation",
                     "pop1_scale", "pop1_contrast")
        p2 = MD.DiffractionPeakPopulation(
            peak_type="gaussian",
            position=truth_of(meta, "pop2_position") * 1.05,
            amplitude=truth_of(meta, "pop2_amplitude") * START_OFFSET,
            width=truth_of(meta, "pop2_width") * 1.3)
        p2.fit_position = p2.fit_amplitude = p2.fit_width = True
        cfg.populations = [p1, p2]
        checks = sd_checks(0, "pop1_mean_size", "pop1_sdeviation", "pop1_scale") + [
            ("pop2 amplitude", "1/cm", "pop2_amplitude",
             lambda c: c.populations[1].amplitude, 5.0),
            ("pop2 position", "1/A", "pop2_position",
             lambda c: c.populations[1].position, 1.0),
            ("pop2 width (sigma)", "1/A", "pop2_width",
             lambda c: c.populations[1].width, 6.0),
        ]

    elif name == "modeling_hardsphere":
        p1 = _sd_pop(meta, "pop1_mean_size", "pop1_sdeviation",
                     "pop1_scale", "pop1_contrast")
        p1.structure_factor = "hard_sphere"
        p1.sf_params = dict(p1.sf_params)
        p1.sf_params["radius"] = truth_of(meta, "sf_radius") * 0.8
        p1.sf_params["volume_fraction"] = truth_of(meta, "sf_volume_fraction") * 0.6
        p1.sf_params_fit = dict(p1.sf_params_fit)
        p1.sf_params_fit["radius"] = True
        p1.sf_params_fit["volume_fraction"] = True
        cfg.populations = [p1]
        checks = sd_checks(0, "pop1_mean_size", "pop1_sdeviation", "pop1_scale") + [
            ("hard-sphere radius", "A", "sf_radius",
             lambda c: c.populations[0].sf_params["radius"], 6.0),
            ("hard-sphere volume fraction", "-", "sf_volume_fraction",
             lambda c: c.populations[0].sf_params["volume_fraction"], 10.0),
        ]

    elif name == "modeling_mass_fractal":
        p1 = MD.MassFractalPopulation(
            Phi=truth_of(meta, "Phi") * START_OFFSET,
            Radius=truth_of(meta, "Radius") * START_OFFSET,
            Dv=truth_of(meta, "Dv") * 0.9,
            Ksi=truth_of(meta, "Ksi") * START_OFFSET,
            Eta=truth_of(meta, "Eta"), Contrast=truth_of(meta, "Contrast"))
        p1.fit_Phi = p1.fit_Radius = p1.fit_Dv = p1.fit_Ksi = True
        p1.fit_Eta = p1.fit_Contrast = False
        cfg.populations = [p1]
        checks = [
            ("Phi", "-", "Phi", lambda c: c.populations[0].Phi, 6.0),
            ("Radius", "A", "Radius", lambda c: c.populations[0].Radius, 4.0),
            ("Dv", "-", "Dv", lambda c: c.populations[0].Dv, 3.0),
            ("Ksi", "A", "Ksi", lambda c: c.populations[0].Ksi, 8.0),
        ]

    elif name == "modeling_surface_fractal":
        p1 = MD.SurfaceFractalPopulation(
            Surface=truth_of(meta, "Surface") * START_OFFSET,
            Ds=truth_of(meta, "Ds") * 0.95,
            Ksi=truth_of(meta, "Ksi") * START_OFFSET,
            Contrast=truth_of(meta, "Contrast"))
        p1.fit_Surface = p1.fit_Ds = p1.fit_Ksi = True
        p1.fit_Contrast = False
        cfg.populations = [p1]
        checks = [
            ("Surface", "1/cm", "Surface", lambda c: c.populations[0].Surface, 6.0),
            ("Ds", "-", "Ds", lambda c: c.populations[0].Ds, 3.0),
            ("Ksi", "A", "Ksi", lambda c: c.populations[0].Ksi, 8.0),
        ]

    elif name == "modeling_guinier_porod":
        p1 = MD.GuinierPorodPopulation(
            G=truth_of(meta, "G") * START_OFFSET,
            Rg1=truth_of(meta, "Rg1") * START_OFFSET,
            s1=truth_of(meta, "s1") * 0.7,
            P=truth_of(meta, "P") * 0.9)
        p1.fit_G = p1.fit_Rg1 = p1.fit_s1 = p1.fit_P = True
        cfg.populations = [p1]
        checks = [
            ("G", "1/cm", "G", lambda c: c.populations[0].G, 5.0),
            ("Rg1", "A", "Rg1", lambda c: c.populations[0].Rg1, 4.0),
            ("s1", "-", "s1", lambda c: c.populations[0].s1, 5.0),
            ("P", "-", "P", lambda c: c.populations[0].P, 3.0),
        ]
    else:
        raise KeyError(name)

    return cfg, checks


def fit_modeling(meta):
    from pyirena.core.modeling import ModelingEngine

    q, I, dI = load(meta)
    cfg, checks = _build_modeling(meta)
    eng = ModelingEngine()
    res = eng.fit(cfg, q, I, dI)
    fitted = res.config

    def _pop_label(p):
        t = getattr(p, "pop_type", "size_dist")
        c = getattr(p, "contrast", None)
        if c is None:
            c = getattr(p, "Contrast", None)
        return f"{t} (contrast {c:g}e20 cm^-4)" if c is not None else t

    pops = ", ".join(_pop_label(p) for p in cfg.populations)
    settings = (f"populations: {pops}; local (TRF) fit, "
                f"size-distribution grid 200 bins, contrast held fixed at the "
                f"value quoted above, start = {START_OFFSET:.2f}x truth")

    for label, unit, key, getter, tol in checks:
        row(meta["name"], meta["tool"], label, unit, truth_of(meta, key),
            float(getter(fitted)), tol, settings)
    row(meta["name"], meta["tool"], "background", "1/cm",
        truth_of(meta, "background"), float(fitted.background), 30.0, settings,
        "a small flat background is a weakly determined nuisance parameter")
    row(meta["name"], meta["tool"], "reduced chi^2", "-", 1.0,
        float(res.reduced_chi_squared), None, settings)
    curve_rows(meta, res.model_I, settings, tol_pct=3.0)


# ===========================================================================
# Data Merge
# ===========================================================================

MERGE_OVERLAP = (0.0060, 0.0400)


def fit_merge(metas):
    from pyirena.core.data_merge import DataMerge, MergeConfig

    m1, m2 = metas["merge_usaxs"], metas["merge_saxs"]
    q1, I1, dI1 = load(m1)
    q2, I2, dI2 = load(m2)
    lo, hi = MERGE_OVERLAP
    cfg = MergeConfig(q_overlap_min=lo, q_overlap_max=hi,
                      fit_scale=True, scale_dataset=2, fit_qshift=False)
    res = DataMerge().optimize(q1, I1, dI1, q2, I2, dI2, cfg)
    settings = (f"DS1 = merge_usaxs (absolute), DS2 = merge_saxs (scaled); "
                f"overlap Q = {lo}-{hi} 1/A, log-log interpolation, scale DS2, "
                f"no Q shift")
    row("merge_usaxs + merge_saxs", "Data Merge", "recovered scale for DS2", "-",
        truth_of(m2, "expected_recovered_scale"), float(res.scale), 1.0, settings,
        "DS2 was generated 1.15x too high, so the correct scale is 1/1.15")
    row("merge_usaxs + merge_saxs", "Data Merge", "fitted background offset",
        "1/cm", 0.0, float(res.background), None, settings,
        "should be ~0: the two branches differ by a pure multiplicative factor")


# ===========================================================================
# Data Manipulation
# ===========================================================================

def fit_manipulation(metas):
    from pyirena.core.data_manipulation import DataManipulation, ScaleConfig

    m_in, m_ref = metas["manipulate_input"], metas["manipulate_reference"]
    q, I, dI = load(m_in)
    a = truth_of(m_in, "a_multiply")
    b = truth_of(m_in, "b_add")
    cfg = ScaleConfig(scale_I=1.0 / a, background=b / a)
    res = DataManipulation.scale(q, I, dI, None, cfg)

    _, I_ref, _ = load(m_ref)
    d = 100.0 * np.abs(res.I - I_ref) / np.abs(I_ref)
    settings = (f"Scale + Background: I_out = (1/{a:g})*I - {b/a:.6g}, "
                f"the exact inverse of the generating operation")
    row("manipulate_input", "Data Manipulation",
        "max abs. deviation from manipulate_reference", "%", 0.0, float(np.max(d)),
        1e-8, settings,
        "the inverse operation must reproduce the reference to machine precision")


# ===========================================================================
# Scattering Contrast  (no data file: a calculator cross-check)
# ===========================================================================

#: Standard atomic weights (IUPAC 2021) and atomic numbers, plus bound coherent
#: neutron scattering lengths b_c in fm from Sears, Neutron News 3 (1992) 26.
#: Used to compute the expected SLDs independently of pyIrena.
_ELEMENTS = {
    "H":  (1.008,   1,  -3.7390),
    "D":  (2.014,   1,   6.6710),
    "C":  (12.011,  6,   6.6460),
    "O":  (15.999,  8,   5.8030),
    "Al": (26.982, 13,   3.4490),
    "Si": (28.085, 14,   4.1491),
    "Fe": (55.845, 26,   9.4500),
}

CONTRAST_CASES = [
    ("SiO2", {"Si": 1, "O": 2}, 2.200, "amorphous silica"),
    ("H2O", {"H": 2, "O": 1}, 1.000, "light water"),
    ("D2O", {"D": 2, "O": 1}, 1.107, "heavy water"),
    ("Al2O3", {"Al": 2, "O": 3}, 3.970, "corundum"),
    ("Fe", {"Fe": 1}, 7.874, "alpha iron"),
    ("C8H8", {"C": 8, "H": 8}, 1.050, "polystyrene"),
]

_R_E = 2.8179403227e-13     # classical electron radius [cm]
_N_A = 6.02214076e23


def _expected_slds(counts, density):
    """X-ray (free electron) and neutron SLD in 1e10 cm^-2, from first principles.

        rho_x = (density N_A / M) Z_tot r_e
        rho_n = (density N_A / M) SUM n_i b_i
    """
    M = sum(n * _ELEMENTS[e][0] for e, n in counts.items())
    Z = sum(n * _ELEMENTS[e][1] for e, n in counts.items())
    b = sum(n * _ELEMENTS[e][2] for e, n in counts.items()) * 1e-13   # fm -> cm
    n_per_cm3 = density * _N_A / M
    return n_per_cm3 * Z * _R_E / 1e10, n_per_cm3 * b / 1e10, M


def run_contrast():
    try:
        from pyirena.core.scattering_contrast import compute_compound, compute_contrast
    except Exception as exc:                                    # pragma: no cover
        NOTES.append(f"Scattering Contrast skipped: {exc}")
        return
    settings = ("free-electron X-ray SLD and bound-coherent neutron SLD; "
                "reference values from IUPAC 2021 atomic weights and "
                "Sears (1992) scattering lengths")
    comps = {}
    for formula, counts, density, label in CONTRAST_CASES:
        try:
            c = compute_compound(formula, density, name=label)
        except Exception as exc:
            NOTES.append(f"Scattering Contrast '{formula}' failed: {exc}")
            return
        comps[formula] = c
        x_exp, n_exp, M_exp = _expected_slds(counts, density)
        row(f"contrast: {formula} ({density:g} g/cm3)", "Scattering Contrast",
            "formula weight", "g/mol", M_exp, float(c.mol_weight), 0.2, settings, label)
        row(f"contrast: {formula} ({density:g} g/cm3)", "Scattering Contrast",
            "X-ray SLD", "1e10 cm^-2", x_exp, float(c.xray_sld), 0.5, settings)
        row(f"contrast: {formula} ({density:g} g/cm3)", "Scattering Contrast",
            "neutron SLD", "1e10 cm^-2", n_exp, float(c.neutron_sld), 1.0, settings)

    for a, b in (("SiO2", "H2O"), ("SiO2", "D2O"), ("Fe", "H2O")):
        r = compute_contrast(comps[a], comps[b])
        xa, na, _ = _expected_slds(dict(CONTRAST_CASES[[c[0] for c in CONTRAST_CASES].index(a)][1]),
                                   CONTRAST_CASES[[c[0] for c in CONTRAST_CASES].index(a)][2])
        xb, nb, _ = _expected_slds(dict(CONTRAST_CASES[[c[0] for c in CONTRAST_CASES].index(b)][1]),
                                   CONTRAST_CASES[[c[0] for c in CONTRAST_CASES].index(b)][2])
        row(f"contrast: {a} vs {b}", "Scattering Contrast",
            "X-ray contrast (delta-rho)^2", "1e20 cm^-4",
            (xa - xb) ** 2, float(r.xray_contrast), 1.0, settings)
        row(f"contrast: {a} vs {b}", "Scattering Contrast",
            "neutron contrast (delta-rho)^2", "1e20 cm^-4",
            (na - nb) ** 2, float(r.neutron_contrast), 2.0, settings)


# ===========================================================================
# Driver and output
# ===========================================================================

DISPATCH = {
    "Size Distribution": fit_sizes,
    "Unified Fit": fit_unified,
    "Unified Fit (slit smearing)": fit_unified,
    "Modeling": fit_modeling,
    "Simple Fits": fit_simple,
    "WAXS Peak Fit": fit_waxs,
}

TOOL_ORDER = ["Size Distribution", "Unified Fit", "Unified Fit (slit smearing)",
              "Modeling", "Simple Fits", "WAXS Peak Fit", "Data Merge",
              "Data Manipulation", "Scattering Contrast"]


def fmt(v):
    if v is None:
        return ""
    if isinstance(v, str):
        return v
    if not np.isfinite(v):
        return ""
    if v == 0:
        return "0"
    a = abs(v)
    return f"{v:.6g}" if 1e-3 <= a < 1e6 else f"{v:.5e}"


def write_csv(path):
    import csv
    cols = ["tool", "dataset", "quantity", "unit", "true_value",
            "pyirena_value", "pyirena_dev_pct", "tolerance_pct", "pyirena_status",
            "irena_value", "irena_dev_pct", "pyirena_minus_irena_pct",
            "pyirena_settings", "comment"]
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(cols)
        for r in ROWS:
            ir_s, ir_dev = irena_for(r)
            hit = IRENA.get((r["dataset"], r["quantity"]))
            gap = ""
            if hit and hit[0]:
                gap = fmt(100.0 * (r["pyirena"] - hit[0]) / hit[0])
            w.writerow([r["tool"], r["dataset"], r["quantity"], r["unit"],
                        fmt(r["truth"]), fmt(r["pyirena"]), fmt(r["dev_pct"]),
                        fmt(r["tol_pct"]), r["status"], ir_s, fmt(ir_dev), gap,
                        r["settings"], r["comment"]])


HEADER = """\
# pyIrena validation results

Every row compares a parameter recovered by **pyIrena** with the **exact value
used to synthesise the data**.  The data files, and the generator that produced
them, are in this folder; see `README.md` for the full ground truth and
`ground_truth.json` for the machine-readable version.

The **Irena** column holds the corresponding value obtained with the Igor Pro
Irena package, analysing the same file with the settings quoted for each
dataset.  Where it is filled, the table is a three-way comparison: both
implementations measured against one common, exactly known reference.  A blank
Irena cell means that quantity has not been analysed in Irena or that Irena does
not report it; `*` marks a quantity Irena defines differently, so the numbers
are not directly comparable.  These values are read from `irena_values.csv`.

Both deviations are `100 (fitted - true) / true`.  `Tol` is the tolerance the pyIrena
regression test enforces; a blank tolerance marks a quantity that is reported
for information but is not expected to be individually determined by the data
(a correlated prefactor, a polynomial coefficient, a reduced chi-squared).

Starting values for every fit were offset from the truth (typically 0.6x) so
that the tables demonstrate convergence rather than assume it.

"""


def write_markdown(path, failures):
    lines = [HEADER]
    lines.append(f"Generated {date.today().isoformat()} by "
                 f"`validationData/run_validation_report.py`.\n")

    n_scored = sum(1 for r in ROWS if r["status"])
    n_pass = sum(1 for r in ROWS if r["status"] == "PASS")
    lines.append(f"**Summary: {n_pass} of {n_scored} scored comparisons within "
                 f"tolerance** ({len(ROWS)} rows in total).\n")

    stats = agreement_stats()
    if stats:
        n, py_med, ir_med, gap_med = stats
        lines.append(
            f"Over the {n} quantities both packages report, the median deviation "
            f"from the known truth is **{py_med:.2f} % for pyIrena and "
            f"{ir_med:.2f} % for Irena**, and the median difference **between the "
            f"two packages is {gap_med:.2f} %**.\n")
    if NOTES:
        lines.append("Notes from this run:\n")
        lines += [f"* {n}" for n in NOTES] + [""]

    for tool in TOOL_ORDER:
        rows = [r for r in ROWS if r["tool"] == tool]
        if not rows:
            continue
        lines.append(f"## {tool}\n")
        for ds in dict.fromkeys(r["dataset"] for r in rows):
            drows = [r for r in rows if r["dataset"] == ds]
            lines.append(f"### `{ds}`\n")
            settings = drows[0]["settings"]
            if settings:
                lines.append(f"*pyIrena settings — use the same in Irena:* {settings}\n")
            lines.append("| quantity | unit | true | pyIrena | dev % | tol % | "
                         "in tol? | **Irena** | **Irena dev %** |")
            lines.append("|---|---|---|---|---|---|---|---|---|")
            for r in drows:
                mark = {"PASS": "yes", "FAIL": "**CHECK**"}.get(r["status"], "-")
                d = r["dev_pct"]
                dev_s = f"{d:+.3f}" if isinstance(d, float) and np.isfinite(d) else ""
                ir_s, ir_dev = irena_for(r)
                ir_dev_s = f"{ir_dev:+.3f}" if np.isfinite(ir_dev) else ""
                lines.append(
                    f"| {r['quantity']} | {r['unit']} | {fmt(r['truth'])} | "
                    f"{fmt(r['pyirena'])} | {dev_s} | "
                    f"{fmt(r['tol_pct'])} | {mark} | {ir_s} | {ir_dev_s} |")
            comments = dict.fromkeys(r["comment"] for r in drows if r["comment"])
            for c in comments:
                lines.append(f"\n> {c}")
            lines.append("")
    if IRENA_NOTES.exists():
        lines.append("---\n")
        lines.append(IRENA_NOTES.read_text().split("-->", 1)[-1].strip())
        lines.append("")
    if failures:
        lines.append("## Rows outside tolerance\n")
        for r in failures:
            lines.append(f"* `{r['dataset']}` — {r['quantity']}: "
                         f"true {fmt(r['truth'])}, pyIrena {fmt(r['pyirena'])} "
                         f"({fmt(r['dev_pct'])} %, tolerance {fmt(r['tol_pct'])} %)")
        lines.append("")
    Path(path).write_text("\n".join(lines) + "\n")


def agreement_stats():
    """(n, median |pyIrena-true|, median |Irena-true|, median |pyIrena-Irena|) in %."""
    py, ir, gap = [], [], []
    for r in ROWS:
        hit = IRENA.get((r["dataset"], r["quantity"]))
        if not hit or hit[0] is None:
            continue
        t, p, i = r["truth"], r["pyirena"], hit[0]
        if not isinstance(t, (int, float)) or t == 0 or not isinstance(p, (int, float)):
            continue
        py.append(abs(100 * (p - t) / t))
        ir.append(abs(100 * (i - t) / t))
        gap.append(abs(100 * (p - i) / i) if i else float("nan"))
    if not py:
        return None
    med = lambda xs: sorted(xs)[len(xs) // 2]                     # noqa: E731
    return len(py), med(py), med(ir), med(gap)


def main():
    load_irena()
    manifest = json.loads((HERE / "ground_truth.json").read_text())
    metas = {m["name"]: m for m in manifest}

    for m in manifest:
        fn = DISPATCH.get(m["tool"])
        if fn is None:
            continue
        try:
            fn(m)
        except Exception as exc:                                # pragma: no cover
            NOTES.append(f"{m['name']}: {type(exc).__name__}: {exc}")
            traceback.print_exc()
        print(f"  done: {m['name']}")

    for label, fn in (("merge", fit_merge), ("manipulation", fit_manipulation)):
        try:
            fn(metas)
        except Exception as exc:                                # pragma: no cover
            NOTES.append(f"{label}: {type(exc).__name__}: {exc}")
            traceback.print_exc()
        print(f"  done: {label}")
    try:
        run_contrast()
    except Exception as exc:                                    # pragma: no cover
        NOTES.append(f"contrast: {type(exc).__name__}: {exc}")
    print("  done: scattering contrast")

    failures = [r for r in ROWS if r["status"] == "FAIL"]
    write_csv(HERE / "VALIDATION_RESULTS.csv")
    write_markdown(HERE / "VALIDATION_RESULTS.md", failures)

    n_scored = sum(1 for r in ROWS if r["status"])
    print(f"\n{len(ROWS)} comparisons, {n_scored} scored, "
          f"{n_scored - len(failures)} within tolerance.")
    stats = agreement_stats()
    if stats:
        n, py_med, ir_med, gap_med = stats
        print(f"{n} quantities also analysed in Irena "
              f"(from {IRENA_CSV.name}): median |pyIrena-true| {py_med:.3f} %, "
              f"|Irena-true| {ir_med:.3f} %, |pyIrena-Irena| {gap_med:.3f} %")
    else:
        print(f"no Irena values found ({IRENA_CSV.name} missing or empty)")
    for r in failures:
        print(f"  OUTSIDE TOLERANCE: {r['dataset']} / {r['quantity']}: "
              f"{fmt(r['dev_pct'])} % (tol {fmt(r['tol_pct'])} %)")
    for n in NOTES:
        print(f"  NOTE: {n}")
    print("\nVALIDATION_RESULTS.md and VALIDATION_RESULTS.csv written.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
