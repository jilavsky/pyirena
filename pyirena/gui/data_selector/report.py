"""
pyirena.gui.data_selector.report — Markdown fit-report builder (_build_report).

Split from the original monolithic data_selector.py (no behavior change).
"""
import logging

log = logging.getLogger(__name__)


import os
from typing import Optional

import numpy as np




def _quality_report_rows(fq: Optional[dict]) -> list:
    """Markdown table rows for robust fit-quality metrics (empty if unavailable).

    Intended to be appended to a tool's '## Fit Quality' table in the report.
    """
    if not fq:
        return []
    rows = []
    s = fq.get('robust_scale_s')
    if s is not None:
        rows.append(f"| σ-scale (robust) | {s:.2f}× |")
        floor = fq.get('realistic_reduced_chi2_floor')
        if floor is not None:
            rows.append(f"| Realistic reduced-χ² floor | {floor:.2f} |")
    mfm = fq.get('max_abs_frac_misfit')
    if mfm is not None:
        rows.append(f"| Max \\|(I−M)/I\\| | {mfm * 100:.1f}% |")
    csr = fq.get('longest_same_sign_run')
    if csr is not None:
        rows.append(f"| Longest same-sign run | {csr} |")
    return rows


def _build_report(file_path: str,
                  data_info: Optional[dict] = None,
                  fit_results: Optional[dict] = None,
                  sizes_results: Optional[dict] = None,
                  simple_fit_results: Optional[dict] = None,
                  waxs_peakfit_results: Optional[dict] = None,
                  modeling_results: Optional[dict] = None,
                  saxs_morph_results: Optional[dict] = None) -> str:
    """
    Build a Markdown report string.

    Args:
        file_path:          Absolute path to the source file.
        data_info:          Dict with keys 'Q', 'I', 'I_error' (optional array).
                            Pass None to omit the data section.
        fit_results:        Dict from load_unified_fit_results().
                            Pass None to omit the unified fit section.
        sizes_results:      Dict from load_sizes_results().
                            Pass None to omit the size distribution section.
        simple_fit_results: Dict from load_simple_fit_results().
                            Pass None to omit the simple fits section.

    Returns:
        Multi-line Markdown string ready to be written to a .md file.
    """
    from datetime import datetime as _dt

    filename = os.path.basename(file_path)
    now = _dt.now().strftime('%Y-%m-%d %H:%M:%S')

    L = []

    # ── Header ───────────────────────────────────────────────────────────────
    L += [
        "# pyIrena Report",
        "",
        "| | |",
        "|---|---|",
        f"| **File** | `{filename}` |",
        f"| **Report generated** | {now} |",
    ]
    if fit_results is not None:
        L += [
            f"| **Unified Fit timestamp** | {fit_results.get('timestamp', 'unknown')} |",
        ]
    if sizes_results is not None:
        L += [
            f"| **Size Dist. timestamp** | {sizes_results.get('timestamp', 'unknown')} |",
        ]
    if simple_fit_results is not None:
        L += [
            f"| **Simple Fits timestamp** | {simple_fit_results.get('timestamp', 'unknown')} |",
        ]
    if waxs_peakfit_results is not None:
        L += [
            f"| **WAXS Peak Fit timestamp** | {waxs_peakfit_results.get('timestamp', 'unknown')} |",
        ]
    if saxs_morph_results is not None:
        L += [
            f"| **SAXS Morph timestamp** | {saxs_morph_results.get('timestamp', 'unknown')} |",
        ]
    if modeling_results is not None:
        L += [
            f"| **Modeling timestamp** | {modeling_results.get('timestamp', 'unknown')} |",
        ]
    L.append("")

    # ── Data summary ─────────────────────────────────────────────────────────
    from pyirena.gui.fmt_utils import eng_fmt as _ef
    if data_info is not None:
        Q = data_info['Q']
        I = data_info['I']
        I_error = data_info.get('I_error')
        L += [
            "## Data Summary",
            "",
            "| Parameter | Value |",
            "|-----------|-------|",
            f"| Q range | {_ef(Q.min())} – {_ef(Q.max())} Å⁻¹ |",
            f"| Intensity range | {_ef(I.min())} – {_ef(I.max())} cm⁻¹ |",
            f"| Data points | {len(Q)} |",
        ]
        if I_error is not None:
            L.append(
                f"| Uncertainty range | {_ef(I_error.min())} – {_ef(I_error.max())} cm⁻¹ |"
            )
        L.append("")

    # ── Fit quality ──────────────────────────────────────────────────────────
    if fit_results is not None:
        chi2      = fit_results['chi_squared']
        bg        = fit_results['background']
        n_levels  = fit_results['num_levels']
        Q_fit     = fit_results['Q']
        residuals = fit_results['residuals']
        levels    = fit_results.get('levels', [])

        L += [
            "## Fit Quality",
            "",
            "| Parameter | Value |",
            "|-----------|-------|",
            f"| Chi-squared (χ²) | {chi2:.4f} |",
            f"| Number of levels | {n_levels} |",
            f"| Background | {_ef(bg)} cm⁻¹ |",
            f"| Q range (fit) | {_ef(Q_fit.min())} – {_ef(Q_fit.max())} Å⁻¹ |",
            f"| Data points (fit) | {len(Q_fit)} |",
            f"| Residuals mean | {np.mean(residuals):.4f} |",
            f"| Residuals std dev | {np.std(residuals):.4f} |",
            f"| Max \\|residual\\| | {np.max(np.abs(residuals)):.4f} |",
        ]
        L += _quality_report_rows(fit_results.get('fit_quality'))
        L.append("")

        # ── Level parameters ──────────────────────────────────────────────
        for i, level in enumerate(levels):
            lnum = i + 1
            L.append(f"## Level {lnum} Parameters")
            L.append("")

            has_mc = any(f'{p}_err' in level for p in ('G', 'Rg', 'B', 'P', 'ETA', 'PACK'))
            if has_mc:
                L += ["| Parameter | Value | Uncertainty (1σ) |",
                      "|-----------|-------|------------------|"]
            else:
                L += ["| Parameter | Value |",
                      "|-----------|-------|"]

            def _row(label, key, unit='', fmt='.4e'):
                val = level.get(key)
                if val is None:
                    return
                if fmt.endswith('f'):
                    val_str = f"{val:{fmt}}{unit}"
                else:
                    val_str = f"{_ef(float(val))}{unit}"
                if has_mc:
                    err = level.get(f'{key}_err', 0.0)
                    if fmt.endswith('f'):
                        err_str = f"± {err:{fmt}}{unit}" if err > 0 else "—"
                    else:
                        err_str = f"± {_ef(float(err))}{unit}" if err > 0 else "—"
                    L.append(f"| {label} | {val_str} | {err_str} |")
                else:
                    L.append(f"| {label} | {val_str} |")

            _row('G',  'G')
            _row('Rg', 'Rg', ' Å')
            _row('B',  'B')
            _row('P',  'P',  fmt='.4f')

            rgcut = level.get('RgCutoff', 0.0)
            if isinstance(rgcut, float) and rgcut > 0.01:
                _row('RgCutoff', 'RgCutoff', ' Å')

            if level.get('correlated', False):
                _row('ETA',  'ETA',  ' Å', fmt='.2f')
                _row('PACK', 'PACK', '',   fmt='.4f')

            _row('Sv',        'Sv',        ' m²/cm³')
            _row('Invariant', 'Invariant', ' cm⁻⁴')

            L.append("")

    # ── Size Distribution results ─────────────────────────────────────────────
    if sizes_results is not None:
        # All scalar metadata is at the top level of the dict returned by
        # load_sizes_results() — not nested under a 'params' key.
        r_grid    = sizes_results.get('r_grid')
        dist      = sizes_results.get('distribution')
        residuals = sizes_results.get('residuals')

        def _sv(key, default=float('nan')):
            """Return scalar from sizes_results, substituting default for None."""
            v = sizes_results.get(key)
            return v if v is not None else default

        chi2        = _sv('chi_squared')
        vf          = _sv('volume_fraction')
        rg          = _sv('rg')
        n_iter      = _sv('n_iterations', 'N/A')
        method      = _sv('method', 'unknown')
        shape       = _sv('shape', 'unknown')
        contrast    = _sv('contrast')
        aspect_ratio = _sv('aspect_ratio')
        n_bins      = _sv('n_bins', 0)
        r_min       = _sv('r_min')
        r_max       = _sv('r_max')
        log_spacing = _sv('log_spacing', False)
        background  = _sv('background')
        error_scale = _sv('error_scale')
        power_law_B = _sv('power_law_B')
        power_law_P = _sv('power_law_P')
        q_fit_min   = sizes_results.get('cursor_q_min')
        q_fit_max   = sizes_results.get('cursor_q_max')

        peak_r = float('nan')
        if dist is not None and r_grid is not None and len(dist) > 0:
            peak_r = float(r_grid[int(np.argmax(dist))])

        def _fmt(v, spec='.4g', suffix=''):
            """Format a value, returning 'N/A' for None/nan."""
            if v is None:
                return 'N/A'
            try:
                fv = float(v)
                if np.isnan(fv):
                    return 'N/A'
            except (TypeError, ValueError):
                return str(v)
            try:
                if spec.endswith('f'):
                    return f"{fv:{spec}}{suffix}"
                sig = int(spec[1:-1]) if len(spec) > 2 else 4
                return _ef(fv, sig=sig) + suffix
            except (TypeError, ValueError):
                return str(v)

        L += [
            "## Size Distribution",
            "",
            "**Fit results:**",
            "",
            "| Parameter | Value |",
            "|-----------|-------|",
            f"| Method | {method} |",
            f"| Particle shape | {shape} |",
            f"| Chi-squared (χ²) | {_fmt(chi2)} |",
            f"| Volume fraction | {_fmt(vf, '.4e')} |",
            f"| Rg | {_fmt(rg)} Å |",
            f"| Peak r | {_fmt(peak_r)} Å |",
            f"| Iterations | {n_iter} |",
        ]
        if residuals is not None:
            L += [
                f"| Residuals mean | {np.mean(residuals):.4f} |",
                f"| Residuals std dev | {np.std(residuals):.4f} |",
            ]
        L += _quality_report_rows(sizes_results.get('fit_quality'))
        L.append("")

        L += [
            "**Grid / model setup:**",
            "",
            "| Parameter | Value |",
            "|-----------|-------|",
            f"| r min | {_fmt(r_min)} Å |",
            f"| r max | {_fmt(r_max)} Å |",
            f"| Bins | {n_bins} |",
            f"| Log spacing | {log_spacing} |",
            f"| Contrast (Δρ)² | {_fmt(contrast)} ×10²⁰ cm⁻⁴ |",
        ]
        try:
            if aspect_ratio is not None and not np.isnan(float(aspect_ratio)) and float(aspect_ratio) != 1.0:
                L.append(f"| Aspect ratio | {_fmt(aspect_ratio, '.4f')} |")
        except (TypeError, ValueError):
            log.debug("suppressed exception", exc_info=True)
        L += [
            f"| Background | {_fmt(background, '.4e')} cm⁻¹ |",
            f"| Error scale | {_fmt(error_scale)} |",
        ]
        try:
            if power_law_B is not None and not np.isnan(float(power_law_B)) and float(power_law_B) != 0.0:
                L.append(f"| Power law B | {_fmt(power_law_B, '.4e')} |")
                L.append(f"| Power law P | {_fmt(power_law_P, '.4f')} |")
        except (TypeError, ValueError):
            log.debug("suppressed exception", exc_info=True)
        if q_fit_min is not None and q_fit_max is not None:
            L.append(f"| Q range (fit) | {_fmt(q_fit_min)} – {_fmt(q_fit_max)} Å⁻¹ |")
        L.append("")

        # Method-specific parameters
        if str(method).lower() == 'maxent':
            sky      = _sv('maxent_sky_background')
            stab     = _sv('maxent_stability')
            max_iter = _sv('maxent_max_iter', 'N/A')
            L += [
                "**MaxEnt parameters:**",
                "",
                "| Parameter | Value |",
                "|-----------|-------|",
                f"| Sky background | {_fmt(sky)} |",
                f"| Stability | {_fmt(stab)} |",
                f"| Max iterations | {max_iter} |",
                "",
            ]
        elif str(method).lower() == 'regularization':
            evalue    = _sv('regularization_evalue')
            min_ratio = _sv('regularization_min_ratio')
            L += [
                "**Regularization parameters:**",
                "",
                "| Parameter | Value |",
                "|-----------|-------|",
                f"| Eigenvalue weight | {_fmt(evalue)} |",
                f"| Min ratio | {_fmt(min_ratio)} |",
                "",
            ]
        elif str(method).lower() == 'tnnls':
            approach = _sv('tnnls_approach_param')
            max_iter = _sv('tnnls_max_iter', 'N/A')
            L += [
                "**TNNLS parameters:**",
                "",
                "| Parameter | Value |",
                "|-----------|-------|",
                f"| Approach parameter | {_fmt(approach)} |",
                f"| Max iterations | {max_iter} |",
                "",
            ]
        elif str(method).lower() in ('montecarlo', 'mcsas'):
            n_rep    = _sv('montecarlo_n_repetitions') or _sv('mcsas_n_repetitions', 'N/A')
            conv     = _sv('montecarlo_convergence') or _sv('mcsas_convergence')
            max_iter = _sv('montecarlo_max_iter') or _sv('mcsas_max_iter', 'N/A')
            L += [
                "**Monte Carlo parameters:**",
                "",
                "| Parameter | Value |",
                "|-----------|-------|",
                f"| Repetitions | {n_rep} |",
                f"| Convergence | {_fmt(conv)} |",
                f"| Max iterations | {max_iter} |",
                "",
            ]

    # ── Simple Fits results ───────────────────────────────────────────────────
    if simple_fit_results is not None:
        sf_model    = simple_fit_results.get('model', 'unknown')
        sf_chi2     = simple_fit_results.get('chi_squared')
        sf_rchi2    = simple_fit_results.get('reduced_chi_squared')
        sf_dof      = simple_fit_results.get('dof')
        sf_q_min    = simple_fit_results.get('q_min')
        sf_q_max    = simple_fit_results.get('q_max')
        sf_complex  = simple_fit_results.get('use_complex_bg', False)
        sf_params   = simple_fit_results.get('params', {})
        sf_std      = simple_fit_results.get('params_std', {})
        sf_derived  = simple_fit_results.get('derived', {})

        def _sf_fmt(v, spec='.4g', suffix=''):
            if v is None:
                return 'N/A'
            try:
                fv = float(v)
                if np.isnan(fv):
                    return 'N/A'
            except (TypeError, ValueError):
                return str(v)
            try:
                if spec.endswith('f'):
                    return f"{fv:{spec}}{suffix}"
                sig = int(spec[1:-1]) if len(spec) > 2 else 4
                return _ef(fv, sig=sig) + suffix
            except (TypeError, ValueError):
                return str(v)

        L += [
            "## Simple Fits",
            "",
            "| Parameter | Value |",
            "|-----------|-------|",
            f"| Model | {sf_model} |",
            f"| Chi-squared (χ²) | {_sf_fmt(sf_chi2)} |",
            f"| Reduced chi² | {_sf_fmt(sf_rchi2)} |",
            f"| DOF | {sf_dof if sf_dof is not None else 'N/A'} |",
            f"| Q range (fit) | {_sf_fmt(sf_q_min)} – {_sf_fmt(sf_q_max)} Å⁻¹ |",
            f"| Complex background | {sf_complex} |",
        ]
        L += _quality_report_rows(simple_fit_results.get('fit_quality'))
        L.append("")

        if sf_params:
            has_std = bool(sf_std)
            if has_std:
                L += ["**Parameters:**", "",
                      "| Parameter | Value | Uncertainty (1σ) |",
                      "|-----------|-------|------------------|"]
            else:
                L += ["**Parameters:**", "",
                      "| Parameter | Value |",
                      "|-----------|-------|"]
            for name, val in sf_params.items():
                val_str = _sf_fmt(val, '.6g')
                if has_std:
                    err = sf_std.get(name)
                    err_str = f"± {_sf_fmt(err, '.3g')}" if err is not None else "—"
                    L.append(f"| {name} | {val_str} | {err_str} |")
                else:
                    L.append(f"| {name} | {val_str} |")
            L.append("")

        if sf_derived:
            L += ["**Derived quantities:**", "",
                  "| Quantity | Value |",
                  "|----------|-------|"]
            for name, val in sf_derived.items():
                L.append(f"| {name} | {_sf_fmt(val, '.6g')} |")
            L.append("")

    # ── WAXS Peak Fit ─────────────────────────────────────────────────────────
    if waxs_peakfit_results is not None:
        wp = waxs_peakfit_results

        def _wp_fmt(v, spec='.4g'):
            if v is None:
                return 'N/A'
            try:
                fv = float(v)
                if np.isnan(fv):
                    return 'N/A'
            except (TypeError, ValueError):
                return str(v)
            try:
                if spec.endswith('f'):
                    return f"{fv:{spec}}"
                sig = int(spec[1:-1]) if len(spec) > 2 else 4
                return _ef(fv, sig=sig)
            except (TypeError, ValueError):
                return str(v)

        L += [
            "## WAXS Peak Fit",
            "",
            "| Parameter | Value |",
            "|-----------|-------|",
            f"| Number of peaks | {wp.get('n_peaks', 0)} |",
            f"| Background shape | {wp.get('bg_shape', 'N/A')} |",
            f"| Chi-squared (χ²) | {_wp_fmt(wp.get('chi_squared'))} |",
            f"| Reduced chi² | {_wp_fmt(wp.get('reduced_chi_squared'))} |",
            f"| DOF | {wp.get('dof', 'N/A')} |",
            f"| Q range (fit) | {_wp_fmt(wp.get('q_min'))} – {_wp_fmt(wp.get('q_max'))} Å⁻¹ |",
        ]
        L += _quality_report_rows(wp.get('fit_quality'))
        L.append("")

        peaks_list = wp.get('peaks', [])
        peaks_std  = wp.get('peaks_std', [])
        for i, pk in enumerate(peaks_list):
            pstd = peaks_std[i] if i < len(peaks_std) else {}
            shape = pk.get('shape', 'Gauss')
            L += [
                f"### Peak {i + 1} ({shape})",
                "",
                "| Parameter | Value | Uncertainty (1σ) |",
                "|-----------|-------|------------------|",
            ]
            for pname in ('A', 'Q0', 'FWHM', 'eta'):
                if pname not in pk:
                    continue
                val = pk[pname]
                val_v = val.get('value') if isinstance(val, dict) else val
                std_v = pstd.get(pname)
                unit = ' Å⁻¹' if pname in ('Q0', 'FWHM') else ''
                std_str = f"± {_wp_fmt(std_v, '.3g')}" if std_v is not None else "—"
                L.append(f"| {pname}{unit} | {_wp_fmt(val_v, '.6g')} | {std_str} |")
            # Derived: integral under the peak (area)
            area_v = pk.get('area')
            area_s = pk.get('area_std')
            if area_v is None:
                # Older HDF5 files (pre-Area) — compute on the fly so reports
                # always carry the derived value.
                from pyirena.core.waxs_peakfit import peak_area, peak_area_std
                area_v = peak_area(shape, pk)
                area_s = peak_area_std(shape, pk, pstd)
            area_std_str = f"± {_wp_fmt(area_s, '.3g')}" if area_s is not None else "—"
            L.append(f"| Area (∫ peak dq) | {_wp_fmt(area_v, '.6g')} | {area_std_str} |")
            L.append("")

    if modeling_results is not None:
        mr = modeling_results

        def _mr_fmt(v, spec='.4g'):
            if v is None:
                return 'N/A'
            try:
                fv = float(v)
                if np.isnan(fv):
                    return 'N/A'
            except (TypeError, ValueError):
                return str(v)
            try:
                if spec.endswith('f'):
                    return f"{fv:{spec}}"
                sig = int(spec[1:-1]) if len(spec) > 2 else 4
                return _ef(fv, sig=sig)
            except (TypeError, ValueError):
                return str(v)

        pops_all = mr.get('populations', [])
        pops = [p for p in pops_all if p.get('enabled', True)]
        L += [
            "## Modeling Results",
            "",
            "| Parameter | Value |",
            "|-----------|-------|",
            f"| Chi-squared (χ²) | {_mr_fmt(mr.get('chi_squared'))} |",
            f"| Background | {_mr_fmt(mr.get('background'))} |",
            f"| Populations (enabled) | {len(pops)} |",
            f"| Q min | {_mr_fmt(mr.get('q_min'))} Å⁻¹ |",
            f"| Q max | {_mr_fmt(mr.get('q_max'))} Å⁻¹ |",
        ]
        L += _quality_report_rows(mr.get('fit_quality'))
        L.append("")
        for pop in pops:
            pt    = pop.get('pop_type', 'size_dist')
            idx   = pop.get('population_index', '?')
            label = pop.get('label', '') or f'P{idx}'
            L += [f"### Population {idx} — {label} ({pt})", ""]
            derived = pop.get('derived', {})
            # UF, peak, and size_dist params are all flat keys on the pop dict
            if pt == 'unified_level':
                rows = [("G [cm⁻¹]", pop.get("G")), ("Rg [Å]", pop.get("Rg")),
                        ("B", pop.get("B")), ("P", pop.get("P")),
                        ("RgCO [Å]", pop.get("RgCO")),
                        ("ETA [Å]", pop.get("ETA")), ("PACK", pop.get("PACK"))]
            elif pt == 'diffraction_peak':
                rows = [("Peak type", pop.get("peak_type")),
                        ("Position Q₀ [Å⁻¹]", pop.get("position")),
                        ("Amplitude [cm⁻¹]", pop.get("amplitude")),
                        ("Width σ [Å⁻¹]", pop.get("width")),
                        ("η (Voigt)", pop.get("eta_voigt"))]
            elif pt == 'guinier_porod':
                rows = [("G [cm⁻¹]", pop.get("G")), ("Rg1 [Å]", pop.get("Rg1")),
                        ("Slope s1", pop.get("s1")), ("Power P", pop.get("P")),
                        ("Rg2 [Å]", pop.get("Rg2")), ("Slope s2", pop.get("s2")),
                        ("RgCO [Å]", pop.get("RgCO") or None),
                        ("ETA [Å]", pop.get("ETA") if pop.get("correlations") else None),
                        ("PACK", pop.get("PACK") if pop.get("correlations") else None)]
                rows = [(k, v) for k, v in rows if v is not None]
            elif pt == 'mass_fractal':
                rows = [("Phi (vol. frac.)", pop.get("Phi")),
                        ("Radius [Å]", pop.get("Radius")),
                        ("Fractal dim. Dv", pop.get("Dv")),
                        ("Ksi [Å]", pop.get("Ksi")),
                        ("Eta", pop.get("Eta")),
                        ("Contrast", pop.get("Contrast"))]
            elif pt == 'surface_fractal':
                rows = [("Surface [cm⁻¹]", pop.get("Surface")),
                        ("Fractal dim. Ds", pop.get("Ds")),
                        ("Ksi [Å]", pop.get("Ksi")),
                        ("Contrast", pop.get("Contrast")),
                        ("Qc [Å⁻¹]", pop.get("Qc") if pop.get("use_porod_transition") else None)]
                rows = [(k, v) for k, v in rows if v is not None]
            else:  # size_dist
                rows = [("Distribution", pop.get("dist_type"))]
                # Raw distribution parameters (mean_size, width, sdeviation, …)
                for pn, pv in (pop.get("dist_params") or {}).items():
                    rows.append((f"Dist: {pn}", pv))
                rows += [("Scale", pop.get("scale")),
                         ("Contrast [10²⁰ cm⁻⁴]", pop.get("contrast")),
                         ("Form factor", pop.get("form_factor"))]
                # Form-factor parameters (SLDs, shell thicknesses, aspect ratio, …)
                for pn, pv in (pop.get("ff_params") or {}).items():
                    rows.append((f"FF: {pn}", pv))
                sf_name = pop.get("structure_factor")
                rows.append(("Structure factor", sf_name))
                # Structure-factor parameters (only meaningful when one is active)
                if sf_name and str(sf_name).lower() != 'none':
                    for pn, pv in (pop.get("sf_params") or {}).items():
                        rows.append((f"SF: {pn}", pv))
                rows += [("Vol. fraction", derived.get("volume_fraction")),
                         ("Mean radius [Å]", derived.get("vol_mean_r"))]
            L += ["| Parameter | Value |", "|-----------|-------|"]
            for rname, rval in rows:
                if rval is not None:
                    L.append(f"| {rname} | {_mr_fmt(rval)} |")
            L.append("")

    if saxs_morph_results is not None:
        sm = saxs_morph_results

        def _sm_fmt(v, spec='.4g'):
            if v is None:
                return 'N/A'
            try:
                fv = float(v)
                if np.isnan(fv):
                    return 'N/A'
            except (TypeError, ValueError):
                return str(v)
            try:
                if spec.endswith('f'):
                    return f"{fv:{spec}}"
                sig = int(spec[1:-1]) if len(spec) > 2 else 4
                return _ef(fv, sig=sig)
            except (TypeError, ValueError):
                return str(v)

        L += [
            "## SAXS Morph Results",
            "",
            "| Parameter | Value |",
            "|-----------|-------|",
            f"| Chi-squared (χ²) | {_sm_fmt(sm.get('chi_squared'))} |",
            f"| Volume fraction φ | {_sm_fmt(sm.get('volume_fraction'))} |",
            f"| Contrast (Δρ)² [10²⁰ cm⁻⁴] | {_sm_fmt(sm.get('contrast'))} |",
            f"| Rg from γ(r) [Å] | {_sm_fmt(sm.get('rg_A'))} |",
            f"| Specific surface area S/V [Å⁻¹] | {_sm_fmt(sm.get('specific_surface_area_inv_A'))} |",
            f"| Voxel cube side | {sm.get('voxel_size', 'N/A')}³ |",
            f"| Box size [Å] | {_sm_fmt(sm.get('box_size_A'))} |",
            f"| Voxel pitch [Å] | {_sm_fmt(sm.get('voxel_pitch_A'))} |",
            "",
        ]
        mm = sm.get('morphology_metrics')
        if mm is not None:
            # mm may be a MorphologyMetrics dataclass or a dict (legacy)
            def _g(name, default=None):
                if hasattr(mm, name):
                    return getattr(mm, name)
                if isinstance(mm, dict):
                    return mm.get(name, default)
                return default

            open_pct = _g('open_porosity_fraction')
            closed_pct = _g('closed_porosity_fraction')
            open_str = (_sm_fmt(open_pct * 100) + ' %'
                        if open_pct is not None else 'N/A')
            closed_str = (_sm_fmt(closed_pct * 100) + ' %'
                          if closed_pct is not None else 'N/A')

            L += [
                "### Morphology of minority phase",
                "",
                "*Single realisation — depends on RNG seed and voxel pitch.*",
                "",
                "| Metric | Value |",
                "|--------|-------|",
                f"| Minority phase value | {_g('minority_phase_value', 'N/A')} "
                f"(φ = {_sm_fmt(_g('minority_volume_fraction'))}) |",
                f"| # of independent clusters | {_g('n_clusters', 'N/A')} |",
                f"| Open porosity (touching boundary) | {open_str} |",
                f"| Closed porosity (isolated) | {closed_str} |",
                f"| Percolating along X | {'yes' if _g('percolating_x') else 'no'} |",
                f"| Percolating along Y | {'yes' if _g('percolating_y') else 'no'} |",
                f"| Percolating along Z | {'yes' if _g('percolating_z') else 'no'} |",
                f"| Euler-Poincaré characteristic χ | {_g('euler_number', 'N/A')} |",
                f"| Pore radius Q1 (25%) [Å] | {_sm_fmt(_g('pore_size_q25_A'))} |",
                f"| Pore radius median (50%) [Å] | {_sm_fmt(_g('pore_size_median_A'))} |",
                f"| Pore radius Q3 (75%) [Å] | {_sm_fmt(_g('pore_size_q75_A'))} |",
                "",
            ]

    L += ["---", "*Generated by pyIrena*", ""]
    return "\n".join(L)
