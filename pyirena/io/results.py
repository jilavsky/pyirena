"""
High-level results loader for pyIrena analysis packages.

Provides :func:`load_result` — a single importable entry point that loads
stored fit results from an NXcanSAS HDF5 file and returns them as a fully
documented dictionary.

Supported analysis packages
---------------------------
``'unified_fit'``
    Unified Fit (Beaucage) results stored under ``entry/unified_fit_results``.

``'size_distribution'``
    Size distribution fitting results stored under ``entry/sizes_results``.

Usage
-----
::

    from pyirena import load_result

    # Load Unified Fit results
    r = load_result("mydata.h5", "unified_fit")
    if r["found"]:
        print(f"chi² = {r['chi_squared']:.4f}")
        for i, lv in enumerate(r["levels"], 1):
            print(f"  Level {i}: Rg = {lv['Rg']:.2f} Å")

    # Load Size Distribution results
    r = load_result("mydata.h5", "size_distribution")
    if r["found"]:
        import numpy as np
        peak_r = r["r_grid"][np.argmax(r["distribution"])]
        print(f"Peak radius = {peak_r:.1f} Å,  Vf = {r['volume_fraction']:.4g}")

    # Non-existent results return an empty structure — no exception raised
    r = load_result("mydata.h5", "size_distribution")
    print(r["found"])       # False
    print(r["r_grid"])      # None
    print(r["chi_squared"]) # None
"""

from __future__ import annotations

from pathlib import Path
from typing import Union


# ──────────────────────────────────────────────────────────────────────────────
# Public API
# ──────────────────────────────────────────────────────────────────────────────

#: Names of the analysis packages recognised by :func:`load_result`.
SUPPORTED_ANALYSES = ('unified_fit', 'size_distribution')


def load_result(
    filepath: Union[str, Path],
    analysis: str,
) -> dict:
    """
    Load stored fit results from an NXcanSAS HDF5 file.

    Returns a fully populated dictionary for the requested analysis package.
    If the file does not contain results for that package (or if the file does
    not exist) the dictionary is still returned — every key is present but set
    to ``None`` (scalars) or ``None`` (arrays), and ``result["found"]`` is
    ``False``.  No exception is raised for missing data.

    Parameters
    ----------
    filepath : str or Path
        Path to the NXcanSAS HDF5 file produced by pyIrena.
    analysis : str
        Which analysis package to load.  One of:

        ``'unified_fit'``
            Unified Fit (Beaucage) results.
        ``'size_distribution'``
            Size distribution fitting results.

    Returns
    -------
    dict
        For ``'unified_fit'`` the dict contains:

        Administrative
            ``found`` *(bool)* — ``True`` if results exist in the file.
            ``timestamp`` *(str | None)* — ISO-format save timestamp.
            ``program`` *(str | None)* — program name (e.g. ``'pyirena.core.unified'``).

        Arrays (``None`` when not found)
            ``Q`` — Q vector [Å⁻¹].
            ``intensity_data`` — measured intensity [cm⁻¹].
            ``intensity_model`` — fitted model intensity [cm⁻¹].
            ``intensity_error`` — measurement uncertainty [cm⁻¹] (may be ``None``
            even when ``found`` is ``True`` if errors were not stored).
            ``residuals`` — normalised residuals (data − model) / error.

        Scalars (``None`` when not found)
            ``chi_squared`` — goodness-of-fit χ².
            ``background`` — flat background value [cm⁻¹].
            ``background_err`` — MC uncertainty on background (``None`` if not run).
            ``num_levels`` — number of Unified Fit levels.

        Level parameters
            ``levels`` *(list of dict)* — one dict per level in order.  Each dict
            may contain:
            ``G``, ``Rg`` [Å], ``B``, ``P``, ``RgCutoff`` [Å], ``ETA`` [Å],
            ``PACK``, ``correlated`` (bool), ``level_number``;
            and optionally MC uncertainty counterparts
            ``G_err``, ``Rg_err``, ``B_err``, ``P_err``, ``ETA_err``, ``PACK_err``.
            ``levels`` is ``[]`` when ``found`` is ``False``.

        For ``'size_distribution'`` the dict contains:

        Administrative
            ``found``, ``timestamp``, ``program`` — same as above.

        Arrays (``None`` when not found or not stored)
            ``Q``, ``intensity_data``, ``intensity_model``, ``intensity_error``,
            ``residuals`` — same meaning as above.
            ``r_grid`` — radius bin centres [Å].
            ``distribution`` — volume-fraction size distribution P(r) [Å⁻¹].
            ``distribution_std`` — per-bin 1σ uncertainty from McSAS repetitions
            (``None`` unless McSAS was the fitting method).

        Fit results (``None`` when not found)
            ``chi_squared``, ``volume_fraction``, ``rg`` [Å], ``n_iterations``,
            ``q_power``.

        Model / grid setup (``None`` when not found)
            ``shape`` — particle shape (``'sphere'`` or ``'spheroid'``).
            ``contrast`` — (Δρ)² [10²⁰ cm⁻⁴].
            ``aspect_ratio`` — spheroid polar/equatorial axis ratio.
            ``r_min``, ``r_max`` — radius grid bounds [Å].
            ``n_bins`` — number of radius bins.
            ``log_spacing`` — ``True`` for logarithmic bin spacing.
            ``background`` — flat background [cm⁻¹].
            ``power_law_B``, ``power_law_P`` — power-law amplitude and exponent.
            ``method`` — fitting method (``'maxent'``, ``'regularization'``,
            ``'tnnls'``, or ``'mcsas'``).
            ``error_scale`` — multiplicative scale applied to measurement errors.

        Method-specific parameters (``None`` when not found or not applicable)
            MaxEnt: ``maxent_sky_background``, ``maxent_stability``,
            ``maxent_max_iter``.
            Regularization: ``regularization_evalue``,
            ``regularization_min_ratio``.
            TNNLS: ``tnnls_approach_param``, ``tnnls_max_iter``.
            McSAS: ``mcsas_n_repetitions``, ``mcsas_convergence``,
            ``mcsas_max_iter``.

        Q ranges (``None`` when not found or not set)
            ``power_law_q_min``, ``power_law_q_max`` — Q range for power-law fit.
            ``background_q_min``, ``background_q_max`` — Q range for background fit.
            ``cursor_q_min``, ``cursor_q_max`` — Q range used for the fit.

    Raises
    ------
    ValueError
        If *analysis* is not one of the supported package names.
    OSError
        If *filepath* exists but cannot be opened as an HDF5 file.

    Examples
    --------
    Batch-process many files and collect volume fractions::

        from pathlib import Path
        from pyirena import load_result

        files = list(Path("data/").glob("*.h5"))
        vf_values = []
        for f in files:
            r = load_result(f, "size_distribution")
            if r["found"]:
                vf_values.append((f.name, r["volume_fraction"]))

    Access level parameters from a Unified Fit::

        from pyirena import load_result

        r = load_result("sample.h5", "unified_fit")
        if r["found"]:
            lv1 = r["levels"][0]
            print(f"Rg = {lv1['Rg']:.1f} ± {lv1.get('Rg_err', 0):.1f} Å")
    """
    analysis = str(analysis).strip().lower()
    if analysis not in SUPPORTED_ANALYSES:
        raise ValueError(
            f"Unknown analysis '{analysis}'. "
            f"Supported values: {SUPPORTED_ANALYSES}"
        )

    filepath = Path(filepath)

    if analysis == 'unified_fit':
        return _load_unified_fit(filepath)
    else:
        return _load_size_distribution(filepath)


# ──────────────────────────────────────────────────────────────────────────────
# Internal helpers
# ──────────────────────────────────────────────────────────────────────────────

def _empty_unified_fit() -> dict:
    """Return the canonical empty structure for a unified_fit result."""
    return {
        # Administrative
        'found':          False,
        'timestamp':      None,
        'program':        None,
        # Scalars
        'chi_squared':    None,
        'background':     None,
        'background_err': None,
        'num_levels':     None,
        # Arrays
        'Q':               None,
        'intensity_data':  None,
        'intensity_model': None,
        'intensity_error': None,
        'residuals':       None,
        # Level parameters
        'levels':          [],
    }


def _empty_size_distribution() -> dict:
    """Return the canonical empty structure for a size_distribution result."""
    return {
        # Administrative
        'found':     False,
        'timestamp': None,
        'program':   None,
        # Arrays
        'Q':                None,
        'intensity_data':   None,
        'intensity_model':  None,
        'intensity_error':  None,
        'residuals':        None,
        'r_grid':           None,
        'distribution':     None,
        'distribution_std': None,
        # Fit results
        'chi_squared':      None,
        'volume_fraction':  None,
        'rg':               None,
        'n_iterations':     None,
        'q_power':          None,
        # Model / grid setup
        'shape':            None,
        'contrast':         None,
        'aspect_ratio':     None,
        'r_min':            None,
        'r_max':            None,
        'n_bins':           None,
        'log_spacing':      None,
        'background':       None,
        'power_law_B':      None,
        'power_law_P':      None,
        'method':           None,
        'error_scale':      None,
        # Method-specific — MaxEnt
        'maxent_sky_background':    None,
        'maxent_stability':         None,
        'maxent_max_iter':          None,
        # Method-specific — Regularization
        'regularization_evalue':    None,
        'regularization_min_ratio': None,
        # Method-specific — TNNLS
        'tnnls_approach_param':     None,
        'tnnls_max_iter':           None,
        # Method-specific — McSAS
        'mcsas_n_repetitions':      None,
        'mcsas_convergence':        None,
        'mcsas_max_iter':           None,
        # Q ranges
        'power_law_q_min':   None,
        'power_law_q_max':   None,
        'background_q_min':  None,
        'background_q_max':  None,
        'cursor_q_min':      None,
        'cursor_q_max':      None,
    }


def _load_unified_fit(filepath: Path) -> dict:
    """Load Unified Fit results; return empty structure on failure."""
    from pyirena.io.nxcansas_unified import load_unified_fit_results

    result = _empty_unified_fit()

    if not filepath.exists():
        return result

    try:
        raw = load_unified_fit_results(filepath)
    except (ValueError, KeyError, OSError):
        # Group not present or file has no Unified Fit data
        return result

    result['found']          = True
    result['timestamp']      = raw.get('timestamp')
    result['program']        = raw.get('program')
    result['chi_squared']    = raw.get('chi_squared')
    result['background']     = raw.get('background')
    result['background_err'] = raw.get('background_err')  # set by MC; may be None
    result['num_levels']     = raw.get('num_levels')
    result['Q']               = raw.get('Q')
    result['intensity_data']  = raw.get('intensity_data')
    result['intensity_model'] = raw.get('intensity_model')
    result['intensity_error'] = raw.get('intensity_error')
    result['residuals']       = raw.get('residuals')
    result['levels']          = raw.get('levels', [])

    return result


def _load_size_distribution(filepath: Path) -> dict:
    """Load size distribution results; return empty structure on failure."""
    from pyirena.io.nxcansas_sizes import load_sizes_results

    result = _empty_size_distribution()

    if not filepath.exists():
        return result

    try:
        raw = load_sizes_results(filepath)
    except (KeyError, OSError):
        # Group not present or file has no size distribution data
        return result

    result['found'] = True

    # Copy every key that exists in the empty template directly from raw,
    # skipping 'found' itself (already set above).
    for key in result:
        if key == 'found':
            continue
        if key in raw:
            result[key] = raw[key]

    return result
