"""
pyirena.batch.convert — Igor pxp/h5xp -> NXcanSAS conversion (igor_to_nexus, pxp_to_nexus).

Split from the original monolithic pyirena/batch.py (no behavior change).
"""

from __future__ import annotations

import logging
import traceback
from pathlib import Path
from typing import TYPE_CHECKING, Dict, List, Optional, Union


from pyirena.logging_setup import ensure_console_output as _ensure_console

if TYPE_CHECKING:
    pass

log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Igor packed experiment (.pxp) import
# ---------------------------------------------------------------------------

def igor_to_nexus(
    igor_file: Union[str, Path],
    output_folder: Optional[Union[str, Path]] = None,
    techniques: Optional[List[str]] = None,
    overwrite: bool = False,
    verbose: bool = False,
) -> Optional[Dict]:
    """Import an Igor Pro packed experiment (.pxp **or** .h5xp) and export
    its reduced data as per-sample NeXus files.

    Auto-detects the format from the file extension and dispatches to the
    appropriate reader. Both formats produce the same output: one
    NXcanSAS ``.h5`` per sample under
    ``<output>/<technique>/<sample>.h5`` (see
    :func:`pyirena.io.pxp_to_nexus.extract_igor_experiment` for the
    folder and wave-name conventions recognised).

    Args:
        igor_file: Path to the input ``.pxp`` (binary Igor packed) or
            ``.h5xp`` (Wavemetrics HDF5 packed) experiment.
        output_folder: Destination directory. Defaults to ``<stem>_data``
            next to the input file. Created if it does not exist.
        techniques: List of techniques to export, e.g. ``["USAXS", "SAXS"]``.
            ``None`` (default) exports every technique present in the file.
        overwrite: If True, existing files in the destination are overwritten.
            If False (default), an unused suffix is appended.
        verbose: Print per-file summary to stdout.

    Returns:
        Dict with keys ``success`` (always True if the function returns),
        ``output_folder``, ``n_written``, ``n_skipped``, ``n_errors``,
        ``files`` (list of per-file dicts: ``source``, ``output``,
        ``technique``, ``n_points``, ``status``, ``message``).
        Returns ``None`` only if the input file cannot be located or its
        extension is unrecognised.

    Example::

        from pyirena.batch import igor_to_nexus
        result = igor_to_nexus("legacy_experiment.pxp", techniques=["USAXS"])
        print(f"Wrote {result['n_written']} files to {result['output_folder']}")

        # h5xp works the same way
        igor_to_nexus("modern_export.h5xp")
    """
    _ensure_console()
    from pyirena.io.pxp_to_nexus import extract_igor_experiment

    src_path = Path(igor_file)
    if not src_path.is_file():
        log.error(f"[pyirena.batch.igor_to_nexus] File not found: {src_path}")
        return None

    try:
        res = extract_igor_experiment(
            path=src_path,
            output_root=output_folder,
            techniques=techniques,
            overwrite=overwrite,
        )
    except ValueError as exc:
        # Unsupported extension — translate to None-return per the
        # batch-API convention (caller already prints the path).
        log.info(f"[pyirena.batch.igor_to_nexus] {exc}")
        return None
    except Exception:
        log.error(f"[pyirena.batch.igor_to_nexus] Error:\n{traceback.format_exc()}")
        return None

    files_list = [
        {
            'source':    f.source_path,
            'output':    str(f.output_path) if f.output_path else None,
            'technique': f.technique,
            'n_points':  f.n_points,
            'status':    f.status,
            'message':   f.message,
        }
        for f in res.files
    ]

    if verbose:
        log.info(f"[pyirena.batch.igor_to_nexus] {res.pxp_path}")
        log.info(f"  output:   {res.output_root}")
        log.info(f"  written:  {res.n_written}")
        log.info(f"  skipped:  {res.n_skipped}")
        log.error(f"  errors:   {res.n_errors}")
        if res.n_unparseable_records:
            log.error(f"  note:     {res.n_unparseable_records} wave record(s) could not be parsed")
        if res.n_igor8_longname_markers:
            log.warning(f"  WARNING:  {res.n_igor8_longname_markers} Igor-8 long-name "
                  f"record(s) seen — some samples likely missing. "
                  f"Re-save as .h5xp from Igor.")

    return {
        'success':     True,
        'output_folder': str(res.output_root),
        'n_written':   res.n_written,
        'n_skipped':   res.n_skipped,
        'n_errors':    res.n_errors,
        'n_unparseable_records': res.n_unparseable_records,
        'n_igor8_longname_markers': res.n_igor8_longname_markers,
        'files':       files_list,
    }


def pxp_to_nexus(
    pxp_file: Union[str, Path],
    output_folder: Optional[Union[str, Path]] = None,
    techniques: Optional[List[str]] = None,
    overwrite: bool = False,
    verbose: bool = False,
) -> Optional[Dict]:
    """Backwards-compatible alias for :func:`igor_to_nexus`.

    .. deprecated:: 0.8.1
       Renamed to :func:`igor_to_nexus` now that .h5xp is also supported.
       The old name is kept so existing scripts continue to work; new
       code should prefer ``igor_to_nexus``.
    """
    _ensure_console()
    return igor_to_nexus(
        igor_file=pxp_file,
        output_folder=output_folder,
        techniques=techniques,
        overwrite=overwrite,
        verbose=verbose,
    )
