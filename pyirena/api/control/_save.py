"""Shared save-target resolution for the control API's ``save_*`` tools.

All six fitting tools answer the same three questions before writing results:
is the target inside ``PYIRENA_DATA_ROOT``, does a file exist there to write
into, and — when saving somewhere new — how is that file seeded so it reopens
as a complete NXcanSAS file rather than a results-only stub. They used to
answer them with six copies of the same block; this module is the single copy.

The third question grew a second answer when sessions stopped always having a
file behind them. A session created by ``open_dataset_from_data`` (the remote
ZMQ service, or any script that has arrays rather than a path) has
``file_path is None``: there is nothing to copy from, so the reduced data is
written out from the session's own arrays instead. Saving such a session
without an ``output_path`` has no sensible target at all, and is refused with
``NO_SOURCE_FILE``.
"""
from __future__ import annotations

from pathlib import Path
from typing import Optional, Tuple

from pyirena.api._paths import PathSecurityError, resolve_safe
from pyirena.api.control.errors import make_error
from pyirena.api.control.session import Session


def resolve_save_target(
    session: Session, output_path: Optional[str]
) -> Tuple[Optional[Path], Optional[dict]]:
    """Resolve where a ``save_*`` tool should write, creating the file if needed.

    Parameters
    ----------
    session : Session
        The session being saved. ``session.file_path`` may be None.
    output_path : str or None
        Explicit target. None means "write back into the session's own file".

    Returns
    -------
    (target, error)
        Exactly one is not None. ``target`` is a path that exists and holds a
        complete NXcanSAS file, ready for ``save_<tool>_results`` to add a
        results group to. ``error`` is a control-API error dict.
    """
    if session.file_path is None and not output_path:
        return None, make_error(
            "This session was created from data, not from a file, so there is "
            "nothing to save back into.",
            suggestion="Pass output_path to write a new NXcanSAS file.",
            code="NO_SOURCE_FILE",
        )

    # Confine the write target to PYIRENA_DATA_ROOT (when set) for both an
    # explicit output_path and the default in-place save.
    try:
        src = resolve_safe(session.file_path, must_exist=False) if session.file_path else None
        target = resolve_safe(output_path, must_exist=False) if output_path else src
    except PathSecurityError as exc:
        return None, make_error(
            str(exc),
            suggestion="Save to a path inside PYIRENA_DATA_ROOT.",
            code="PATH_NOT_ALLOWED",
        )

    # Saving to a *new* location must yield a complete, re-openable NXcanSAS
    # file — not a results-only stub. With a source file, seed it from there
    # (reduced data + metadata, stale results stripped); the original is never
    # modified. Without one, write the session's arrays out as a fresh file.
    if target != src and not target.exists():
        try:
            if src is not None:
                from pyirena.io._nxcansas_common import copy_and_strip_results
                copy_and_strip_results(src, target)
            else:
                from pyirena.io.nxcansas_unified import create_nxcansas_file
                target.parent.mkdir(parents=True, exist_ok=True)
                create_nxcansas_file(
                    filepath=target,
                    q=session.q,
                    intensity=session.intensity,
                    error=session.error,
                    dq=session.dq,
                    sample_name=session.label or target.stem,
                )
        except Exception as exc:
            return None, make_error(
                f"Could not create output file '{target}': {exc}",
                suggestion="Check the source file exists and the target is writable.",
                code="SAVE_ERROR",
            )

    return target, None
