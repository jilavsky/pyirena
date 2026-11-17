#!/usr/bin/env python
"""Regenerate the documentation screenshots.

Copied from jilavsky-docs/templates/package-docs. Usage (from the repository
root, with the GUI extras installed)::

    python scripts/make_doc_screenshots.py                 # every shot
    python scripts/make_doc_screenshots.py --only <page>-03-<what>
    python scripts/make_doc_screenshots.py --only <page>   # one page
    python scripts/make_doc_screenshots.py --list
    python scripts/make_doc_screenshots.py --keep          # keep the scratch folder

Each window is opened offscreen (``QT_QPA_PLATFORM=offscreen``) in the light
theme, at a fixed size, on a scratch copy of a tracked file. The shot sets
the state, grabs the window with ``QWidget.grab()`` and writes
``<page dir>/img/<page>/<id>.png`` next to the page that shows it.

Your settings are not touched (``_package.isolate_settings``); the scratch
folder is deleted at the end unless ``--keep`` is given. Not part of the
test suite. Fonts come from the operating system: regenerate a page's shots on
one machine. Shot IDs and what each must show are in
``planning/docs-screenshots/SHOTLIST.md``.
"""

from __future__ import annotations

import argparse
import os
import shutil
import sys
import tempfile
import time
from pathlib import Path

# Must happen before anything imports Qt or the package.
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
_SCRATCH = Path(tempfile.mkdtemp(prefix="docshots-"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from doc_shots import _package  # noqa: E402

_package.isolate_settings(_SCRATCH)

import doc_shots  # noqa: E402
from doc_shots._common import Context, save_png, selected  # noqa: E402

SIZE_WARN = 300_000  # bytes; PLAN.md §9


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument(
        "--only",
        nargs="+",
        metavar="ID_OR_PAGE",
        help="shot ids or page names to regenerate (default: all)",
    )
    ap.add_argument("--list", action="store_true", help="list registered shots and exit")
    ap.add_argument("--keep", action="store_true", help="keep the scratch folder (debugging)")
    args = ap.parse_args(argv)
    try:
        return _run(args)
    finally:
        if args.keep:
            print(f"scratch kept: {_SCRATCH}")
        else:
            shutil.rmtree(_SCRATCH, ignore_errors=True)


def _run(args) -> int:
    doc_shots.load_all()
    if args.list:
        for s in selected(None):
            print(f"{s.shot_id:45s} {s.path.relative_to(_package.REPO)}")
        return 0

    shots = selected(args.only)
    app = _package.make_app()
    ctx = Context(app, _SCRATCH)
    failed = []
    t_all = time.monotonic()
    for s in shots:
        t0 = time.monotonic()
        # Fresh settings per shot: windows save state on close, and the next
        # shot must start from the coded defaults again.
        _package.isolate_settings(Path(tempfile.mkdtemp(dir=_SCRATCH, prefix="shot-")))
        try:
            image = s.func(ctx)
            size = save_png(image, s.path)
        except Exception as exc:  # report and continue with the next shot
            import traceback

            traceback.print_exc()
            failed.append(s.shot_id)
            print(f"FAILED {s.shot_id}: {exc}")
            continue
        finally:
            # A widget the shot already closed may still be listed (its deferred
            # delete has not run); closing it again can raise. Do not force the
            # deferred deletes here: that crashed PySide6 in Bernardyn's run.
            # A shot that needs them gone calls _common.flush_deletes itself.
            for w in app.topLevelWidgets():
                try:
                    w.close()
                    w.deleteLater()
                except Exception as exc:  # one bad widget must not end the run
                    print(f"clean-up after {s.shot_id}: {exc!r}")
            ctx.settle(5)
        flag = "  (> 300 KB)" if size > SIZE_WARN else ""
        print(f"{s.shot_id:45s} {size / 1024:6.0f} KB  {time.monotonic() - t0:5.1f} s{flag}")
    print(f"{len(shots) - len(failed)}/{len(shots)} shots in {time.monotonic() - t_all:.0f} s")
    if failed:
        print("failed:", ", ".join(failed))
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
