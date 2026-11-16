#!/usr/bin/env python
"""Regenerate the documentation screenshots in docs/tools/img/.

Usage (from the repository root, with the GUI extras installed)::

    python scripts/make_doc_screenshots.py                 # every shot
    python scripts/make_doc_screenshots.py --only unified-fit-03-guinier-fit
    python scripts/make_doc_screenshots.py --only unified-fit   # one page
    python scripts/make_doc_screenshots.py --list

Each panel is opened offscreen (``QT_QPA_PLATFORM=offscreen``) with pyIrena's
light theme, at a fixed window size, on a copy of a tracked dataset from
``validationData/`` (``testData/`` where a tool needs real data). The script
sets the state (runs the fit, selects a tab), grabs the window with
``QWidget.grab()`` and writes ``docs/tools/img/<page>/<id>.png``. Numbered
callouts on "The window" shots are drawn from the widgets' own geometry.

Your settings are not touched: ``PYIRENA_STATE_DIR`` points at a temporary
folder and window-geometry persistence is off. Not part of the test suite.
Fonts come from the operating system, so shots made on macOS and on Windows
differ slightly; regenerate a page's shots on one machine.

The shot ids, and what each shot must show, are listed in
``planning/docs-screenshots/SHOTLIST.md``.
"""

from __future__ import annotations

import argparse
import os
import sys
import tempfile
import time
from pathlib import Path

# Must happen before anything imports Qt or pyIrena's state manager.
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
_SCRATCH = Path(tempfile.mkdtemp(prefix="pyirena-docshots-"))
os.environ["PYIRENA_STATE_DIR"] = str(_SCRATCH / "state")

sys.path.insert(0, str(Path(__file__).resolve().parent))

import doc_shots  # noqa: E402
from doc_shots._common import Context, make_app, save_png, selected  # noqa: E402

SIZE_WARN = 300_000  # bytes; PLAN.md §9


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--only", nargs="+", metavar="ID_OR_PAGE",
                    help="shot ids or page names to regenerate (default: all)")
    ap.add_argument("--list", action="store_true", help="list registered shots and exit")
    args = ap.parse_args(argv)

    doc_shots.load_all()
    if args.list:
        for s in selected(None):
            print(f"{s.shot_id:45s} {s.path}")
        return 0

    shots = selected(args.only)
    app = make_app()
    ctx = Context(app, _SCRATCH)
    failed = []
    t_all = time.monotonic()
    for s in shots:
        t0 = time.monotonic()
        # A fresh settings folder per shot: panels save their state on close,
        # and the next shot must start from the coded defaults again.
        os.environ["PYIRENA_STATE_DIR"] = tempfile.mkdtemp(dir=_SCRATCH, prefix="state-")
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
            for w in app.topLevelWidgets():
                w.close()
                w.deleteLater()
            ctx.settle(5)
        flag = "  (> 300 KB)" if size > SIZE_WARN else ""
        print(f"{s.shot_id:45s} {size / 1024:6.0f} KB  {time.monotonic() - t0:5.1f} s{flag}")
    print(f"{len(shots) - len(failed)}/{len(shots)} shots in {time.monotonic() - t_all:.0f} s")
    if failed:
        print("failed:", ", ".join(failed))
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
