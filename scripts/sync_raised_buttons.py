"""
Keep the vendored copies of ``raised_buttons.py`` in step with pyIrena's.

``pyirena/gui/raised_buttons.py`` is the source of truth for the raised
push-button look shared by Jan's Qt packages. Each package carries its own
byte-identical copy, so it installs and runs on its own. This script checks
or refreshes those copies.

Usage (from anywhere)::

    python scripts/sync_raised_buttons.py            # report status (default)
    python scripts/sync_raised_buttons.py --write    # refresh existing copies
    python scripts/sync_raised_buttons.py --write --add   # also create missing ones

``--root`` points at the folder holding the sibling repos (default: the
parent of this pyIrena checkout). Exit status is 1 when a copy has drifted
and ``--write`` was not given, so the check can run from CI or a hook.

Adding a package: append its relative path to ``TARGETS``.
"""

from __future__ import annotations

import argparse
import re
import shutil
import sys
from pathlib import Path

PYIRENA_ROOT = Path(__file__).resolve().parents[1]
SOURCE = PYIRENA_ROOT / "pyirena" / "gui" / "raised_buttons.py"

#: Where each sibling repo keeps its copy, relative to ``--root``.
#: Put it next to that package's ``_qt.py`` shim.
TARGETS: dict[str, str] = {
    "Matilda": "Matilda/matilda/gui/raised_buttons.py",
    "Bernardyn": "Bernardyn/bernardyn/gui/raised_buttons.py",
    "AIDA": "Aida/src/aida/ui/qt/raised_buttons.py",
    "DataReporter": "DataReporter/src/datareporter/gui/raised_buttons.py",
    "Anisa": "Anisa/anisa/gui/raised_buttons.py",
}

_VERSION = re.compile(r'^__version__\s*=\s*["\']([^"\']+)', re.M)


def _version(path: Path) -> str:
    m = _VERSION.search(path.read_text(encoding="utf-8"))
    return m.group(1) if m else "?"


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[1])
    ap.add_argument("--root", type=Path, default=PYIRENA_ROOT.parent,
                    help="folder containing the sibling repos")
    ap.add_argument("--write", action="store_true", help="refresh drifted copies")
    ap.add_argument("--add", action="store_true",
                    help="with --write, also create copies where none exist yet")
    args = ap.parse_args(argv)

    src_bytes = SOURCE.read_bytes()
    src_ver = _version(SOURCE)
    print(f"source  pyirena  v{src_ver}  {SOURCE}")
    drifted = 0
    for name, rel in TARGETS.items():
        dst = args.root / rel
        repo = args.root / rel.split("/")[0]
        if not repo.is_dir():
            print(f"skip    {name:12s} repo not found at {repo}")
            continue
        if not dst.exists():
            if args.write and args.add and dst.parent.is_dir():
                shutil.copyfile(SOURCE, dst)
                print(f"added   {name:12s} v{src_ver}  {dst}")
            else:
                print(f"absent  {name:12s} (not adopted yet; --write --add to create)")
            continue
        if dst.read_bytes() == src_bytes:
            print(f"ok      {name:12s} v{src_ver}")
            continue
        old = _version(dst)
        if args.write:
            shutil.copyfile(SOURCE, dst)
            print(f"updated {name:12s} v{old} -> v{src_ver}  {dst}")
        else:
            drifted += 1
            print(f"DRIFT   {name:12s} v{old} (source v{src_ver})  {dst}")
    if drifted:
        print(f"\n{drifted} copy/copies differ from pyIrena's. Review, then run with --write.")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
