"""Print the full folder tree of a .pxp file, with optional filter.

Usage:
    python scripts/diag_tree.py /path/to/Big.pxp
    python scripts/diag_tree.py /path/to/Big.pxp --filter AMC
    python scripts/diag_tree.py /path/to/Big.pxp --max-depth 4
"""
from __future__ import annotations

import argparse
from pathlib import Path

from pyirena.io.pxp_to_nexus import _load_pxp_filesystem


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("pxp", help="Path to the .pxp file")
    ap.add_argument("--filter", default=None,
                    help="Only print folder paths whose name contains this substring "
                         "(case-insensitive). Also prints their parents for context.")
    ap.add_argument("--max-depth", type=int, default=4,
                    help="Maximum tree depth to print (default: 4)")
    args = ap.parse_args()

    pxp = Path(args.pxp)
    print(f"Loading {pxp} ({pxp.stat().st_size:,} bytes)\n")
    filesystem, n_skipped = _load_pxp_filesystem(pxp)
    print(f"Loaded; {n_skipped} record(s) unparseable\n")

    filt = args.filter.lower() if args.filter else None

    def walk(d: dict, path: list[str], depth: int) -> bool:
        """Return True if this subtree contained anything to print."""
        printed = False
        for k in sorted(d.keys()):
            v = d[k]
            if not isinstance(v, dict):
                continue
            new_path = path + [k]
            full = "root:" + ":".join(new_path)
            n_waves = sum(1 for x in v.values() if not isinstance(x, dict))
            n_subs = sum(1 for x in v.values() if isinstance(x, dict))

            this_matches = filt is None or filt in k.lower()
            child_matched = False
            if depth + 1 < args.max_depth:
                # Recurse to know if children matched (so we print parents
                # of matches even if the parent name itself doesn't match).
                child_matched = _walk_check(v, filt, args.max_depth, depth + 1)

            if this_matches or child_matched:
                indent = "  " * depth
                print(f"{indent}{k}/  ({n_subs} subfolders, {n_waves} waves)  [{full}]")
                printed = True
                if depth + 1 < args.max_depth:
                    walk(v, new_path, depth + 1)
        return printed

    walk(filesystem["root"], [], 0)
    return 0


def _walk_check(d: dict, filt: str | None, max_depth: int, depth: int) -> bool:
    """Recursive existence check: any folder name in this subtree matches?"""
    if filt is None:
        return True
    if depth >= max_depth:
        return False
    for k, v in d.items():
        if not isinstance(v, dict):
            continue
        if filt in k.lower():
            return True
        if _walk_check(v, filt, max_depth, depth + 1):
            return True
    return False


if __name__ == "__main__":
    raise SystemExit(main())
