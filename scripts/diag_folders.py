"""Diagnose which USAXS/SAXS/WAXS sample folders the pxp importer sees,
and which of those have the expected DSM_*/R_* wave triple.

Usage:
    python scripts/diag_folders.py /path/to/Big.pxp [--technique USAXS]
"""
from __future__ import annotations

import argparse
from pathlib import Path

from pyirena.io.pxp_to_nexus import (
    _load_pxp_filesystem, _walk_tree, _pick_waves,
    WAVE_PICKERS, TECHNIQUE_FOLDERS,
)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("pxp", help="Path to the .pxp file")
    ap.add_argument("--technique", choices=["USAXS", "SAXS", "WAXS"], default=None,
                    help="Only show folders for this technique (default: all)")
    ap.add_argument("--show-waves", action="store_true",
                    help="For each folder, list every wave name it contains")
    args = ap.parse_args()

    pxp = Path(args.pxp)
    print(f"Loading {pxp} ({pxp.stat().st_size:,} bytes)\n")

    filesystem, n_skipped = _load_pxp_filesystem(pxp)
    print(f"Loaded; {n_skipped} record(s) unparseable\n")

    n_with_waves = 0
    n_missing_triple = 0
    n_total = 0

    for rel_path, technique, folder in _walk_tree(filesystem["root"], None, ()):
        if args.technique and technique != args.technique:
            continue
        n_total += 1
        igor_path = "root:" + ":".join(rel_path)

        picked = _pick_waves(folder, technique)
        wave_names = sorted(k for k, v in folder.items() if not isinstance(v, dict))

        if picked is not None:
            n_with_waves += 1
            status = "OK "
            extra = f"{len(picked[0])} pts"
        else:
            n_missing_triple += 1
            # Inspect the picker requirements to show what's missing
            req = WAVE_PICKERS.get(technique, [])
            missing = []
            for q, i, e, dq in req:
                have = [w for w in (q, i, e) if w in folder]
                if len(have) == 3:
                    # All 3 present in folder dict, but _pick_waves still
                    # returned None — must be a sub-folder collision
                    extra = f"all 3 present but _pick_waves failed (sub-folder collision?)"
                    break
                missing.append(f"{q}/{i}/{e}")
            else:
                extra = f"missing wave triple: tried {missing}"
            status = "-- "

        print(f"  {status} {technique:6s} {igor_path:70s} {extra}")
        if args.show_waves:
            print(f"             waves ({len(wave_names)}): {wave_names}")

    print(f"\nTotals: {n_with_waves} importable, {n_missing_triple} missing triple, "
          f"{n_total} total sample folders (in selected techniques)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
