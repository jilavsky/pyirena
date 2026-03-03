"""
batch_merge.py — Command-line batch merge for two folders of SAS data.

Finds matched file pairs in --folder1 (DS1, lower-Q) and --folder2 (DS2,
higher-Q), then merges each pair using the pyirena Data Merge engine.
A JSON config file sets the overlap Q range and optimisation options; if
omitted, the overlap is auto-detected and all default settings are used.

Usage
-----
Run from any directory::

    python batch_merge.py --folder1 saxs/ --folder2 waxs/ --config merge.json

Arguments
---------
--folder1 DIR   Folder of DS1 files (lower-Q, absolute scale).  Required.
--folder2 DIR   Folder of DS2 files (higher-Q).  Required.
--config JSON   JSON config file with Q range and optimisation settings.
                Produced by the GUI's "Save JSON Config" button.  Optional.
--output DIR    Output folder for merged files.  Default: <folder1>_merged/
                next to folder1.
--ext1 EXT      File extension filter for folder1 (default: .h5).
--ext2 EXT      File extension filter for folder2 (default: same as --ext1).
--match         Match files by shared prefix + trailing number (default).
                Pass --no-match to process files in alphabetical order instead.
--verbose       Print per-file optimisation results (default: on).
--quiet         Suppress per-file output; only print the final summary.

Examples
--------
Merge matched file pairs, auto-detect overlap::

    python batch_merge.py --folder1 saxs/ --folder2 waxs/

Use a saved config file and a specific output folder::

    python batch_merge.py \\
        --folder1 saxs/ --folder2 waxs/ \\
        --config merge_config.json \\
        --output merged/

Process .dat text files in alphabetical order (no name-based matching)::

    python batch_merge.py \\
        --folder1 usaxs/ --folder2 saxs/ \\
        --ext1 .dat --no-match
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _collect_files(folder: Path, ext: str) -> list[Path]:
    """Return sorted list of files with the given extension in *folder*."""
    files = sorted(folder.glob(f"*{ext}"))
    if not files:
        # Try case-insensitive upper extension too (.H5, .DAT, …)
        files = sorted(folder.glob(f"*{ext.upper()}"))
    return files


def _match_pairs(
    files1: list[Path],
    files2: list[Path],
) -> list[tuple[Path, Path]]:
    """
    Match files by (prefix-before-first-underscore, last-integer-in-stem).

    Uses the same algorithm as DataMerge.match_files() but works directly on
    Path objects.  Returns a sorted list of (ds1_path, ds2_path) tuples.
    Unmatched files are silently skipped.
    """
    import re

    def _key(p: Path) -> tuple[str, str]:
        stem = p.stem
        prefix = stem.split('_')[0]
        nums = re.findall(r'\d+', stem)
        last_num = nums[-1] if nums else ''
        return prefix, last_num

    idx2 = {_key(f): f for f in files2}
    pairs: list[tuple[Path, Path]] = []
    seen_keys: set[tuple[str, str]] = set()

    for f1 in files1:
        k = _key(f1)
        if k in idx2 and k not in seen_keys:
            pairs.append((f1, idx2[k]))
            seen_keys.add(k)
        elif k not in idx2:
            print(f"  [warn] No match in folder2 for: {f1.name}")

    return pairs


def _zip_pairs(
    files1: list[Path],
    files2: list[Path],
) -> list[tuple[Path, Path]]:
    """Pair files in list order (first with first, etc.)."""
    n = min(len(files1), len(files2))
    if len(files1) != len(files2):
        print(
            f"  [warn] folder1 has {len(files1)} file(s), "
            f"folder2 has {len(files2)} — processing first {n} pair(s)."
        )
    return list(zip(files1[:n], files2[:n]))


def _print_summary(rows: list[dict]) -> None:
    """Print a fixed-width summary table."""
    if not rows:
        print("\nNo files processed.")
        return

    col_file = max((len(r['file1']) for r in rows), default=8)
    col_file = max(col_file, 8)

    header = (
        f"{'DS1 file':<{col_file}}  "
        f"{'Status':<7}  "
        f"{'Scale':>10}  "
        f"{'BG':>12}  "
        f"{'chi²':>10}  "
        f"Output"
    )
    sep = '-' * len(header)
    print()
    print(sep)
    print(header)
    print(sep)

    for r in rows:
        status = 'OK' if r['success'] else 'FAIL'
        scale  = f"{r['scale']:.4g}"  if r['scale']  is not None else '—'
        bg     = f"{r['bg']:.4g}"     if r['bg']     is not None else '—'
        chi2   = f"{r['chi2']:.4g}"   if r['chi2']   is not None else '—'
        out    = r['output'] or r.get('error', '—')
        print(
            f"{r['file1']:<{col_file}}  "
            f"{status:<7}  "
            f"{scale:>10}  "
            f"{bg:>12}  "
            f"{chi2:>10}  "
            f"{out}"
        )

    print(sep)
    n_ok   = sum(1 for r in rows if r['success'])
    n_fail = len(rows) - n_ok
    print(f"Total: {len(rows)}  |  OK: {n_ok}  |  Failed: {n_fail}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Batch-merge two folders of SAS data using the pyirena Data Merge engine.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )

    parser.add_argument('--folder1', required=True,
                        help='DS1 folder (lower-Q, absolute scale)')
    parser.add_argument('--folder2', required=True,
                        help='DS2 folder (higher-Q)')
    parser.add_argument('--config', default=None,
                        help='JSON merge config file (optional)')
    parser.add_argument('--output', default=None,
                        help='Output folder (default: <folder1>_merged/)')
    parser.add_argument('--ext1', default='.h5',
                        help='File extension for folder1 (default: .h5)')
    parser.add_argument('--ext2', default=None,
                        help='File extension for folder2 (default: same as --ext1)')
    parser.add_argument('--match', dest='match', action='store_true', default=True,
                        help='Match files by name (default)')
    parser.add_argument('--no-match', dest='match', action='store_false',
                        help='Pair files in alphabetical order instead of matching by name')
    parser.add_argument('--quiet', action='store_true', default=False,
                        help='Suppress per-file verbose output')

    args = parser.parse_args()

    folder1 = Path(args.folder1)
    folder2 = Path(args.folder2)

    # Validate folders
    for folder, label in [(folder1, '--folder1'), (folder2, '--folder2')]:
        if not folder.is_dir():
            print(f"Error: {label} '{folder}' is not a directory.", file=sys.stderr)
            sys.exit(1)

    ext1 = args.ext1 if args.ext1.startswith('.') else '.' + args.ext1
    ext2 = (args.ext2 if args.ext2 else args.ext1)
    if not ext2.startswith('.'):
        ext2 = '.' + ext2

    # Collect files
    files1 = _collect_files(folder1, ext1)
    files2 = _collect_files(folder2, ext2)

    if not files1:
        print(f"Error: no *{ext1} files found in {folder1}", file=sys.stderr)
        sys.exit(1)
    if not files2:
        print(f"Error: no *{ext2} files found in {folder2}", file=sys.stderr)
        sys.exit(1)

    # Build pairs
    if args.match:
        print(f"Matching files in '{folder1.name}' with '{folder2.name}' by name…")
        pairs = _match_pairs(files1, files2)
    else:
        print(f"Pairing files in '{folder1.name}' and '{folder2.name}' in order…")
        pairs = _zip_pairs(files1, files2)

    if not pairs:
        print("No file pairs found — nothing to do.", file=sys.stderr)
        sys.exit(1)

    print(f"Found {len(pairs)} pair(s) to merge.\n")

    # Determine output folder
    if args.output:
        output_folder = Path(args.output)
    else:
        output_folder = folder1.parent / (folder1.name + '_merged')

    output_folder.mkdir(parents=True, exist_ok=True)
    print(f"Output folder: {output_folder}\n")

    # Import merge_data (deferred so --help works without pyirena installed)
    try:
        from pyirena.batch import merge_data
    except ImportError as exc:
        print(f"Error: could not import pyirena.batch: {exc}", file=sys.stderr)
        print("Install pyirena with:  pip install -e .", file=sys.stderr)
        sys.exit(1)

    # Process pairs
    summary_rows: list[dict] = []

    for i, (f1, f2) in enumerate(pairs, 1):
        print(f"[{i}/{len(pairs)}] {f1.name}  +  {f2.name}")
        try:
            result = merge_data(
                f1, f2,
                config_file=args.config,
                save_to_nexus=True,
                output_folder=output_folder,
                verbose=(not args.quiet),
            )
        except Exception as exc:
            print(f"  ERROR: {exc}")
            summary_rows.append({
                'file1':   f1.name,
                'success': False,
                'scale':   None,
                'bg':      None,
                'chi2':    None,
                'output':  None,
                'error':   str(exc),
            })
            continue

        if result is None:
            print("  ERROR: merge_data returned None (check messages above)")
            summary_rows.append({
                'file1':   f1.name,
                'success': False,
                'scale':   None,
                'bg':      None,
                'chi2':    None,
                'output':  None,
                'error':   'merge_data returned None',
            })
        else:
            status = 'OK' if result['success'] else 'NOT CONVERGED'
            out_path = result.get('output_path') or '(not saved)'
            print(f"  {status}  scale={result['scale']:.4g}  "
                  f"BG={result['background']:.4g}  "
                  f"chi²={result['chi_squared']:.4g}  -> {out_path}")
            summary_rows.append({
                'file1':   f1.name,
                'success': result['success'],
                'scale':   result['scale'],
                'bg':      result['background'],
                'chi2':    result['chi_squared'],
                'output':  Path(out_path).name if out_path != '(not saved)' else out_path,
            })

    _print_summary(summary_rows)

    n_fail = sum(1 for r in summary_rows if not r['success'])
    sys.exit(1 if n_fail else 0)


if __name__ == '__main__':
    main()
