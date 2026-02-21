"""
plot_sas_tree.py — Batch SAXS plotter for a directory tree.

Traverses the current working directory (or a given root) up to 6 levels deep,
finds all HDF5/H5/HDF5 SAS data files, and generates one JPEG image per
directory that contains matching files.  The image filename is built from the
concatenated directory names relative to the root
(e.g. ``topFolder_subFolder_deepFolder.jpg``).

Usage
-----
Run from the directory containing your data (the JPEG files are saved there)::

    python plot_sas_tree.py

Optional arguments::

    python plot_sas_tree.py \\
        --root /path/to/data \\
        --max-per-graph 8 \\
        --content data unified_fit size_distribution \\
        --dpi 200

Arguments
---------
--root          Root directory to traverse.  Defaults to CWD.
--max-per-graph Maximum files per JPEG image (default: no limit).
--content       Space-separated list of content types to include.
                Choices: data, unified_fit, size_distribution.
                Default: data
--dpi           JPEG resolution in dots per inch.  Default: 150.
--output-dir    Where to write JPEG files.  Defaults to --root.
--depth         Maximum folder depth below root (1–6).  Default: 6.
--show          Display each graph interactively (matplotlib GUI).

Examples
--------
Plot I(Q) data only from CWD::

    python plot_sas_tree.py

Include unified-fit overlay, 8 datasets per graph::

    python plot_sas_tree.py --content data unified_fit --max-per-graph 8

Full content, saving to a separate figures folder::

    python plot_sas_tree.py \\
        --content data unified_fit size_distribution \\
        --output-dir figures/ \\
        --dpi 200
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path


# ---------------------------------------------------------------------------
# Supported HDF5 extensions (text files are skipped — they carry no fit results)
# ---------------------------------------------------------------------------
_HDF5_EXTENSIONS = {'.h5', '.hdf5', '.hdf'}

_VALID_CONTENT = {'data', 'unified_fit', 'size_distribution'}


def _collect_dirs(root: Path, max_depth: int) -> list[Path]:
    """Return all directories under *root* (inclusive) up to *max_depth* levels deep."""
    dirs = []
    for dirpath, dirnames, filenames in os.walk(root):
        current = Path(dirpath)
        depth = len(current.relative_to(root).parts)
        if depth > max_depth:
            dirnames.clear()   # prune the walk — don't go deeper
            continue
        dirs.append(current)
        # Sort sub-dirs for reproducible order
        dirnames.sort()
    return dirs


def _hdf5_files_in(directory: Path) -> list[Path]:
    """Return sorted list of HDF5 SAS files directly inside *directory*."""
    files = [
        p for p in directory.iterdir()
        if p.is_file() and p.suffix.lower() in _HDF5_EXTENSIONS
    ]
    return sorted(files)


def _dir_stem(root: Path, directory: Path) -> str:
    """
    Build an output filename stem from the directory path relative to *root*.

    ``root/topFolder/subFolder`` → ``topFolder_subFolder``
    ``root`` itself             → ``root`` (basename of root)
    """
    rel = directory.relative_to(root)
    parts = list(rel.parts)
    if not parts:
        # The root directory itself
        return root.name or 'root'
    return '_'.join(parts)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description='Batch-generate SAXS JPEG images from a directory tree.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        '--root',
        type=Path,
        default=None,
        help='Root directory to traverse.  Defaults to CWD.',
    )
    parser.add_argument(
        '--max-per-graph',
        type=int,
        default=None,
        metavar='N',
        help='Maximum number of files per JPEG image.  Default: all in one image.',
    )
    parser.add_argument(
        '--content',
        nargs='+',
        default=['data'],
        choices=sorted(_VALID_CONTENT),
        metavar='TYPE',
        help=(
            'Content type(s) to include: data, unified_fit, size_distribution. '
            'Default: data'
        ),
    )
    parser.add_argument(
        '--dpi',
        type=int,
        default=150,
        help='JPEG resolution in dots-per-inch.  Default: 150.',
    )
    parser.add_argument(
        '--output-dir',
        type=Path,
        default=None,
        metavar='DIR',
        help='Directory where JPEG files are written.  Defaults to --root.',
    )
    parser.add_argument(
        '--depth',
        type=int,
        default=6,
        choices=range(1, 7),
        metavar='1-6',
        help='Maximum folder depth below root to search.  Default: 6.',
    )
    parser.add_argument(
        '--show',
        action='store_true',
        help='Display each graph interactively after saving.',
    )

    args = parser.parse_args(argv)

    # ── Resolve root ─────────────────────────────────────────────────────────
    root = (args.root or Path.cwd()).resolve()
    if not root.is_dir():
        print(f"ERROR: root directory not found: {root}", file=sys.stderr)
        return 1

    # ── Resolve output directory ──────────────────────────────────────────────
    output_dir = (args.output_dir or root).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    # ── Import plot_saxs ──────────────────────────────────────────────────────
    try:
        from pyirena.plotting.plot_saxs import plot_saxs
    except ImportError as exc:
        print(
            f"ERROR: could not import pyirena.  Make sure pyirena is installed:\n"
            f"  pip install pyirena[plotting]\n"
            f"  {exc}",
            file=sys.stderr,
        )
        return 1

    # ── Traverse directory tree ───────────────────────────────────────────────
    dirs = _collect_dirs(root, args.depth)

    total_images = 0
    skipped_dirs = 0

    for directory in dirs:
        files = _hdf5_files_in(directory)
        if not files:
            skipped_dirs += 1
            continue

        stem = _dir_stem(root, directory)

        # When max_per_graph splits into multiple chunks, append _part001 etc.
        # plot_saxs handles chunking internally and returns list of paths.
        # We need a custom output_dir per-directory so stems don't collide.
        # Strategy: use a temp subdir, then rename with our stem.
        import tempfile, shutil

        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            try:
                created = plot_saxs(
                    files=files,
                    content=args.content,
                    output_dir=tmp_path,
                    max_per_graph=args.max_per_graph,
                    dpi=args.dpi,
                    show=args.show,
                )
            except Exception as exc:
                print(f"  [SKIP] {directory}: {exc}", file=sys.stderr)
                skipped_dirs += 1
                continue

            for i, src in enumerate(created):
                if len(created) == 1:
                    dest_name = f'{stem}.jpg'
                else:
                    dest_name = f'{stem}_part{i + 1:03d}.jpg'
                dest = output_dir / dest_name
                shutil.move(str(src), str(dest))
                print(f'[plot_sas_tree] Saved: {dest}')
                total_images += 1

    # ── Summary ───────────────────────────────────────────────────────────────
    print(f'\nDone.  {total_images} image(s) created, {skipped_dirs} director(y/ies) skipped.')
    return 0


if __name__ == '__main__':
    sys.exit(main())
