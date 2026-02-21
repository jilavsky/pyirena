"""
plot_saxs — command-line and importable multi-file SAXS plotting utility.

Creates log-log I(Q) plots (with optional size distribution panels) from
NXcanSAS HDF5 files and/or text data files, saving each graph as a JPEG.

Usage
-----
::

    from pyirena import plot_saxs

    # Single file — raw data only
    plot_saxs("sample.h5")

    # Single file — data + unified fit overlay
    plot_saxs("sample.h5", content=["data", "unified_fit"])

    # Folder of files — all content types, 8 files per graph
    plot_saxs(
        "data/",
        content=["data", "unified_fit", "size_distribution"],
        max_per_graph=8,
        output_dir="figures/",
        dpi=200,
    )

    # List of specific files
    plot_saxs(
        ["run001.h5", "run002.h5", "run003.h5"],
        content=["data", "size_distribution"],
    )
"""

from __future__ import annotations

import math
from pathlib import Path
from typing import List, Optional, Sequence, Tuple, Union


# ---------------------------------------------------------------------------
# Supported file globs (for folder input)
# ---------------------------------------------------------------------------
_GLOB_PATTERNS = ('*.h5', '*.hdf5', '*.hdf', '*.txt', '*.dat')

# Supported content keys
_VALID_CONTENT = frozenset({'data', 'unified_fit', 'size_distribution'})


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def plot_saxs(
    files: Union[str, Path, Sequence[Union[str, Path]]],
    content: Sequence[str] = ('data',),
    output_dir: Optional[Union[str, Path]] = None,
    max_per_graph: Optional[int] = None,
    dpi: int = 150,
    show: bool = False,
) -> List[Path]:
    """Create log-log I(Q) graphs from SAS data files, saved as JPEG images.

    Parameters
    ----------
    files : str, Path, list of str/Path, or folder Path
        Input file(s).  If a folder is given every ``*.h5``, ``*.hdf5``,
        ``*.hdf``, ``*.txt``, and ``*.dat`` file inside is collected
        (non-recursive, sorted by name).
    content : sequence of str
        What to overlay on each graph.  Any combination of:

        ``'data'``
            Experimental I(Q) scatter plot (works with all file types).
        ``'unified_fit'``
            Model I(Q) from Unified Fit results (HDF5 only).  Plotted as a
            line on the same I(Q) axes.
        ``'size_distribution'``
            P(r) size distribution from Size Distribution fit results
            (HDF5 only).  Shown in a second panel beside the I(Q) graph.

        Default: ``('data',)``
    output_dir : str or Path, optional
        Directory where JPEG images are written.  Created if it does not
        exist.  Defaults to the parent folder of the first input file.
    max_per_graph : int, optional
        Maximum number of files placed on a single graph.  When there are
        more files than this limit they are split into multiple JPEG images.
        ``None`` (default) puts all files on one graph.
    dpi : int
        JPEG resolution in dots per inch.  Default 150.
    show : bool
        If ``True`` call ``plt.show()`` to display each graph interactively.
        Default ``False``.

    Returns
    -------
    list of Path
        Absolute paths to the created JPEG files (one per graph chunk).

    Raises
    ------
    ImportError
        If matplotlib is not installed.
    ValueError
        If *content* contains unrecognised keys or no valid files are found.

    Examples
    --------
    Batch-process a folder and export one figure per 10 files::

        from pyirena import plot_saxs
        plot_saxs("data/", content=["data", "unified_fit"], max_per_graph=10)

    Check which JPEGs were created::

        paths = plot_saxs("sample.h5", content=["data", "unified_fit"])
        print(paths)
    """
    try:
        import matplotlib
        matplotlib.use('Agg')   # non-interactive backend for headless use
        import matplotlib.pyplot as plt
        import matplotlib.colors as mcolors
    except ImportError as exc:
        raise ImportError(
            "matplotlib is required for plot_saxs.  Install it with:\n"
            "  pip install matplotlib"
        ) from exc

    import numpy as np

    # ── Validate content ────────────────────────────────────────────────────
    content = list(content)
    unknown = set(content) - _VALID_CONTENT
    if unknown:
        raise ValueError(
            f"Unknown content key(s): {unknown}.  "
            f"Valid values: {sorted(_VALID_CONTENT)}"
        )
    if not content:
        raise ValueError("content must contain at least one item.")

    show_iq   = 'data' in content or 'unified_fit' in content
    show_dist = 'size_distribution' in content

    # ── Resolve file list ───────────────────────────────────────────────────
    file_paths = _resolve_files(files)
    if not file_paths:
        raise ValueError("No valid input files found.")

    # ── Output directory ────────────────────────────────────────────────────
    if output_dir is None:
        out_dir = file_paths[0].parent
    else:
        out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # ── Split into chunks ───────────────────────────────────────────────────
    if max_per_graph is not None and max_per_graph > 0:
        chunks = [
            file_paths[i: i + max_per_graph]
            for i in range(0, len(file_paths), max_per_graph)
        ]
    else:
        chunks = [file_paths]

    output_paths: List[Path] = []

    for chunk_idx, chunk in enumerate(chunks):
        n = len(chunk)
        colors = _gen_colors_mpl(n)

        # ── Build figure layout ─────────────────────────────────────────────
        if show_iq and show_dist:
            fig, (ax_iq, ax_dist) = plt.subplots(
                1, 2, figsize=(14, 6),
                gridspec_kw={'width_ratios': [1.4, 1]},
            )
        elif show_iq:
            fig, ax_iq = plt.subplots(1, 1, figsize=(8, 6))
            ax_dist = None
        else:
            # size_distribution only
            fig, ax_dist = plt.subplots(1, 1, figsize=(7, 6))
            ax_iq = None

        fig.subplots_adjust(hspace=0.35, wspace=0.35)

        # ── Load and plot each file ─────────────────────────────────────────
        for idx, fpath in enumerate(chunk):
            color = colors[idx]
            label = fpath.name

            # -- Raw data (I(Q)) --
            if 'data' in content and ax_iq is not None:
                raw = _load_raw(fpath)
                if raw is not None:
                    q, I = raw
                    ax_iq.scatter(
                        q, I,
                        s=4, color=color,
                        label=label, alpha=0.85,
                        linewidths=0,
                    )

            # -- Unified Fit model I(Q) --
            if 'unified_fit' in content and ax_iq is not None:
                uf = _load_unified_fit(fpath)
                if uf is not None:
                    q_m, I_m = uf
                    fit_lbl = label if 'data' not in content else None
                    ax_iq.plot(
                        q_m, I_m,
                        color=color, linewidth=1.8, alpha=0.9,
                        label=fit_lbl, linestyle='--',
                    )

            # -- Size distribution P(r) --
            if 'size_distribution' in content and ax_dist is not None:
                sd = _load_size_dist(fpath)
                if sd is not None:
                    r, dist = sd
                    ax_dist.plot(
                        r, dist,
                        color=color, linewidth=1.6, alpha=0.9,
                        label=label,
                    )

        # ── Axes formatting ─────────────────────────────────────────────────
        if ax_iq is not None:
            ax_iq.set_xscale('log')
            ax_iq.set_yscale('log')
            ax_iq.set_xlabel('Q (Å⁻¹)', fontsize=12)
            ax_iq.set_ylabel('Intensity (cm⁻¹)', fontsize=12)
            ax_iq.set_title('I(Q)', fontsize=13)
            ax_iq.grid(True, which='both', alpha=0.3, linewidth=0.5)
            _add_legend(ax_iq)

        if ax_dist is not None:
            ax_dist.set_xscale('log')
            ax_dist.set_xlabel('r (Å)', fontsize=12)
            ax_dist.set_ylabel('P(r)  (vol. frac. Å⁻¹)', fontsize=12)
            ax_dist.set_title('Size Distribution', fontsize=13)
            ax_dist.grid(True, which='both', alpha=0.3, linewidth=0.5)
            _add_legend(ax_dist)

        # ── Overall title ───────────────────────────────────────────────────
        suffix = f' ({chunk_idx + 1}/{len(chunks)})' if len(chunks) > 1 else ''
        content_str = ' + '.join(content)
        fig.suptitle(f'pyIrena — {content_str}{suffix}', fontsize=11, y=1.01)

        # ── Save ────────────────────────────────────────────────────────────
        stem = _make_stem(chunk, chunk_idx, len(chunks))
        out_path = out_dir / f'{stem}.jpg'
        fig.savefig(
            str(out_path),
            dpi=dpi,
            format='jpeg',
            bbox_inches='tight',
        )
        output_paths.append(out_path)
        print(f"[plot_saxs] Saved: {out_path}")

        if show:
            plt.show()

        plt.close(fig)

    return output_paths


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _resolve_files(
    files: Union[str, Path, Sequence[Union[str, Path]]]
) -> List[Path]:
    """Expand input into a sorted list of Path objects."""
    if isinstance(files, (str, Path)):
        p = Path(files)
        if p.is_dir():
            result: List[Path] = []
            for pat in _GLOB_PATTERNS:
                result.extend(p.glob(pat))
            return sorted(set(result))
        return [p]
    # Sequence
    return sorted({Path(f) for f in files})


def _load_raw(fpath: Path):
    """Return (Q, I) arrays from a data file, or None on failure."""
    try:
        ext = fpath.suffix.lower()
        if ext in ('.txt', '.dat'):
            from pyirena.io.hdf5 import readTextFile
            data = readTextFile(str(fpath.parent), fpath.name)
        else:
            from pyirena.io.hdf5 import readGenericNXcanSAS
            data = readGenericNXcanSAS(str(fpath.parent), fpath.name)
        if data is None:
            return None
        import numpy as np
        q = np.asarray(data['Q'], dtype=float)
        I = np.asarray(data['Intensity'], dtype=float)
        mask = (q > 0) & (I > 0)
        return q[mask], I[mask]
    except Exception:
        return None


def _load_unified_fit(fpath: Path):
    """Return (Q, I_model) arrays from HDF5 Unified Fit results, or None."""
    from pyirena.io.results import load_result
    try:
        r = load_result(fpath, 'unified_fit')
        if not r['found']:
            return None
        import numpy as np
        q = np.asarray(r['Q'], dtype=float)
        I = np.asarray(r['intensity_model'], dtype=float)
        mask = (q > 0) & (I > 0)
        return q[mask], I[mask]
    except Exception:
        return None


def _load_size_dist(fpath: Path):
    """Return (r_grid, distribution) arrays from HDF5 size dist results, or None."""
    from pyirena.io.results import load_result
    try:
        r = load_result(fpath, 'size_distribution')
        if not r['found'] or r['r_grid'] is None or r['distribution'] is None:
            return None
        import numpy as np
        r_arr = np.asarray(r['r_grid'], dtype=float)
        d_arr = np.asarray(r['distribution'], dtype=float)
        return r_arr, d_arr
    except Exception:
        return None


def _gen_colors_mpl(n: int) -> list:
    """Return *n* evenly-spaced colours from blue to red (matplotlib format)."""
    import matplotlib.pyplot as plt
    import numpy as np
    if n <= 0:
        return []
    if n == 1:
        return ['#3070b0']  # single medium blue
    cmap = plt.get_cmap('gist_rainbow_r')
    return [cmap(i / (n - 1)) for i in range(n)]


def _add_legend(ax) -> None:
    """Add a legend only when there are labelled artists."""
    handles, labels = ax.get_legend_handles_labels()
    if handles:
        ax.legend(
            handles, labels,
            fontsize=8,
            loc='best',
            framealpha=0.8,
            markerscale=2.0,
        )


def _make_stem(chunk: List[Path], chunk_idx: int, n_chunks: int) -> str:
    """Create a JPEG filename stem from the chunk."""
    if n_chunks == 1 and len(chunk) == 1:
        return chunk[0].stem
    if n_chunks == 1:
        # Multiple files, one graph: use first stem + count
        return f'{chunk[0].stem}_and_{len(chunk) - 1}_more'
    return f'plot_saxs_part{chunk_idx + 1:03d}'
