"""
_nxcansas_common.py — shared NXcanSAS output helpers.

Helpers used by both the Data Merge and Data Manipulation I/O modules
(``nxcansas_data_merge`` / ``nxcansas_data_manipulation``) when writing
processed 1-D data back to an NXcanSAS file:

- stripping non-positive / non-finite data points,
- copying an input file while removing existing pyirena result groups,
- replacing the Q/I/Idev/Qdev arrays in-place,
- appending a Qdev dataset to a freshly created file.
"""
from __future__ import annotations

import logging
import shutil
from datetime import datetime
from pathlib import Path
from typing import Optional

import h5py
import numpy as np

from pyirena.io.hdf5 import find_matching_groups
from pyirena.io.schema import TOOL_REGISTRY

log = logging.getLogger(__name__)

# Every pyirena result group, stripped when a tool seeds a new output file from
# a source file.  Derived from TOOL_REGISTRY rather than listed by hand: a
# hand-maintained copy silently fell two tools behind (SAXS Morph and Fractals
# results were being carried into derived data files), which is what
# pyirena/tests/test_tool_registration.py now guards.
PYIRENA_RESULT_GROUPS = [schema['group'] for schema in TOOL_REGISTRY.values()]


def strip_nonpositive_intensities(
    q: np.ndarray,
    I: np.ndarray,
    dI: np.ndarray,
    dQ: Optional[np.ndarray],
) -> tuple[np.ndarray, np.ndarray, np.ndarray, Optional[np.ndarray], int]:
    """Drop any points with non-positive or non-finite Q or I.

    Subtraction and merge operations can produce negative intensities (e.g.
    at sample-buffer boundaries, or at the edges of DS1 where a flat
    background is subtracted).  Such points are not physically meaningful
    and break downstream analysis (in particular log plots and the Size
    Distribution fit).  Igor Pro's equivalent routines delete these points
    before saving; we do the same.

    Returns the cleaned arrays and the number of points removed.
    """
    q  = np.asarray(q,  dtype=float)
    I  = np.asarray(I,  dtype=float)
    dI = np.asarray(dI, dtype=float)
    mask = np.isfinite(q) & np.isfinite(I) & (q > 0) & (I > 0)
    n_removed = int(np.sum(~mask))
    if n_removed == 0:
        return q, I, dI, dQ, 0
    q  = q[mask]
    I  = I[mask]
    dI = dI[mask]
    if dQ is not None:
        dQ = np.asarray(dQ, dtype=float)[mask]
    return q, I, dI, dQ, n_removed


def copy_and_strip_results(src: Path, dst: Path) -> None:
    """Copy *src* to *dst* then delete known pyirena result groups."""
    shutil.copy2(src, dst)
    with h5py.File(dst, 'a') as f:
        for grp_path in PYIRENA_RESULT_GROUPS:
            if grp_path in f:
                del f[grp_path]


def replace_nxcansas_data(
    filepath: Path,
    q: np.ndarray,
    I: np.ndarray,
    dI: np.ndarray,
    dQ: Optional[np.ndarray],
    context: str = "data",
) -> None:
    """Overwrite Q/I/Idev/Qdev arrays in the first NXcanSAS sasdata group.

    *context* is used only in the error message (e.g. "merged data",
    "manipulated data") so callers can identify which tool failed.
    """
    with h5py.File(filepath, 'a') as f:
        # Locate the sasdata group dynamically (same approach as the reader)
        sasdata_paths = find_matching_groups(
            f,
            required_attributes={'canSAS_class': 'SASdata'},
            required_items={},
        )
        if not sasdata_paths:
            # Fallback: look for NXdata group
            sasdata_paths = find_matching_groups(
                f,
                required_attributes={'NX_class': 'NXdata'},
                required_items={},
            )
        if not sasdata_paths:
            raise RuntimeError(
                f"Could not locate a sasdata/NXdata group in {filepath}. "
                f"Cannot replace {context}."
            )

        sasdata = f[sasdata_paths[0]]

        # Replace each array dataset; delete-then-recreate to allow size change
        for name, data, attrs in [
            ('Q',    q,  {'units': '1/angstrom', 'long_name': 'Q'}),
            ('I',    I,  {'units': '1/cm',       'long_name': 'Intensity'}),
            ('Idev', dI, {'units': '1/cm',       'long_name': 'Uncertainties'}),
        ]:
            if name in sasdata:
                del sasdata[name]
            ds = sasdata.create_dataset(name, data=data)
            for k, v in attrs.items():
                ds.attrs[k] = v

        # I.uncertainties attribute
        if 'I' in sasdata:
            sasdata['I'].attrs['uncertainties'] = 'Idev'

        # A stale scalar slit length (dQl) from the source file no longer
        # matches the rewritten Q/I (point count may have changed).  Drop it
        # here; callers that produce slit-smeared output re-add it explicitly
        # via append_dql() after this call.
        if 'dQl' in sasdata:
            del sasdata['dQl']

        # Qdev (Q resolution) — add only if provided
        if dQ is not None:
            if 'Qdev' in sasdata:
                del sasdata['Qdev']
            ds_qdev = sasdata.create_dataset('Qdev', data=dQ)
            ds_qdev.attrs['units'] = '1/angstrom'
            ds_qdev.attrs['long_name'] = 'Q resolution'
            sasdata['Q'].attrs['resolutions'] = 'Qdev'
        else:
            # Remove stale Qdev if the new data has none
            if 'Qdev' in sasdata:
                del sasdata['Qdev']
            if 'Q' in sasdata and 'resolutions' in sasdata['Q'].attrs:
                del sasdata['Q'].attrs['resolutions']

        # Update file timestamp
        f.attrs['file_time'] = datetime.now().isoformat()


def append_dq(filepath: Path, dQ: np.ndarray, sample_name: str) -> None:
    """Add a Qdev dataset to the sasdata group of a freshly-created NXcanSAS file."""
    with h5py.File(filepath, 'a') as f:
        sasdata_paths = find_matching_groups(
            f,
            required_attributes={'canSAS_class': 'SASdata'},
            required_items={},
        )
        if sasdata_paths:
            sasdata = f[sasdata_paths[0]]
            if 'Qdev' not in sasdata:
                ds = sasdata.create_dataset('Qdev', data=dQ)
                ds.attrs['units'] = '1/angstrom'
                ds.attrs['long_name'] = 'Q resolution'
                if 'Q' in sasdata:
                    sasdata['Q'].attrs['resolutions'] = 'Qdev'


def append_dql(filepath: Path, slit_length: float) -> None:
    """Mark the first sasdata group as slit smeared by writing a scalar ``dQl``.

    Adds the slit (half-)length as a ``dQl`` dataset and includes ``dQl`` in
    the ``Q@resolutions`` attribute (SasView / NXcanSAS convention, alongside a
    per-point ``dQw``/``Qdev`` if present).  This is what downstream pyirena
    readers use to auto-detect slit smearing (see ``readGenericNXcanSAS``).
    No-op when ``slit_length <= 0``.
    """
    if not slit_length or slit_length <= 0:
        return
    with h5py.File(filepath, 'a') as f:
        sasdata_paths = find_matching_groups(
            f, required_attributes={'canSAS_class': 'SASdata'}, required_items={},
        )
        if not sasdata_paths:
            return
        sasdata = f[sasdata_paths[0]]
        if 'dQl' in sasdata:
            del sasdata['dQl']
        ds = sasdata.create_dataset('dQl', data=float(slit_length))
        ds.attrs['units'] = '1/angstrom'
        ds.attrs['long_name'] = 'Slit length (half-height)'
        if 'Q' in sasdata:
            existing = sasdata['Q'].attrs.get('resolutions', '')
            tokens = [t.strip() for t in str(existing).split(',') if t.strip()]
            # Per-point width token first (dQw preferred; keep Qdev if that's
            # what is present), then the slit length dQl.
            if 'Qdev' in sasdata and 'Qdev' not in tokens and 'dQw' not in tokens:
                tokens = ['Qdev'] + tokens
            if 'dQl' not in tokens:
                tokens.append('dQl')
            sasdata['Q'].attrs['resolutions'] = ','.join(tokens)


def drop_smr_entries(filepath: Path) -> int:
    """Delete any ``*_SMR`` (slit-smeared twin) entry groups from *filepath*.

    Matilda writes both a desmeared default entry and a ``<name>_SMR`` slit-
    smeared copy.  When a tool rewrites only the default entry (merge /
    manipulation), the ``_SMR`` twin would otherwise survive with stale,
    now-inconsistent data.  Removing it prevents a later
    ``prefer_slit_smeared`` load from silently returning the wrong curve.

    Returns the number of groups removed.
    """
    removed = 0
    with h5py.File(filepath, 'a') as f:
        entry = f.get('entry')
        # Top-level entries may also be siblings of 'entry' in some layouts;
        # scan both the root and the 'entry' group for *_SMR children.
        scopes = []
        if isinstance(entry, h5py.Group):
            scopes.append(entry)
        scopes.append(f)
        for scope in scopes:
            for key in list(scope.keys()):
                if key.endswith('_SMR') and isinstance(scope.get(key), h5py.Group):
                    del scope[key]
                    removed += 1
    return removed


# ---------------------------------------------------------------------------
# Source metadata propagation (Data Merge)
# ---------------------------------------------------------------------------
#
# A merged file is seeded from DS1, so DS1's metadata survives by being copied
# with the rest of the file.  DS2 contributed only numbers, which meant that
# merging USAXS + SAXS threw away everything recorded about the SAXS
# measurement.  The groups below are carried across under a technique suffix
# so a USAXS + SAXS + WAXS chain ends up with ``metadata`` (USAXS, from DS1),
# ``metadata_saxs`` and ``metadata_waxs`` side by side.  See GitHub issue #21.

#: Groups carried from a merge input into the output.
PROPAGATED_METADATA_GROUPS = ('metadata', 'instrument', 'sample')

#: Datasets at or above this many bytes are skipped when copying those groups.
#: SAXS/WAXS files keep the raw 2-D detector image under
#: ``entry/instrument/detector/data`` — ~10 MB on a current Pilatus/Eiger, and
#: meaningless in a merged 1-D curve file.  Matilda drops the same dataset when
#: it reduces (``convertSWAXS._geometry_from_dicts`` callers do
#: ``del instrument_dict['detector']['data']``).
BULK_DATASET_BYTES = 1_000_000


def _child_ci(group: h5py.Group, name: str) -> Optional[str]:
    """Return *group*'s child whose name matches *name* ignoring case.

    Beamline files are not consistent about this: USAXS reductions write
    ``entry/metadata`` while SAXS/WAXS files carry the raw ``entry/Metadata``.
    An exact-match lookup silently finds nothing on half the inputs.
    """
    lowered = name.lower()
    for key in group:
        if key.lower() == lowered:
            return key
    return None


def _nxentry(f: h5py.File) -> Optional[h5py.Group]:
    """Return the file's NXentry group, or None."""
    key = _child_ci(f, 'entry')
    if key is not None and isinstance(f[key], h5py.Group):
        return f[key]
    for key in f:
        obj = f[key]
        if isinstance(obj, h5py.Group):
            nx = obj.attrs.get('NX_class')
            if isinstance(nx, bytes):
                nx = nx.decode('utf-8', 'replace')
            if nx == 'NXentry':
                return obj
    return None


def detect_technique(filepath: Path) -> str:
    """Identify which instrument produced *filepath*.

    Returns ``'usaxs'``, ``'saxs'``, ``'waxs'``, or ``'unknown'`` when the file
    is not HDF5, cannot be opened, or carries none of the markers.

    The SAXS/WAXS discriminator is the one Matilda already relies on to pick a
    detector geometry (``convertSWAXS._geometry_from_dicts``): a SAXS file's
    metadata has ``pin_ccd_tilt_x``, a WAXS file's has ``waxs_ccd_tilt_x``.
    A USAXS reduction has neither and has a ``flyScan`` group.  Nothing is
    inferred from the filename — a renamed file would then mislabel its own
    metadata, which is worse than labelling it ``unknown``.
    """
    filepath = Path(filepath)
    try:
        if not h5py.is_hdf5(filepath):
            return 'unknown'
    except Exception:
        return 'unknown'

    try:
        with h5py.File(filepath, 'r') as f:
            entry = _nxentry(f)
            if entry is None:
                return 'unknown'

            md_key = _child_ci(entry, 'metadata')
            if md_key is not None and isinstance(entry[md_key], h5py.Group):
                md = entry[md_key]
                if _child_ci(md, 'pin_ccd_tilt_x') is not None:
                    return 'saxs'
                if _child_ci(md, 'waxs_ccd_tilt_x') is not None:
                    return 'waxs'

            if _child_ci(entry, 'flyScan') is not None:
                return 'usaxs'
    except Exception:
        log.debug("technique detection failed for %s", filepath, exc_info=True)
        return 'unknown'

    return 'unknown'


def _copy_pruned(src: h5py.Group, dst_parent: h5py.Group, dst_name: str,
                 skipped: list) -> h5py.Group:
    """Recursively copy *src* to ``dst_parent[dst_name]``, skipping bulk arrays.

    Datasets of ``BULK_DATASET_BYTES`` or more are not copied; their paths are
    appended to *skipped* so the destination group can record what was left
    behind rather than quietly presenting a partial instrument group as whole.
    """
    dst = dst_parent.create_group(dst_name)
    for key, value in src.attrs.items():
        dst.attrs[key] = value

    def _recurse(s: h5py.Group, d: h5py.Group) -> None:
        for key in s:
            try:
                obj = s[key]
            except Exception:
                log.debug("unreadable item %s/%s — skipped", s.name, key, exc_info=True)
                continue
            if isinstance(obj, h5py.Group):
                sub = d.create_group(key)
                for akey, avalue in obj.attrs.items():
                    sub.attrs[akey] = avalue
                _recurse(obj, sub)
            elif isinstance(obj, h5py.Dataset):
                if obj.nbytes >= BULK_DATASET_BYTES:
                    skipped.append(obj.name)
                    continue
                try:
                    s.copy(key, d, name=key)
                except Exception:
                    log.debug("could not copy %s", obj.name, exc_info=True)

    _recurse(src, dst)
    return dst


def copy_metadata_groups(
    src_path: Path,
    dst_path: Path,
    suffix: Optional[str] = None,
    groups: tuple = PROPAGATED_METADATA_GROUPS,
) -> dict:
    """Carry *src_path*'s metadata/instrument/sample groups into *dst_path*.

    Each group is written under the NXentry of *dst_path* as
    ``<name>_<suffix>`` — ``metadata_saxs``, ``instrument_saxs``, … — or under
    its plain ``<name>`` when *suffix* is None, which is what a fresh output
    file with no DS1 groups of its own wants.  An existing destination group of
    the same name is replaced, so re-saving a merge is idempotent.

    Datasets of ``BULK_DATASET_BYTES`` or more are skipped (see
    :func:`_copy_pruned`); each destination group records the source file in
    ``pyirena_source_file`` and anything left out in ``pyirena_skipped``.

    Returns ``{'copied': [names], 'skipped': [source paths]}``.  A source that
    is not HDF5, has no NXentry, or has none of *groups* is not an error — the
    result simply reports nothing copied.
    """
    src_path = Path(src_path)
    dst_path = Path(dst_path)
    result = {'copied': [], 'skipped': []}

    try:
        if not h5py.is_hdf5(src_path):
            return result
    except Exception:
        return result

    try:
        with h5py.File(src_path, 'r') as src, h5py.File(dst_path, 'a') as dst:
            src_entry = _nxentry(src)
            dst_entry = _nxentry(dst)
            if src_entry is None or dst_entry is None:
                return result

            for name in groups:
                src_key = _child_ci(src_entry, name)
                if src_key is None or not isinstance(src_entry[src_key], h5py.Group):
                    continue

                dst_name = f"{name}_{suffix}" if suffix else name

                # Only ever replace a group we wrote ourselves on an earlier
                # save (that is what makes re-saving idempotent).  Any other
                # occupant is somebody else's data and must not be deleted:
                # create_nxcansas_file() names the *subentry* after the sample,
                # so a file called sample.h5 has its merged sasdata sitting at
                # entry/sample — exactly where DS1's NXsample group wants to
                # go.  Deleting that would throw away the merged curve.
                existing = dst_entry.get(dst_name) if dst_name in dst_entry else None
                if existing is not None:
                    if (isinstance(existing, h5py.Group)
                            and 'pyirena_source_file' in existing.attrs):
                        del dst_entry[dst_name]
                    else:
                        alt = f"{dst_name}_src"
                        if alt in dst_entry:
                            log.warning(
                                "Not propagating %s from %s: both %s and %s are "
                                "already taken in %s.",
                                name, src_path, dst_name, alt, dst_path,
                            )
                            continue
                        log.warning(
                            "%s already exists in %s and was not written by "
                            "pyirena — propagating %s's %s as %s instead.",
                            dst_name, dst_path, src_path, name, alt,
                        )
                        dst_name = alt

                skipped: list = []
                grp = _copy_pruned(src_entry[src_key], dst_entry, dst_name, skipped)
                grp.attrs['pyirena_source_file'] = str(src_path)
                grp.attrs['pyirena_source_group'] = src_entry[src_key].name
                if skipped:
                    grp.attrs['pyirena_skipped'] = [str(p) for p in skipped]

                result['copied'].append(dst_name)
                result['skipped'].extend(skipped)
    except Exception:
        log.warning("Could not propagate metadata from %s into %s",
                    src_path, dst_path, exc_info=True)

    return result
