"""
pyirena/io/pxp_to_nexus.py — Import legacy Igor Pro packed experiments
(``.pxp``) and export their reduced 1-D SAS data as per-sample NXcanSAS
HDF5 files.

This bridges old Igor-based analysis into the pyirena (Python / NeXus)
world. Each ``.pxp`` typically holds many samples organised by technique
(``root:USAXS:<sample>``, ``root:SAXS:<sample>``, ``root:WAXS:<sample>``).
This module walks that tree, recognises the standard wave-name conventions
used by the USAXS reduction pipeline, parses the rich wave-note metadata
(when present), and writes one ``.h5`` per sample using the existing
:func:`pyirena.io.nxcansas_unified.create_nxcansas_file`.

Wave-name conventions recognised
--------------------------------
USAXS (desmeared, primary):
    ``DSM_Qvec``, ``DSM_Int``, ``DSM_Error``  (+ optional ``DSM_dQ``)

SAXS / WAXS (reduced):
    ``R_Qvec``, ``R_Int``, ``R_Error``        (+ optional ``w_<folder>``
    as dQ if present in the same folder)

Samples that contain none of the above (e.g. blanks holding only raw
detector counts) are silently skipped — they appear as ``skipped`` in
the per-file summary.

Folder-name conventions recognised
----------------------------------
Top-level Igor data-folder names that map to a technique:

================ ==========
Folder name      Technique
================ ==========
``USAXS``        USAXS
``SAXS``         SAXS
``WAXS``         WAXS
``Imported SAXS`` SAXS (legacy ASCII-imported)
``SAS``          SAXS (user convention)
================ ==========

The mapping is data-driven (see :data:`TECHNIQUE_FOLDERS`) and easy to
extend. Folders not in the map are walked recursively — any descendant
folder that contains a recognised wave triple is still exported, and its
position in the tree is preserved in the output directory layout. This
means in-situ experiments organised as
``root:SAXS:<heater_run>:<step>:<sample>`` will round-trip into
``out/SAXS/<heater_run>/<step>/<sample>.h5``.

Wave-note parsing
-----------------
Igor wave notes use the convention ``key=value;key=value;``. Some pipelines
(notably the APS USAXS instrument) inject NeXus-style metadata wrapped in
sentinel markers — ``NXSampleStart`` / ``NXSampleEnd``,
``NXInstrumentStart`` / ``NXInstrumentEnd``, ``NXMetadataStart`` /
``NXMetadataEnd``. The parser extracts these into sectioned dicts that
map cleanly to NXcanSAS sample/instrument groups; anything outside the
sections (or in a note without sentinels) is preserved verbatim under
``entry/notes/``.

Typical usage
-------------
Headless::

    from pyirena.io.pxp_to_nexus import extract_pxp_to_nexus

    summaries = extract_pxp_to_nexus("legacy.pxp")
    for s in summaries:
        print(s.output_path, s.n_points, s.technique)

Custom output folder, only USAXS::

    extract_pxp_to_nexus("legacy.pxp",
                         output_root="/tmp/nexus_out",
                         techniques=["USAXS"])
"""

from __future__ import annotations

import logging
import re
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable, Sequence

import numpy as np

from pyirena.io.nxcansas_unified import create_nxcansas_file


__all__ = [
    "extract_pxp_to_nexus",
    "ExtractedFile",
    "ExtractionResult",
    "TECHNIQUE_FOLDERS",
    "WAVE_PICKERS",
]


logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Configuration tables (data-driven, easy to extend)
# ---------------------------------------------------------------------------

#: Map of Igor folder name → canonical technique string. Folder name matching
#: is case-insensitive and exact. Add aliases here when users follow other
#: conventions (e.g. importing ASCII into a folder called ``"Imported SAXS"``).
TECHNIQUE_FOLDERS: dict[str, str] = {
    "USAXS":         "USAXS",
    "SAXS":          "SAXS",
    "WAXS":          "WAXS",
    "Imported SAXS": "SAXS",
    "SAS":           "SAXS",
    "Imported":      "SAXS",
}


#: Wave-name picker per technique. Each entry is a list of ``(Q_name, I_name,
#: Err_name, dQ_name_or_None)`` tuples tried in order. The first tuple whose
#: Q, I, and Err waves all exist in a sample folder wins. dQ is optional —
#: if its name is None or not present, the resulting file will not contain
#: a Qdev dataset.
WAVE_PICKERS: dict[str, list[tuple[str, str, str, str | None]]] = {
    "USAXS": [
        # Desmeared primary: matches the chosen "DSM_* only" policy.
        ("DSM_Qvec", "DSM_Int", "DSM_Error", "DSM_dQ"),
    ],
    "SAXS": [
        # R_* is the unified Q/I/error triple emitted by the SAXS reduction.
        ("R_Qvec", "R_Int", "R_Error", None),
    ],
    "WAXS": [
        ("R_Qvec", "R_Int", "R_Error", None),
    ],
}


# ---------------------------------------------------------------------------
# Result types
# ---------------------------------------------------------------------------

@dataclass
class ExtractedFile:
    """One row in the extraction summary.

    Attributes
    ----------
    source_path:
        Igor data-folder path inside the pxp, e.g. ``"root:SAXS:Sample01"``.
    output_path:
        Filesystem path of the NeXus file written, or ``None`` if skipped.
    technique:
        Canonical technique string (``"USAXS"`` / ``"SAXS"`` / ``"WAXS"``)
        or ``"UNKNOWN"`` for skipped folders.
    n_points:
        Number of Q points in the exported dataset (0 if skipped).
    status:
        ``"ok"`` | ``"skipped"`` | ``"error"``.
    message:
        Human-readable explanation (reason for skip / error text / "").
    """
    source_path: str
    output_path: Path | None
    technique: str
    n_points: int
    status: str
    message: str = ""


@dataclass
class ExtractionResult:
    """Overall result of an :func:`extract_pxp_to_nexus` call."""
    pxp_path: Path
    output_root: Path
    files: list[ExtractedFile] = field(default_factory=list)
    n_unparseable_records: int = 0

    @property
    def n_written(self) -> int:
        return sum(1 for f in self.files if f.status == "ok")

    @property
    def n_skipped(self) -> int:
        return sum(1 for f in self.files if f.status == "skipped")

    @property
    def n_errors(self) -> int:
        return sum(1 for f in self.files if f.status == "error")


# ---------------------------------------------------------------------------
# Defensive .pxp loader
# ---------------------------------------------------------------------------

def _load_pxp_filesystem(pxp_path: Path) -> tuple[dict, int]:
    """Load a .pxp file and return ``(filesystem, n_skipped_records)``.

    Unlike :func:`igor2.packed.load`, this never aborts on a single bad
    wave record: the offending record is skipped (counted) and traversal
    continues, so a legacy experiment with one weird wave still yields a
    fully usable folder tree. Folder structure depends only on
    ``FolderStartRecord`` / ``FolderEndRecord`` records, which parse
    cheaply and reliably.
    """
    try:
        from igor2.packed import (
            setup_packed_file_record_header, _RECORD_TYPE,
            _UnknownRecord, _UnusedRecord, PACKEDRECTYPE_MASK,
            _byte_order, _need_to_reorder_bytes,
        )
        from igor2.record.folder import FolderStartRecord, FolderEndRecord
        from igor2.record.wave import WaveRecord
    except ImportError as e:
        raise ImportError(
            "The 'igor2' package is required to import Igor packed experiments. "
            "Install it with:  pip install igor2"
        ) from e

    records: list[Any] = []
    byte_order: str | None = None
    header_struct = setup_packed_file_record_header()
    skipped = 0

    with open(pxp_path, "rb") as f:
        while True:
            b = f.read(header_struct.size)
            if len(b) == 0:
                break
            if len(b) < header_struct.size:
                break
            header = header_struct.unpack_from(b)
            if header["version"] and not byte_order:
                need = _need_to_reorder_bytes(header["version"])
                byte_order = _byte_order(need)
                if need:
                    header_struct = setup_packed_file_record_header(byte_order=byte_order)
                    header = header_struct.unpack_from(b)
            data = bytes(f.read(header["numDataBytes"]))
            if len(data) < header["numDataBytes"]:
                logger.warning("truncated record in %s — stopping", pxp_path)
                break
            rtype_id = header["recordType"] & PACKEDRECTYPE_MASK
            rtype = _RECORD_TYPE.get(rtype_id, _UnknownRecord)
            try:
                rec = rtype(header, data, byte_order=byte_order)
                records.append(rec)
            except Exception as exc:
                skipped += 1
                records.append(None)
                logger.debug("skipped %s record (rtype_id=%s): %s",
                             rtype.__name__, rtype_id, exc)

    # Build a folder tree from the surviving records.
    filesystem = {"root": {}}
    stack: list[tuple[str, dict]] = [("root", filesystem["root"])]

    for rec in records:
        if rec is None:
            continue
        if isinstance(rec, FolderStartRecord):
            name = rec.null_terminated_text
            if isinstance(name, bytes):
                name = name.decode("utf-8", errors="replace")
            parent = stack[-1][1]
            parent.setdefault(name, {})
            stack.append((name, parent[name]))
        elif isinstance(rec, FolderEndRecord):
            if len(stack) > 1:
                stack.pop()
        elif isinstance(rec, WaveRecord):
            try:
                bname = rec.wave["wave"]["wave_header"]["bname"]
                if isinstance(bname, bytes):
                    bname = bname.decode("utf-8", errors="replace").rstrip("\x00")
                stack[-1][1][bname] = rec
            except Exception:
                # Unreadable wave name — drop silently; structure is the priority.
                pass

    return filesystem, skipped


def _wave_array(wave_rec) -> np.ndarray | None:
    """Extract the numpy data array from an igor2 WaveRecord."""
    try:
        data = wave_rec.wave["wave"]["wData"]
        arr = np.asarray(data)
        # Flatten 1-D wraps; reject text-typed or multi-dim arrays.
        if arr.dtype.kind in ("U", "S", "O"):
            return None
        return np.asarray(arr.ravel(), dtype=np.float64)
    except Exception:
        return None


def _wave_note(wave_rec) -> str:
    """Decode a wave note to a string (empty if missing)."""
    try:
        note = wave_rec.wave["wave"].get("note", b"")
        if isinstance(note, bytes):
            return note.decode("utf-8", errors="replace")
        return str(note)
    except Exception:
        return ""


# ---------------------------------------------------------------------------
# Wave-note → NeXus metadata parser
# ---------------------------------------------------------------------------

# Section markers used by the APS USAXS pipeline. Each entry maps
# (start_marker, end_marker) → metadata key prefix that
# `create_nxcansas_file` understands.
_NOTE_SECTIONS: list[tuple[str, str, str]] = [
    ("NXSampleStart",    "NXSampleEnd",    "sample."),
    ("NXInstrumentStart", "NXInstrumentEnd", "instrument."),
    ("NXUserStart",      "NXUserEnd",      "user."),
    ("NXMetadataStart",  "NXMetadataEnd",  ""),   # general — goes to entry/notes
]


def _parse_kv_pairs(text: str) -> dict[str, str]:
    """Parse Igor-style ``key=value;key=value;`` into a dict.

    Values are stripped; trailing empty entries are dropped. Values that
    look numeric are left as strings — the caller decides what to coerce.
    """
    out: dict[str, str] = {}
    for chunk in text.split(";"):
        chunk = chunk.strip()
        if not chunk or "=" not in chunk:
            continue
        k, _, v = chunk.partition("=")
        k = k.strip()
        v = v.strip()
        if k:
            out[k] = v
    return out


def _coerce(val: str) -> Any:
    """Best-effort string→number coercion. Leaves unparseable strings alone."""
    if val == "" or val.lower() in ("nan", "inf", "-inf"):
        return val
    try:
        if "." in val or "e" in val.lower():
            return float(val)
        return int(val)
    except ValueError:
        return val


def parse_wave_note(note: str) -> dict[str, Any]:
    """Convert an Igor wave note into a NXcanSAS-friendly metadata dict.

    Returns a dict whose keys use the conventions understood by
    :func:`pyirena.io.nxcansas_unified.create_nxcansas_file` —
    ``sample.<key>``, ``instrument.<key>``, plus top-level ``title`` etc.

    The parser handles three layers:

    1. **Section markers** (``NXSampleStart`` … ``NXSampleEnd`` etc.) —
       contents are flattened into ``sample.<key>`` keys.
    2. **Top-level fields outside sections** (e.g. ``DATAFILE=..``,
       ``DATE=..``, ``COMMENT=..``) — mapped to ``title``, ``start_time``,
       and ``notes/<original_key>``.
    3. **Bare units / long_name** (e.g. ``units=arb;long_name=Intensity;``) —
       these are wave-level annotations, not sample metadata, so they are
       dropped from the output dict.
    """
    if not note or "=" not in note:
        return {}

    meta: dict[str, Any] = {}

    # First, peel off each known section. We work on a mutable copy of the
    # note and remove sections as we go so the remainder can be parsed
    # as flat top-level pairs.
    remainder = note
    for start, end, prefix in _NOTE_SECTIONS:
        i_start = remainder.find(start)
        if i_start < 0:
            continue
        i_end = remainder.find(end, i_start + len(start))
        if i_end < 0:
            continue
        section_body = remainder[i_start + len(start): i_end]
        # Strip leading/trailing semicolons that bracket the section.
        section_body = section_body.strip(";").strip()
        section_kv = _parse_kv_pairs(section_body)
        for k, v in section_kv.items():
            meta[f"{prefix}{k}" if prefix else k] = _coerce(v)
        # Remove the entire ``StartMarker..EndMarker`` chunk including
        # the markers themselves.
        remainder = (remainder[:i_start].rstrip(";") + ";"
                     + remainder[i_end + len(end):].lstrip(";"))

    # Remove section sentinels that have no end-marker (defensive).
    for start, end, _ in _NOTE_SECTIONS:
        remainder = remainder.replace(start, "").replace(end, "")
    # Remove the marker pair labels themselves (e.g. ``Nexus_attributesStartHere``)
    for marker in ("Nexus_attributesStartHere", "Nexus_attributesEndHere"):
        remainder = remainder.replace(marker, "")

    # Now parse what's left as flat key=value pairs. Map a few well-known
    # keys to NXcanSAS-canonical names; route everything else into
    # ``notes/<key>`` to preserve provenance without polluting sample/.
    flat = _parse_kv_pairs(remainder)

    # Wave-level annotation keys that should not become entry-level metadata.
    _IGNORE = {"units", "long_name", "Units", "Wname"}

    for k, v in flat.items():
        if k in _IGNORE:
            continue
        kl = k.lower()
        if kl == "datafile":
            meta.setdefault("title", v)
            meta[f"notes.{k}"] = v
        elif kl == "date":
            meta.setdefault("start_time", v)
            meta[f"notes.{k}"] = v
        elif kl == "comment":
            # Use comment as title only if no explicit title already.
            meta.setdefault("title", v)
            meta[f"notes.{k}"] = v
        else:
            meta[f"notes.{k}"] = _coerce(v)

    return meta


# ---------------------------------------------------------------------------
# Filename / folder-name sanitisers
# ---------------------------------------------------------------------------

_INVALID_FS_CHARS = re.compile(r'[<>:"/\\|?*\x00-\x1f]')


def _safe_filename(name: str) -> str:
    """Sanitise an Igor folder name for use as a filesystem path component."""
    cleaned = _INVALID_FS_CHARS.sub("_", name).strip(" .")
    return cleaned or "unnamed"


def _unique_path(path: Path) -> Path:
    """Return *path* if it doesn't exist, else append ``_2``, ``_3``, … to the stem."""
    if not path.exists():
        return path
    n = 2
    while True:
        candidate = path.with_stem(f"{path.stem}_{n}")
        if not candidate.exists():
            return candidate
        n += 1


# ---------------------------------------------------------------------------
# Tree traversal
# ---------------------------------------------------------------------------

def _classify_folder(folder_name: str) -> str | None:
    """Return the canonical technique string for a top-level folder, or None."""
    if folder_name in TECHNIQUE_FOLDERS:
        return TECHNIQUE_FOLDERS[folder_name]
    # Case-insensitive fallback.
    for k, v in TECHNIQUE_FOLDERS.items():
        if folder_name.lower() == k.lower():
            return v
    return None


def _pick_waves(folder: dict, technique: str
                ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray | None, str] | None:
    """Try each picker for *technique*; return (Q, I, err, dq, note) on first hit."""
    pickers = WAVE_PICKERS.get(technique, [])
    for q_name, i_name, e_name, dq_name in pickers:
        if q_name not in folder or i_name not in folder or e_name not in folder:
            continue
        q_rec = folder[q_name]
        i_rec = folder[i_name]
        e_rec = folder[e_name]
        if any(isinstance(r, dict) for r in (q_rec, i_rec, e_rec)):
            continue   # one of them is a sub-folder, not a wave
        q = _wave_array(q_rec)
        intensity = _wave_array(i_rec)
        err = _wave_array(e_rec)
        if q is None or intensity is None or err is None:
            continue
        # Trim to common length defensively (real files match, but be safe)
        n = min(len(q), len(intensity), len(err))
        if n < 2:
            continue
        q, intensity, err = q[:n], intensity[:n], err[:n]

        dq: np.ndarray | None = None
        if dq_name and dq_name in folder and not isinstance(folder[dq_name], dict):
            dq_arr = _wave_array(folder[dq_name])
            if dq_arr is not None and len(dq_arr) >= n:
                dq = dq_arr[:n]

        # Wave note comes from the intensity wave (richest metadata).
        note = _wave_note(i_rec)
        return q, intensity, err, dq, note
    return None


def _walk_tree(node: dict, parent_technique: str | None, rel_path: tuple[str, ...]
               ) -> Iterable[tuple[tuple[str, ...], str, dict]]:
    """Yield ``(rel_path, technique, folder_dict)`` for each sample-bearing folder.

    Walks *node* depth-first. At each folder, attempts wave extraction with
    the inherited technique. Folders that don't yield data are still walked
    into (their children may). Top-level folders that match
    ``TECHNIQUE_FOLDERS`` set/override the inherited technique for the
    subtree.
    """
    for name, child in node.items():
        if not isinstance(child, dict):
            continue
        # Top-level folders set the technique. (Also accept re-declaration
        # mid-tree, e.g. root:Custom:USAXS:..., for robustness.)
        own_tech = _classify_folder(name)
        new_technique = own_tech if own_tech is not None else parent_technique
        new_rel = rel_path + (name,)

        # Try wave extraction here (only meaningful if we know the technique).
        if new_technique is not None and parent_technique is not None:
            # Already inside a technique subtree — try to harvest waves
            yield (new_rel, new_technique, child)
        elif new_technique is not None and own_tech is not None:
            # Entered a technique-root folder; do NOT yield the folder itself
            # (it holds sample sub-folders, not waves). Just descend.
            pass

        # Recurse — children may hold more samples (in-situ sub-runs, etc.)
        yield from _walk_tree(child, new_technique, new_rel)


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def extract_pxp_to_nexus(
    pxp_path: str | Path,
    output_root: str | Path | None = None,
    techniques: Sequence[str] | None = None,
    overwrite: bool = False,
) -> ExtractionResult:
    """Extract reduced 1-D SAS data from a packed Igor experiment into NeXus files.

    Parameters
    ----------
    pxp_path:
        Path to the input ``.pxp`` file.
    output_root:
        Destination folder. If ``None`` (default), a sibling folder named
        ``<pxp_stem>_data`` is created next to the input file. The folder
        is created if it does not exist. Inside, one sub-folder per
        technique (``USAXS/``, ``SAXS/``, ``WAXS/``) is created on demand.
    techniques:
        Optional iterable of technique strings to keep
        (``["USAXS", "SAXS", "WAXS"]``). ``None`` (default) keeps every
        technique present in the file.
    overwrite:
        If ``True``, existing files in the destination are overwritten.
        If ``False`` (default), an unused suffix (``_2``, ``_3``, …) is
        appended to keep both copies.

    Returns
    -------
    ExtractionResult
        Container with per-file summaries (:class:`ExtractedFile`) and
        counts of unparseable records. Use ``.n_written`` /
        ``.n_skipped`` / ``.n_errors`` for quick reporting.

    Notes
    -----
    The extractor is intentionally conservative: it only writes a NeXus
    file when it finds the full Q/I/error triple defined in
    :data:`WAVE_PICKERS` for the relevant technique. Sample folders that
    hold only raw detector counts (typical for "blank" measurements)
    silently appear in the result as ``status="skipped"``.

    The Igor folder hierarchy is **mirrored** in the output:
    ``root:SAXS:in_situ_run_42:step_03:Sa1`` becomes
    ``out/SAXS/in_situ_run_42/step_03/Sa1.h5``. This preserves the
    user's organisation of in-situ experiments.
    """
    pxp_path = Path(pxp_path)
    if not pxp_path.is_file():
        raise FileNotFoundError(f"PXP file not found: {pxp_path}")

    if output_root is None:
        output_root = pxp_path.with_name(f"{pxp_path.stem}_data")
    output_root = Path(output_root)
    output_root.mkdir(parents=True, exist_ok=True)

    keep_techniques: set[str] | None = None
    if techniques is not None:
        keep_techniques = {t.upper() for t in techniques}

    logger.info("Loading %s…", pxp_path)
    filesystem, n_skipped_records = _load_pxp_filesystem(pxp_path)
    logger.info("Loaded %d folders; %d wave records were unparseable",
                _count_folders(filesystem["root"]), n_skipped_records)

    result = ExtractionResult(
        pxp_path=pxp_path,
        output_root=output_root,
        n_unparseable_records=n_skipped_records,
    )

    # Walk and write
    for rel_path, technique, folder in _walk_tree(filesystem["root"], None, ()):
        if keep_techniques and technique.upper() not in keep_techniques:
            continue

        igor_path = "root:" + ":".join(rel_path)

        picked = _pick_waves(folder, technique)
        if picked is None:
            # This folder is in a known technique subtree but doesn't carry
            # the expected wave triple. Skip silently (probably an
            # intermediate / blank folder).
            #
            # Don't emit a skipped row for nested non-leaf folders (would
            # spam the summary). Only record skips for the leaf-looking
            # ones: folders that contain at least one wave.
            has_any_wave = any(not isinstance(v, dict) for v in folder.values())
            if has_any_wave:
                result.files.append(ExtractedFile(
                    source_path=igor_path,
                    output_path=None,
                    technique=technique,
                    n_points=0,
                    status="skipped",
                    message=f"no recognised {technique} wave triple in folder",
                ))
            continue

        q, intensity, err, dq, note = picked
        meta = parse_wave_note(note)

        # Build the output path: <root>/<technique>/<rel_path_without_root_tech>/<sample>.h5
        # rel_path is e.g. ('SAXS', 'Sa10_brRigNP_05mg__0013')  or
        #                 ('SAXS', 'heater_run', 'step03', 'Sa1')
        # The first element is the technique-root folder; strip it.
        sub_components = [_safe_filename(c) for c in rel_path[1:]]
        if not sub_components:
            sub_components = ["unnamed"]
        sample_name = sub_components[-1]
        sample_dir = output_root / technique
        if len(sub_components) > 1:
            sample_dir = sample_dir.joinpath(*sub_components[:-1])
        sample_dir.mkdir(parents=True, exist_ok=True)

        out_path = sample_dir / f"{sample_name}.h5"
        if not overwrite:
            out_path = _unique_path(out_path)

        try:
            create_nxcansas_file(
                filepath=out_path,
                q=q,
                intensity=intensity,
                error=err,
                dq=dq,
                sample_name=sample_name,
                metadata=meta,
            )
        except Exception as exc:
            result.files.append(ExtractedFile(
                source_path=igor_path,
                output_path=out_path,
                technique=technique,
                n_points=len(q),
                status="error",
                message=str(exc),
            ))
            logger.exception("Failed to write %s", out_path)
            continue

        result.files.append(ExtractedFile(
            source_path=igor_path,
            output_path=out_path,
            technique=technique,
            n_points=len(q),
            status="ok",
        ))

    logger.info("Extraction finished: %d written, %d skipped, %d errors",
                result.n_written, result.n_skipped, result.n_errors)
    return result


def _count_folders(node: dict) -> int:
    """Recursively count folders in the filesystem tree."""
    n = 0
    for v in node.values():
        if isinstance(v, dict):
            n += 1 + _count_folders(v)
    return n


# ---------------------------------------------------------------------------
# CLI entry point (for quick testing without the GUI)
# ---------------------------------------------------------------------------

def _main(argv: Sequence[str] | None = None) -> int:
    """``python -m pyirena.io.pxp_to_nexus <file.pxp> [output_dir] [TECH ...]``"""
    import argparse

    parser = argparse.ArgumentParser(
        description="Extract reduced SAS data from an Igor packed experiment "
                    "into per-sample NXcanSAS HDF5 files.",
    )
    parser.add_argument("pxp", help="Path to the .pxp file")
    parser.add_argument("-o", "--output", help="Output folder (default: <pxp>_data)")
    parser.add_argument("-t", "--technique", action="append",
                        help="Only export this technique (USAXS/SAXS/WAXS). "
                             "Repeat for multiple; omit for all.")
    parser.add_argument("--overwrite", action="store_true",
                        help="Overwrite existing output files")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="Verbose logging")

    args = parser.parse_args(argv)

    logging.basicConfig(
        level=logging.INFO if args.verbose else logging.WARNING,
        format="%(levelname)s  %(message)s",
    )

    result = extract_pxp_to_nexus(
        args.pxp,
        output_root=args.output,
        techniques=args.technique,
        overwrite=args.overwrite,
    )

    print(f"\nSource:   {result.pxp_path}")
    print(f"Output:   {result.output_root}")
    print(f"Written:  {result.n_written}")
    print(f"Skipped:  {result.n_skipped}")
    print(f"Errors:   {result.n_errors}")
    if result.n_unparseable_records:
        print(f"Note:     {result.n_unparseable_records} wave record(s) could not be parsed")

    if args.verbose:
        print("\nPer-file summary:")
        for f in result.files:
            tag = {"ok": "OK ", "skipped": "-- ", "error": "ERR"}.get(f.status, "?")
            extra = f" [{f.message}]" if f.message else ""
            out = str(f.output_path) if f.output_path else "(no output)"
            print(f"  {tag} {f.technique:6s} {f.source_path:60s} -> {out}{extra}")

    return 0 if result.n_errors == 0 else 1


if __name__ == "__main__":
    sys.exit(_main())
