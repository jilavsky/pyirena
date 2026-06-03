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
    "extract_h5xp_to_nexus",
    "extract_igor_experiment",
    "ExtractedFile",
    "ExtractionResult",
    "TECHNIQUE_FOLDERS",
    "WAVE_PICKERS",
    "WAVE_PICKERS_H5XP",
    "parse_wave_note",
]


logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Configuration tables (data-driven, easy to extend)
# ---------------------------------------------------------------------------

#: Top-level Igor folder names that are infrastructure / configuration
#: and never carry sample data. Skipped entirely (their subtrees are not
#: walked). Add here if your experiment uses other reserved folders.
INFRASTRUCTURE_FOLDERS: frozenset[str] = frozenset({
    "Packages",          # Igor's own package data — never sample data
    "SavedSampleSets",   # Indra/Irena sample-set bookkeeping
    "WindowMacros",
    "WMTemp",
})


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


#: Wave-name picker per technique for **.pxp** experiments. Each entry is a
#: list of ``(Q_name, I_name, Err_name, dQ_name_or_None)`` tuples tried in
#: order. The first tuple whose Q, I, and Err waves all exist in a sample
#: folder wins. dQ is optional — if its name is None or not present, the
#: resulting file will not contain a Qdev dataset.
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


#: Wave-name picker per technique for **.h5xp** files (Wavemetrics' HDF5
#: packed-experiment format). Same tuple shape as :data:`WAVE_PICKERS`,
#: but the literal token ``"<folder>"`` inside a name is substituted with
#: the actual sample folder name at lookup time. This matches the two
#: conventions produced by :mod:`pyirena.io.h5xp_writer`:
#:
#: * lowercase ``q_<folder>`` / ``r_<folder>`` / ``s_<folder>`` /
#:   ``dq_<folder>`` — the per-sample-folder default of
#:   :func:`pyirena.io.h5xp_writer.write_iq_data`
#: * uppercase ``Q`` / ``R`` / ``S`` / ``dQ`` — older Igor-side export
#:   helpers and some hand-built h5xp files
#:
#: Both are tried in order; the first match wins.
WAVE_PICKERS_H5XP: dict[str, list[tuple[str, str, str, str | None]]] = {
    "USAXS": [
        ("q_<folder>", "r_<folder>", "s_<folder>", "dq_<folder>"),
        ("Q", "R", "S", "dQ"),
    ],
    "SAXS": [
        ("q_<folder>", "r_<folder>", "s_<folder>", "dq_<folder>"),
        ("Q", "R", "S", "dQ"),
    ],
    "WAXS": [
        ("q_<folder>", "r_<folder>", "s_<folder>", "dq_<folder>"),
        ("Q", "R", "S", "dQ"),
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
    """Overall result of an :func:`extract_pxp_to_nexus` call.

    Attributes
    ----------
    pxp_path, output_root, files, n_unparseable_records:
        Self-explanatory.
    n_igor8_longname_markers:
        Count of Igor-8-only record types (26 and 33) seen in the
        ``.pxp``. A non-zero value means the experiment contains folder
        or wave names longer than 31 characters that ``igor2`` cannot
        decode — and therefore the importer probably missed entire
        samples whose folder names are too long. Callers should warn
        the user and suggest re-saving the experiment as ``.h5xp``.
        Always 0 for h5xp imports (which have no such limitation).
    """
    pxp_path: Path
    output_root: Path
    files: list[ExtractedFile] = field(default_factory=list)
    n_unparseable_records: int = 0
    n_igor8_longname_markers: int = 0

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

#: Igor binary-wave version that igor2 0.5.x doesn't recognise. v7 was added
#: in Igor Pro 8 to support wave names / dimension labels longer than 31
#: bytes; binary layout is identical to v5 up to the optional sections, so
#: a v7 wave can be parsed as v5 if we patch the version field and ignore
#: the trailing long-name block.
_IGOR_WAVE_VERSION_V5 = 5
_IGOR_WAVE_VERSION_V7 = 7


def _patch_v7_wave_to_v5(data: bytes, byte_order: str) -> bytes:
    """Rewrite the leading version field of a binary-wave payload from
    7 to 5 so igor2 0.5.x will accept it.

    A v7 wave has the same BinHeader5 + WaveHeader5 + wave-data layout as
    v5; the only additions are:
    - ``optionsSize1`` (4 bytes at offset 60) is reinterpreted as
      ``longWaveNameSize`` (2 bytes) + ``optionsSize1`` (2 bytes)
    - a long-wave-name block is appended after the string-indices section
    The numeric wave data — which is all we need for SAS extraction —
    sits at the same offset as in v5 and reads identically. The trailing
    long-name section is harmless because BinHeader5's recorded wfmSize
    tells igor2 when to stop reading the wave data, and our caller never
    inspects what's after the wData section.
    """
    if len(data) < 2:
        return data
    import struct
    ver_fmt = "<h" if byte_order in ("<", "=") else ">h"
    try:
        version = struct.unpack_from(ver_fmt, data, 0)[0]
    except struct.error:
        return data
    if version != _IGOR_WAVE_VERSION_V7:
        return data
    return struct.pack(ver_fmt, _IGOR_WAVE_VERSION_V5) + data[2:]


#: Packed record-type IDs introduced by Igor Pro 8 to support folder/
#: wave names longer than 31 chars. ``igor2 0.5.x`` doesn't know what
#: they are, so it treats them as ``_UnknownRecord`` and skips them.
#: We just count them so the importer can warn the user that the
#: experiment likely has folders we can't see, and recommend re-saving
#: as ``.h5xp`` (which has no name-length limit).
_IGOR8_LONGNAME_RECORD_TYPES: frozenset[int] = frozenset({26, 33})


def _load_pxp_filesystem(pxp_path: Path) -> tuple[dict, int, int]:
    """Load a .pxp file and return
    ``(filesystem, n_skipped_records, n_igor8_longname_markers)``.

    Unlike :func:`igor2.packed.load`, this never aborts on a single bad
    wave record: the offending record is skipped (counted) and traversal
    continues, so a legacy experiment with one weird wave still yields a
    fully usable folder tree. Folder structure depends only on
    ``FolderStartRecord`` / ``FolderEndRecord`` records, which parse
    cheaply and reliably.

    The third return value counts how many ``Unknown`` records of the
    Igor-8 long-name marker types (26, 33) were seen. A non-zero count
    means the experiment was saved by Igor Pro 8 or later **and**
    contains folder or wave names longer than 31 characters, which
    ``igor2 0.5.x`` cannot decode. The high-level summary uses this to
    tell the user to re-save the experiment as ``.h5xp`` (which has no
    such limit).
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
    skipped = 0
    igor8_longname_markers = 0

    # Always force an explicit byte order on the record header. The
    # default ``setup_packed_file_record_header()`` (no argument) uses
    # native byte order *with native alignment padding*, which on some
    # platforms (notably macOS Python builds) inflates the 8-byte logical
    # header to 16 bytes — and immediately misreads numDataBytes from
    # what it thinks is the trailing pad of a record that doesn't have
    # one. Probe the first two bytes to choose endianness, then build
    # the struct with explicit byte order so size is reliably 8 and
    # there's no alignment padding.
    with open(pxp_path, "rb") as f:
        peek = f.read(2)
    if len(peek) < 2:
        return {"root": {}}, 0, 0
    # Igor packed file records start with `short recordType`. Valid
    # record-type IDs are 1..11 (plus high bits), so the first byte is
    # always small; on little-endian it occupies byte 0, on big-endian
    # byte 1. peek[0] != 0 → little-endian; peek[0] == 0 → big-endian.
    initial_byte_order = "<" if peek[0] != 0 else ">"
    header_struct = setup_packed_file_record_header(byte_order=initial_byte_order)
    byte_order: str = initial_byte_order

    # Hard upper bound on a single record. No legitimate Igor wave or
    # procedure-file record approaches this; anything larger means we
    # drifted out of alignment (e.g. on a record type we don't handle)
    # and the next 4 bytes are misread as a record size. Bail rather
    # than allocate gigabytes and OOM.
    file_size = pxp_path.stat().st_size
    _MAX_RECORD_BYTES = 256 * 1024 * 1024  # 256 MB

    with open(pxp_path, "rb") as f:
        while True:
            pos = f.tell()
            b = f.read(header_struct.size)
            if len(b) == 0:
                break
            if len(b) < header_struct.size:
                break
            header = header_struct.unpack_from(b)
            # Refine byte-order detection if the first record carries a
            # non-zero version field (matches upstream igor2 behaviour).
            if header["version"]:
                need = _need_to_reorder_bytes(header["version"])
                refined = _byte_order(need)
                if refined != byte_order:
                    byte_order = refined
                    header_struct = setup_packed_file_record_header(byte_order=byte_order)
                    header = header_struct.unpack_from(b)

            num_bytes = header["numDataBytes"]
            # Sanity-check the record-size field. A negative or absurdly
            # large value almost always means the loop is out of
            # alignment (we read past the end of a record type we don't
            # handle). Stop cleanly so we don't OOM trying to allocate
            # gigabytes of garbage.
            if (num_bytes < 0
                    or num_bytes > _MAX_RECORD_BYTES
                    or pos + header_struct.size + num_bytes > file_size):
                logger.warning(
                    "%s: implausible record size %d at offset 0x%x "
                    "(file is %d bytes). Stopping load; %d records read so far.",
                    pxp_path.name, num_bytes, pos, file_size, len(records),
                )
                break

            data = bytes(f.read(num_bytes))
            if len(data) < num_bytes:
                logger.warning("truncated record in %s — stopping", pxp_path)
                break
            rtype_id = header["recordType"] & PACKEDRECTYPE_MASK
            if rtype_id in _IGOR8_LONGNAME_RECORD_TYPES:
                igor8_longname_markers += 1
            rtype = _RECORD_TYPE.get(rtype_id, _UnknownRecord)
            try:
                rec = rtype(header, data, byte_order=byte_order)
                records.append(rec)
            except Exception as exc:
                # If this looks like a v7 wave (long-name format from
                # Igor Pro 8/10), patch the version field to 5 and retry.
                # igor2 0.5.x doesn't know about v7 even though the
                # binary layout is identical up to the trailing
                # long-name block.
                if rtype is WaveRecord:
                    patched = _patch_v7_wave_to_v5(data, byte_order)
                    if patched is not data:
                        try:
                            rec = WaveRecord(header, patched, byte_order=byte_order)
                            records.append(rec)
                            continue
                        except Exception as exc2:
                            exc = exc2  # fall through to "skipped" path below

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

    return filesystem, skipped, igor8_longname_markers


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


def _parse_kv_pairs(text: str, sep: str = "=") -> dict[str, str]:
    """Parse Igor-style ``key<sep>value;key<sep>value;`` into a dict.

    The two flavours seen in the wild:

    * **.pxp files** (Igor itself + APS USAXS reduction) → ``key=value;``
    * **.h5xp files** (Wavemetrics-documented HDF5 export format, and
      what :mod:`pyirena.io.h5xp_writer` emits) → ``key:value;``

    Pass ``sep=":"`` for h5xp notes. Quoted values
    (``key:"some;string;with;semicolons"``) are honoured for h5xp,
    matching what :func:`pyirena.io.h5xp_writer.make_wave_note` writes.

    Values are stripped of surrounding whitespace and matching outer
    double-quotes; trailing empty chunks are dropped. Numeric coercion is
    the caller's job (see :func:`_coerce`).
    """
    out: dict[str, str] = {}
    # Split on ";" but respect double-quoted segments so a quoted value
    # containing a semicolon is not chopped in half.
    chunks: list[str] = []
    buf: list[str] = []
    in_quote = False
    for ch in text:
        if ch == '"':
            in_quote = not in_quote
            buf.append(ch)
        elif ch == ";" and not in_quote:
            chunks.append("".join(buf))
            buf = []
        else:
            buf.append(ch)
    if buf:
        chunks.append("".join(buf))

    for chunk in chunks:
        chunk = chunk.strip()
        if not chunk or sep not in chunk:
            continue
        k, _, v = chunk.partition(sep)
        k = k.strip()
        v = v.strip()
        # Strip matching outer double-quotes
        if len(v) >= 2 and v[0] == '"' and v[-1] == '"':
            v = v[1:-1]
        if k:
            out[k] = v
    return out


def _detect_kv_separator(note: str) -> str:
    """Return ``"="`` or ``":"`` for a given note. Defaults to ``"="``.

    Decision rule: count delimiters at the top level (outside quotes).
    Whichever wins is the format. Ties go to ``"="`` because that's the
    pxp default and the historical case.
    """
    n_eq = n_co = 0
    in_quote = False
    for ch in note:
        if ch == '"':
            in_quote = not in_quote
        elif not in_quote:
            if ch == "=":
                n_eq += 1
            elif ch == ":":
                n_co += 1
    return ":" if n_co > n_eq else "="


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
    if not note or ("=" not in note and ":" not in note):
        return {}

    # h5xp notes use "key:value;", pxp notes use "key=value;". Auto-detect
    # so callers can stay agnostic.
    sep = _detect_kv_separator(note)

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
        section_kv = _parse_kv_pairs(section_body, sep=sep)
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

    # Now parse what's left as flat key/value pairs. Map a few well-known
    # keys to NXcanSAS-canonical names; route everything else into
    # ``notes/<key>`` to preserve provenance without polluting sample/.
    flat = _parse_kv_pairs(remainder, sep=sep)

    # Wave-level annotation keys that should not become entry-level metadata.
    _IGNORE = {"units", "long_name", "Units", "Wname", "IgorFolder"}

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


def _pick_waves(folder: dict, technique: str,
                pickers_table: dict | None = None,
                folder_name: str | None = None,
                ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray | None, str] | None:
    """Try each picker for *technique*; return (Q, I, err, dq, note) on first hit.

    *pickers_table* defaults to :data:`WAVE_PICKERS` (pxp conventions).
    Pass :data:`WAVE_PICKERS_H5XP` for h5xp imports.

    The literal token ``"<folder>"`` inside any picker entry is replaced
    with *folder_name* before lookup (used by the h5xp picker to match
    the ``q_<folder>``/``r_<folder>``/``s_<folder>`` naming produced by
    :func:`pyirena.io.h5xp_writer.write_iq_data`).
    """
    if pickers_table is None:
        pickers_table = WAVE_PICKERS
    pickers = pickers_table.get(technique, [])

    def _resolve(name: str | None) -> str | None:
        if name is None or folder_name is None:
            return name
        return name.replace("<folder>", folder_name)

    for q_name_raw, i_name_raw, e_name_raw, dq_name_raw in pickers:
        q_name  = _resolve(q_name_raw)
        i_name  = _resolve(i_name_raw)
        e_name  = _resolve(e_name_raw)
        dq_name = _resolve(dq_name_raw)
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
               ) -> Iterable[tuple[tuple[str, ...], list[str], dict]]:
    """Yield ``(rel_path, candidate_techniques, folder_dict)`` for each
    sample-bearing folder discovered in *node*.

    *candidate_techniques* is a list because some Igor experiments
    organise data by sample rather than by technique — a single folder
    at root level can contain USAXS, SAXS, **and** WAXS wave triples
    side by side (the per-sample layout seen e.g. when in-situ time
    series are saved as one folder per timepoint). In that case the
    caller must try every technique's picker against the folder and
    emit one output file per technique whose triple is present.

    Walking rules:
      * Inside a known technique subtree → yield with
        ``[parent_technique]`` so only that one picker is tried.
      * At root level (no parent technique) for a folder that holds at
        least one wave → yield with ``["USAXS", "SAXS", "WAXS"]``.
      * Technique-root folders (``root:USAXS:`` etc.) are not yielded
        themselves; only their children are.

    Top-level folder names matching :data:`TECHNIQUE_FOLDERS` set the
    technique for the subtree (case-insensitive). Re-declaring a
    technique mid-tree is also accepted, for files that nest e.g.
    ``root:Custom:USAXS:Sample``.
    """
    for name, child in node.items():
        if not isinstance(child, dict):
            continue

        # Skip known infrastructure folders entirely (and don't descend
        # into them) so their sub-folders don't pollute the summary with
        # "no wave triple" skip rows.
        if parent_technique is None and name in INFRASTRUCTURE_FOLDERS:
            continue

        own_tech = _classify_folder(name)
        new_technique = own_tech if own_tech is not None else parent_technique
        new_rel = rel_path + (name,)

        if new_technique is not None and parent_technique is not None:
            # Already inside a technique subtree — exactly one candidate.
            yield (new_rel, [new_technique], child)
        elif own_tech is not None:
            # Entered a technique-root folder; do not yield the folder
            # itself (it holds sample sub-folders). Just descend.
            pass
        elif parent_technique is None:
            # Root-level folder with no technique parent. If it contains
            # any waves, treat it as a sample-by-sample folder that may
            # carry up to three technique-specific wave triples.
            has_any_wave = any(not isinstance(v, dict) for v in child.values())
            if has_any_wave:
                yield (new_rel, ["USAXS", "SAXS", "WAXS"], child)

        # Recurse — children may hold more samples (in-situ sub-runs, etc.)
        yield from _walk_tree(child, new_technique, new_rel)


# ---------------------------------------------------------------------------
# h5xp loader
# ---------------------------------------------------------------------------

class _H5Wave:
    """Wave adapter that quacks like an igor2 ``WaveRecord`` for the
    extractor pipeline.

    Wraps an open ``h5py.Dataset`` so :func:`_wave_array` and
    :func:`_wave_note` can treat both pxp wave records and h5xp datasets
    uniformly. Reads the data eagerly (small 1-D arrays) so we don't have
    to keep file handles open during traversal.
    """

    __slots__ = ("wave",)

    def __init__(self, dataset: "h5py.Dataset"):
        # Mirror the nested dict shape that igor2.record.wave.WaveRecord
        # exposes, so _wave_array / _wave_note work without branching.
        try:
            data = np.asarray(dataset[()])
        except Exception:
            data = None
        note_attr = dataset.attrs.get("IGORWaveNote", b"")
        if isinstance(note_attr, np.ndarray):
            note_attr = note_attr.item() if note_attr.size else b""
        if isinstance(note_attr, bytes):
            note_attr = note_attr.decode("utf-8", errors="replace")
        else:
            note_attr = str(note_attr)
        bname = dataset.name.rsplit("/", 1)[-1]
        self.wave = {
            "wave": {
                "wave_header": {"bname": bname},
                "wData": data,
                "note": note_attr,
            }
        }


def _load_h5xp_filesystem(h5xp_path: Path) -> tuple[dict, int]:
    """Load the ``/Packed Data/`` tree from an h5xp file as the same
    nested-dict structure that :func:`_load_pxp_filesystem` produces.

    Only the ``Packed Data/`` subtree is walked — ``History``,
    ``Recreation``, ``Symbolic Paths`` etc. are ignored because they
    never contain reduced 1-D scattering data. The ``Packed Data/Results/``
    subtree (collected scalar tables written by the h5xp batch extractor)
    is also skipped: it holds per-parameter waves indexed by sample,
    which don't match the per-sample-folder shape this extractor expects.

    Returns ``(filesystem, 0)`` — h5xp uses standard HDF5 so there is no
    "unparseable record" failure mode equivalent to the pxp one.
    """
    import h5py

    filesystem: dict = {"root": {}}

    with h5py.File(h5xp_path, "r") as f:
        if "Packed Data" not in f:
            return filesystem, 0
        packed = f["Packed Data"]
        for tech_name in packed:
            tech_grp = packed[tech_name]
            if not isinstance(tech_grp, h5py.Group):
                continue
            if tech_name == "Results":
                continue   # scalar collection — not per-sample data

            tech_dict: dict = {}
            _h5xp_collect_into(tech_grp, tech_dict)
            if tech_dict:
                filesystem["root"][tech_name] = tech_dict

    return filesystem, 0


def _h5xp_collect_into(group: "h5py.Group", out: dict) -> None:
    """Recurse into *group*, mirroring the h5py group/dataset tree into
    a nested dict of dicts/_H5Wave objects that the extractor pipeline
    can consume.
    """
    import h5py

    for name, item in group.items():
        if isinstance(item, h5py.Group):
            child: dict = {}
            _h5xp_collect_into(item, child)
            out[name] = child
        elif isinstance(item, h5py.Dataset):
            # Skip non-numeric / scalar datasets that aren't waves
            # (e.g. text waves carrying SampleNames). The picker drops
            # anything that isn't a 1-D numeric array anyway.
            out[name] = _H5Wave(item)


# ---------------------------------------------------------------------------
# Shared extraction engine
# ---------------------------------------------------------------------------

def _extract_filesystem_to_nexus(
    source_path: Path,
    filesystem: dict,
    output_root: Path,
    keep_techniques: set[str] | None,
    overwrite: bool,
    pickers_table: dict,
    n_unparseable_records: int = 0,
    n_igor8_longname_markers: int = 0,
) -> ExtractionResult:
    """Walk an already-loaded *filesystem* dict and write per-sample NeXus
    files. Shared by both :func:`extract_pxp_to_nexus` and
    :func:`extract_h5xp_to_nexus` — the only difference between them is
    how the filesystem dict was built and which picker table to use.
    """
    result = ExtractionResult(
        pxp_path=source_path,
        output_root=output_root,
        n_unparseable_records=n_unparseable_records,
        n_igor8_longname_markers=n_igor8_longname_markers,
    )

    for rel_path, candidate_techniques, folder in _walk_tree(
        filesystem["root"], None, ()
    ):
        igor_path = "root:" + ":".join(rel_path)
        sample_folder_name = rel_path[-1] if rel_path else None

        # Determine whether the first path component is a technique
        # root (USAXS / SAXS / WAXS) so we know which path components
        # are "sample subfolders" for output-path mirroring.
        first_is_technique_root = (
            len(rel_path) > 0 and _classify_folder(rel_path[0]) is not None
        )

        any_matched = False
        for technique in candidate_techniques:
            if keep_techniques and technique.upper() not in keep_techniques:
                continue
            picked = _pick_waves(folder, technique,
                                  pickers_table=pickers_table,
                                  folder_name=sample_folder_name)
            if picked is None:
                continue
            any_matched = True

            q, intensity, err, dq, note = picked
            meta = parse_wave_note(note)

            # Output path layout:
            #   <root>/<technique>/<sub_components>/<sample_name>.h5
            #
            # If the folder was found *inside* a technique-root parent
            # (rel_path = ('USAXS', 'Sample01')), drop the technique
            # element and mirror the remaining path. Otherwise (root
            # level, e.g. ('AMC_..._0033',)), mirror the whole path
            # under the matched technique.
            if first_is_technique_root:
                source_components = rel_path[1:]
            else:
                source_components = rel_path
            sub_components = [_safe_filename(c) for c in source_components]
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

        if not any_matched:
            # Folder had waves but none of the candidate techniques' wave
            # triples matched. Only record skips for "leaf-ish" folders
            # (those that contain waves) to keep the summary readable
            # for deeply nested trees that don't carry data.
            has_any_wave = any(not isinstance(v, dict) for v in folder.values())
            if has_any_wave:
                tech_label = (candidate_techniques[0]
                              if len(candidate_techniques) == 1
                              else "/".join(candidate_techniques))
                result.files.append(ExtractedFile(
                    source_path=igor_path,
                    output_path=None,
                    technique=tech_label,
                    n_points=0,
                    status="skipped",
                    message=f"no recognised wave triple in folder "
                            f"(tried {', '.join(candidate_techniques)})",
                ))

    logger.info("Extraction finished: %d written, %d skipped, %d errors",
                result.n_written, result.n_skipped, result.n_errors)
    return result


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
    filesystem, n_skipped_records, n_igor8_markers = _load_pxp_filesystem(pxp_path)
    logger.info(
        "Loaded %d folders; %d wave records were unparseable; "
        "%d Igor-8 long-name marker(s) seen",
        _count_folders(filesystem["root"]), n_skipped_records, n_igor8_markers,
    )
    if n_igor8_markers > 0:
        logger.warning(
            "%s contains %d Igor-8 long-name record(s) (folders/waves with "
            "names > 31 chars). igor2 cannot decode these, so some samples "
            "will be missing from the output. Re-save the experiment as "
            ".h5xp from Igor and re-import to capture them.",
            pxp_path.name, n_igor8_markers,
        )

    return _extract_filesystem_to_nexus(
        source_path=pxp_path,
        filesystem=filesystem,
        output_root=output_root,
        keep_techniques=keep_techniques,
        overwrite=overwrite,
        pickers_table=WAVE_PICKERS,
        n_unparseable_records=n_skipped_records,
        n_igor8_longname_markers=n_igor8_markers,
    )


def extract_h5xp_to_nexus(
    h5xp_path: str | Path,
    output_root: str | Path | None = None,
    techniques: Sequence[str] | None = None,
    overwrite: bool = False,
) -> ExtractionResult:
    """Extract reduced 1-D SAS data from an Igor h5xp packed experiment
    into per-sample NeXus files.

    Same shape and contract as :func:`extract_pxp_to_nexus`, but reads
    the Wavemetrics HDF5 packed-experiment format (``.h5xp``) instead of
    the legacy binary ``.pxp`` format. h5xp is what pyirena's own
    :mod:`pyirena.io.h5xp_writer` produces and what modern Igor 9+
    exports — so this is the natural import path for files that already
    came from Python or were saved in the new format.

    Parameters
    ----------
    h5xp_path:
        Path to the input ``.h5xp`` file.
    output_root, techniques, overwrite:
        See :func:`extract_pxp_to_nexus`.

    Returns
    -------
    ExtractionResult
        Same container type, but ``n_unparseable_records`` will always
        be 0 because h5xp uses standard HDF5 and never fails per-record.

    Notes
    -----
    Only the ``/Packed Data/`` subtree is walked. h5xp metadata
    sections (``History``, ``Recreation``, ``Symbolic Paths``,
    ``Packed Procedure Files``) are ignored — they never contain
    reduced 1-D data. The ``Packed Data/Results/`` group (per-parameter
    collected-value waves) is also skipped because it doesn't follow
    the per-sample-folder shape this extractor needs.

    The wave-name convention is ``Q``/``R``/``S``(+``dQ``) per
    :data:`WAVE_PICKERS_H5XP` — matches what
    :func:`pyirena.io.h5xp_writer.write_iq_data` emits.
    """
    h5xp_path = Path(h5xp_path)
    if not h5xp_path.is_file():
        raise FileNotFoundError(f"h5xp file not found: {h5xp_path}")

    if output_root is None:
        output_root = h5xp_path.with_name(f"{h5xp_path.stem}_data")
    output_root = Path(output_root)
    output_root.mkdir(parents=True, exist_ok=True)

    keep_techniques: set[str] | None = None
    if techniques is not None:
        keep_techniques = {t.upper() for t in techniques}

    logger.info("Loading %s…", h5xp_path)
    filesystem, _ = _load_h5xp_filesystem(h5xp_path)
    logger.info("Loaded %d folders from h5xp", _count_folders(filesystem["root"]))

    return _extract_filesystem_to_nexus(
        source_path=h5xp_path,
        filesystem=filesystem,
        output_root=output_root,
        keep_techniques=keep_techniques,
        overwrite=overwrite,
        pickers_table=WAVE_PICKERS_H5XP,
    )


def extract_igor_experiment(
    path: str | Path,
    output_root: str | Path | None = None,
    techniques: Sequence[str] | None = None,
    overwrite: bool = False,
) -> ExtractionResult:
    """Auto-dispatch to :func:`extract_pxp_to_nexus` or
    :func:`extract_h5xp_to_nexus` based on the file extension.

    Use this from GUIs / CLIs where the caller doesn't want to care
    which Igor format the user picked.

    Recognised extensions (case-insensitive): ``.pxp``, ``.pxt``
    (Igor packed binary) → pxp route; ``.h5xp``, ``.hxp``, ``.h5`` whose
    root contains a ``/Packed Data/`` group → h5xp route.
    """
    p = Path(path)
    ext = p.suffix.lower()
    if ext in (".pxp", ".pxt"):
        return extract_pxp_to_nexus(p, output_root, techniques, overwrite)
    if ext in (".h5xp", ".hxp"):
        return extract_h5xp_to_nexus(p, output_root, techniques, overwrite)
    # .h5 fallback: peek at the structure
    if ext in (".h5", ".hdf5"):
        try:
            import h5py
            with h5py.File(p, "r") as f:
                if "Packed Data" in f:
                    return extract_h5xp_to_nexus(p, output_root, techniques, overwrite)
        except Exception:
            pass
    raise ValueError(
        f"{p.name}: unsupported Igor experiment format. "
        f"Expected .pxp, .pxt, .h5xp, or .hxp."
    )


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
    """``python -m pyirena.io.pxp_to_nexus <file.pxp|.h5xp> [-o DIR] [-t TECH ...]``

    Accepts both ``.pxp`` (binary Igor packed) and ``.h5xp`` (Wavemetrics
    HDF5 packed-experiment) inputs; the format is dispatched automatically
    on the file extension.
    """
    import argparse

    parser = argparse.ArgumentParser(
        description="Extract reduced SAS data from an Igor packed experiment "
                    "(.pxp or .h5xp) into per-sample NXcanSAS HDF5 files.",
    )
    parser.add_argument("input", help="Path to the .pxp or .h5xp file")
    parser.add_argument("-o", "--output", help="Output folder (default: <input>_data)")
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

    result = extract_igor_experiment(
        args.input,
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
    if result.n_igor8_longname_markers:
        print(
            f"\n*** WARNING: this .pxp contains "
            f"{result.n_igor8_longname_markers} Igor-8 long-name record(s) "
            f"(folders or waves with names > 31 chars).\n"
            f"  igor2 cannot decode these, so some samples are MISSING from "
            f"the output above.\n"
            f"  Workaround: in Igor Pro, save the experiment as .h5xp instead "
            f"and re-import."
        )

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
