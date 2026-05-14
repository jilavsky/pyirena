"""
pyirena/io/h5xp_writer.py — Write Igor Pro HDF5 packed experiment (.h5xp) files.

h5xp is Wavemetrics' documented HDF5 experiment format (supported since Igor
Pro 9.0).  A third-party writer can create files that Igor opens natively as
a full experiment, including waves (arrays), wave notes (key:value metadata),
and the folder hierarchy that maps to Igor's data-folder tree.

References
----------
Wavemetrics HDF5 packed experiment specification:
  https://docs.wavemetrics.com/content/wavemetrics-igor-pro-hdf5-packed-experiments

File organisation produced by this module
-----------------------------------------
::

    / (root)
    ├── Packed Data/
    │   ├── SAXS/
    │   │   └── <sample_name>/
    │   │       ├── Q     — scattering vector [Å⁻¹]
    │   │       ├── R     — reduced intensity I(Q)
    │   │       ├── S     — statistical uncertainty σ(Q)  (if available)
    │   │       └── dQ    — Q resolution σ_Q              (if available)
    │   ├── WAXS/
    │   │   └── <sample_name>/     (same QRS structure)
    │   └── Results/
    │       └── <technique>/
    │           ├── SampleNames    — text wave, one entry per file
    │           └── <param>        — float64 wave, one value per file
    ├── History/
    ├── Miscellaneous/
    ├── Packed Procedure Files/
    ├── Recreation/
    └── Symbolic Paths/

Wave notes
----------
Every wave carries an ``IGORWaveNote`` attribute with metadata encoded as
``key:value;key2:value2;`` pairs (string values are double-quoted).  The
pyirena writer uses this to store:

* For QRS data waves: NXcanSAS metadata (sample name, description, date, …).
* For result waves: fit parameters (Rg, G, chi²…) so Igor macros can read
  them without opening the HDF5 structure.

Usage
-----
::

    import numpy as np
    from pyirena.io.h5xp_writer import create_h5xp, write_iq_data, write_result_wave

    with create_h5xp("output.h5xp") as f:
        write_iq_data(f, "my_sample", q, intensity, error=sigma,
                      wave_note={"SampleName": "my_sample", "Date": "2026-05-14"})

        write_result_wave(f, "my_sample", "UnifiedFitIntensity", q_model, I_model,
                          params={"Rg_L1": 25.3, "G_L1": 1.2e-3, "chi_squared": 0.004})

    # For appending to an existing file:
    with open_h5xp("output.h5xp") as f:
        write_iq_data(f, "second_sample", q2, I2)
"""

from __future__ import annotations

import time
from contextlib import contextmanager
from pathlib import Path
from typing import Any, Sequence

import h5py
import numpy as np

# ---------------------------------------------------------------------------
# Igor Pro constants
# ---------------------------------------------------------------------------

# Igor wave type codes (IGORWaveType attribute).
_WT_TEXT    = np.int32(0)    # text wave (variable-length strings)
_WT_FLOAT32 = np.int32(2)    # single precision float (NT_FP32)
_WT_FLOAT64 = np.int32(4)    # double precision float (NT_FP64)

# Igor timestamps: seconds since 1904-01-01 00:00:00 (local time).
# Unix timestamps: seconds since 1970-01-01 00:00:00 UTC.
# Offset = (1970 - 1904) years in seconds = 24107 days × 86400 s/day.
_IGOR_EPOCH_OFFSET: float = 24107 * 86400  # 2 082 844 800

# Minimum Igor version required to load h5xp files.
_REQUIRED_VERSION = "9.0"
_FILE_VERSION     = "9.05"

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _igor_now() -> np.float64:
    """Current wall-clock time as an Igor Pro local timestamp."""
    return np.float64(time.time() + _IGOR_EPOCH_OFFSET)


def make_wave_note(params: dict[str, Any]) -> str:
    """Encode a dict as an Igor wave-note string ``key:value;key2:value2;``.

    String values are double-quoted if they contain ``;`` or ``:``.
    Numeric values are formatted with up to 8 significant figures.
    """
    parts: list[str] = []
    for k, v in params.items():
        key = str(k).replace(":", "_").replace(";", "_")
        if isinstance(v, str):
            needs_quotes = ";" in v or ":" in v
            val = f'"{v}"' if needs_quotes else v
        elif isinstance(v, float) or isinstance(v, np.floating):
            val = f"{v:.8g}"
        elif isinstance(v, (int, np.integer)):
            val = str(int(v))
        elif isinstance(v, bool):
            val = "1" if v else "0"
        else:
            val = str(v)
        parts.append(f"{key}:{val}")
    return ";".join(parts) + ";" if parts else ""


def _write_wave(
    parent: h5py.Group,
    name: str,
    data: np.ndarray,
    wave_note: str = "",
    wave_type: np.int32 | None = None,
) -> h5py.Dataset:
    """Write one Igor wave into *parent* group.

    Parameters
    ----------
    parent:
        Open HDF5 group to write the dataset into.
    name:
        Wave (dataset) name inside *parent*.
    data:
        numpy array.  float64 arrays are written as-is; others are cast to
        float64 unless *wave_type* is supplied.
    wave_note:
        Igor wave-note string (``key:value;`` pairs).  Written as UTF-8 bytes.
    wave_type:
        Override ``IGORWaveType`` attribute.  If None, inferred from dtype.
    """
    if name in parent:
        del parent[name]

    ds = parent.create_dataset(name, data=data)
    ts = _igor_now()
    ds.attrs["IGORWaveType"]                 = wave_type if wave_type is not None else _WT_FLOAT64
    ds.attrs["IGORWaveCreationDateLocal"]    = ts
    ds.attrs["IGORWaveModificationDateLocal"] = ts
    if wave_note:
        ds.attrs["IGORWaveNote"] = np.bytes_(wave_note)   # fixed-length; Igor rejects vlen strings
    return ds


def _ensure_group(f: h5py.File, path: str) -> h5py.Group:
    """Create HDF5 group at *path* if it does not exist; return it."""
    if path not in f:
        f.create_group(path)
    return f[path]


def _sample_folder(category: str, folder_name: str) -> str:
    """Return the Igor packed-data path for one sample."""
    return f"Packed Data/{category}/{folder_name}"


# ---------------------------------------------------------------------------
# Public API — file creation / opening
# ---------------------------------------------------------------------------

def create_h5xp(path: str | Path, overwrite: bool = False) -> h5py.File:
    """Create a new Igor Pro h5xp file and return an open ``h5py.File``.

    The caller is responsible for closing the file (use as a context manager
    or call ``f.close()`` explicitly).  See also :func:`open_h5xp` for
    appending to an existing file.

    Parameters
    ----------
    path:
        Destination path; must end with ``.h5xp`` (not enforced).
    overwrite:
        If True, silently overwrite an existing file.  If False (default),
        raise ``FileExistsError`` if the file already exists.

    Returns
    -------
    h5py.File
        Open in ``"w"`` mode.
    """
    p = Path(path)
    if p.exists() and not overwrite:
        raise FileExistsError(
            f"{p} already exists. Pass overwrite=True to replace it."
        )
    f = h5py.File(p, "w")
    _write_root_attrs(f)
    _create_skeleton(f)
    return f


@contextmanager
def open_h5xp(path: str | Path):
    """Open an existing h5xp file for appending (context manager).

    Yields an open ``h5py.File`` in ``"a"`` mode; closes on exit.

    Example
    -------
    ::

        with open_h5xp("results.h5xp") as f:
            write_iq_data(f, "sample2", q, I)
    """
    f = h5py.File(Path(path), "a")
    try:
        yield f
    finally:
        f.close()


def _write_root_attrs(f: h5py.File) -> None:
    # String attributes must be fixed-length bytes (not Python str / variable-length).
    # Version numbers are float64, not strings.
    # This mirrors the exact HDF5 types written by Igor Pro 9.
    f.attrs["IGORPlatform"]                 = np.bytes_("Windows")
    f.attrs["IGORArchitecture"]             = np.bytes_("Intel")
    f.attrs["IGORAppBits"]                  = np.int32(64)
    f.attrs["IGORSystemTextEncodingName"]   = np.bytes_("UTF-8")
    f.attrs["IGORSystemTextEncodingCode"]   = np.int32(4)   # 4 = UTF-8 on Windows
    f.attrs["IGORFileVersion"]              = np.bytes_(_FILE_VERSION)
    f.attrs["IGORVersion"]                  = np.float64(float(_FILE_VERSION))
    f.attrs["IGORRequiredVersion"]          = np.float64(float(_REQUIRED_VERSION))
    f.attrs["IGORRequiredVersionReason"]    = np.bytes_(
        "Versions of Igor Pro before 9.00 do not support loading experiments "
        "from HDF5 experiment files."
    )
    f.attrs["IGORBuildNumber"]              = np.int32(0)
    f.attrs["IGORRequiredBuildNumber"]      = np.int32(0)
    f.attrs["IGORDefaultFont"]              = np.bytes_("Arial")


def _create_skeleton(f: h5py.File) -> None:
    """Create the top-level groups and required datasets that a valid h5xp must contain."""
    for grp in ("Packed Data", "History", "Miscellaneous",
                "Packed Procedure Files", "Recreation", "Symbolic Paths"):
        if grp not in f:
            f.create_group(grp)

    # History — scalar fixed-length bytes string (empty = no commands run yet)
    if "History/History" not in f:
        f.create_dataset("History/History",
                         data=np.bytes_(b"\r"))

    # Recreation Procedures — minimal header; Igor uses this to rebuild window layout
    if "Recreation/Recreation Procedures" not in f:
        rec = (
            '// Platform=Windows, IGORVersion=9.050, architecture=Intel, '
            'systemTextEncoding="UTF-8", historyTextEncoding="UTF-8", '
            'procwinTextEncoding="UTF-8", recreationTextEncoding="UTF-8", build=0\r'
            '#pragma TextEncoding = "UTF-8"\r'
            'Silent 101 // use | as bitwise or -- not comment.\r'
        )
        f.create_dataset("Recreation/Recreation Procedures",
                         data=np.bytes_(rec.encode("utf-8")))

    # Packed Procedure Files — minimal Igor procedure window boilerplate
    if "Packed Procedure Files/Procedure" not in f:
        proc = (
            '#pragma TextEncoding = "UTF-8"\r'
            '#pragma rtGlobals=3\t\t// Use modern global access method and strict wave access\r'
            '#pragma DefaultTab={3,20,4}\t// Set default tab width in Igor Pro 9 and later\r'
        )
        ds = f.create_dataset("Packed Procedure Files/Procedure",
                              data=np.bytes_(proc.encode("utf-8")))
        ds.attrs["IGORWindowType"]  = np.bytes_(b"Built-in Procedure Window")
        ds.attrs["IGORWindowTitle"] = np.bytes_(b"Procedure")


# ---------------------------------------------------------------------------
# Public API — writing I(Q) data
# ---------------------------------------------------------------------------

def write_iq_data(
    f: h5py.File,
    folder_name: str,
    q: np.ndarray,
    intensity: np.ndarray,
    error: np.ndarray | None = None,
    dq: np.ndarray | None = None,
    wave_note: dict[str, Any] | str = "",
    category: str = "SAXS",
) -> str:
    """Write I(Q) data as Igor QRS-named waves.

    Writes three (or four) waves into ``Packed Data/<category>/<folder_name>/``:

    ====  ============================================
    Wave  Contents
    ====  ============================================
    Q     Scattering vector (Å⁻¹)
    R     Reduced intensity I(Q)
    S     Statistical uncertainty σ(Q)  (if *error*)
    dQ    Q resolution σ_Q               (if *dq*)
    ====  ============================================

    The same wave note is attached to every wave in the folder.

    Parameters
    ----------
    f:
        Open h5py.File (from :func:`create_h5xp` or :func:`open_h5xp`).
    folder_name:
        Igor data-folder name (typically the HDF5 source file stem).
    q, intensity:
        1-D numpy arrays; both must have the same length.
    error:
        Optional 1-D uncertainty array (same length as *q*).
    dq:
        Optional Q-resolution array (same length as *q*).
    wave_note:
        Metadata as a dict (converted via :func:`make_wave_note`) or a
        pre-formatted ``"key:value;"`` string.
    category:
        ``"SAXS"`` (default) or ``"WAXS"``.  Controls which sub-folder of
        ``Packed Data/`` the waves are placed in.

    Returns
    -------
    str
        HDF5 path of the created sample folder.
    """
    note_str = make_wave_note(wave_note) if isinstance(wave_note, dict) else wave_note
    folder_path = _sample_folder(category, folder_name)
    grp = _ensure_group(f, folder_path)

    q_arr  = np.asarray(q,         dtype=np.float64)
    I_arr  = np.asarray(intensity, dtype=np.float64)

    _write_wave(grp, "Q", q_arr,  note_str)
    _write_wave(grp, "R", I_arr,  note_str)
    if error is not None:
        _write_wave(grp, "S", np.asarray(error, dtype=np.float64), note_str)
    if dq is not None:
        _write_wave(grp, "dQ", np.asarray(dq,    dtype=np.float64), note_str)

    return folder_path


# ---------------------------------------------------------------------------
# Public API — writing model / result arrays
# ---------------------------------------------------------------------------

def write_result_wave(
    f: h5py.File,
    folder_name: str,
    igor_wave_name: str,
    x: np.ndarray,
    y: np.ndarray,
    params: dict[str, Any] | None = None,
    category: str = "SAXS",
) -> None:
    """Write a model or result array as an Igor wave alongside the QRS data.

    The X array is stored as ``<igor_wave_name>_X`` and the Y array as
    ``<igor_wave_name>``.  Fit parameters (Rg, chi² …) go into the wave note
    of the Y wave.

    Parameters
    ----------
    f:
        Open h5py.File.
    folder_name:
        Same sample folder used in :func:`write_iq_data`.
    igor_wave_name:
        Canonical Igor wave name (e.g. ``"UnifiedFitIntensity"``).  Must be a
        valid HDF5 dataset name.
    x, y:
        1-D numpy arrays; both must have the same length.
    params:
        Scalar fit parameters to embed in the wave note.
    category:
        ``"SAXS"`` or ``"WAXS"``.
    """
    note_str = make_wave_note(params or {})
    folder_path = _sample_folder(category, folder_name)
    grp = _ensure_group(f, folder_path)

    _write_wave(grp, igor_wave_name,        np.asarray(y, dtype=np.float64), note_str)
    _write_wave(grp, igor_wave_name + "_X", np.asarray(x, dtype=np.float64), "")


# ---------------------------------------------------------------------------
# Public API — writing the cross-file Results table
# ---------------------------------------------------------------------------

def write_results_table(
    f: h5py.File,
    technique: str,
    param_name: str,
    values: Sequence[float],
    sample_names: Sequence[str],
    units: str = "",
) -> None:
    """Append a parameter-vs-sample array to the Results table.

    Creates or updates ``Packed Data/Results/<technique>/``:

    * ``SampleNames`` — text wave with one entry per sample (string array).
    * ``<param_name>`` — float64 wave with the parameter values.

    The wave note of the parameter wave includes ``units``.

    Parameters
    ----------
    f:
        Open h5py.File.
    technique:
        Igor technique/tool name (e.g. ``"Unified Fit"``, ``"Size Distribution"``).
        Used as the sub-folder name.
    param_name:
        Parameter name in Igor-compatible syntax (e.g. ``"Rg_L1"``, ``"Vf"``).
    values:
        Sequence of floats, one per sample.  ``float("nan")`` for missing data.
    sample_names:
        Sequence of strings matching the *folder_name* values used in
        :func:`write_iq_data`, one per sample.
    units:
        Optional unit string embedded in the wave note.
    """
    # Sanitise technique name for use as an HDF5 group name
    tech_safe = technique.replace("/", "-").replace("\\", "-")
    folder_path = f"Packed Data/Results/{tech_safe}"
    grp = _ensure_group(f, folder_path)

    # --- SampleNames text wave -------------------------------------------
    names_arr = np.array(sample_names, dtype=object)
    if "SampleNames" in grp:
        del grp["SampleNames"]
    ds_names = grp.create_dataset("SampleNames", data=names_arr,
                                   dtype=h5py.string_dtype())
    ts = _igor_now()
    ds_names.attrs["IGORWaveType"]                  = _WT_TEXT
    ds_names.attrs["IGORWaveCreationDateLocal"]     = ts
    ds_names.attrs["IGORWaveModificationDateLocal"] = ts

    # --- Parameter float wave --------------------------------------------
    note_str = make_wave_note({"units": units} if units else {})
    vals_arr = np.array(values, dtype=np.float64)
    _write_wave(grp, param_name, vals_arr, note_str)


# ---------------------------------------------------------------------------
# Public API — convenience: extract NXcanSAS metadata as wave note dict
# ---------------------------------------------------------------------------

def nxcansas_metadata(h5_file: h5py.File, entry_path: str = "entry") -> dict[str, Any]:
    """Read basic NXcanSAS metadata from *h5_file* and return as a dict.

    Suitable for passing directly to *wave_note* in :func:`write_iq_data`.

    Parameters
    ----------
    h5_file:
        Open h5py.File pointing to an NXcanSAS file.
    entry_path:
        HDF5 path of the NXentry group (default ``"entry"``).

    Returns
    -------
    dict
        Keys: SampleName, Description, Title, StartTime, Instrument.
        Missing fields are omitted silently.
    """
    meta: dict[str, Any] = {}
    entry = h5_file.get(entry_path)
    if entry is None:
        return meta

    def _str(ds_path: str) -> str | None:
        ds = entry.get(ds_path)
        if ds is None:
            return None
        v = ds[()]
        if isinstance(v, bytes):
            return v.decode("utf-8", errors="replace")
        if isinstance(v, np.ndarray):
            v = v.flat[0]
            if isinstance(v, bytes):
                return v.decode("utf-8", errors="replace")
        return str(v)

    for key, path in [
        ("SampleName",  "sample/name"),
        ("Description", "sample/description"),
        ("Title",       "title"),
        ("StartTime",   "start_time"),
        ("Instrument",  "instrument/name"),
    ]:
        val = _str(path)
        if val:
            meta[key] = val

    return meta
