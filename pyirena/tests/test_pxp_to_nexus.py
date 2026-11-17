"""
Unit + integration tests for the pxp → NeXus importer.

The IBW-v5 binary format is intricate (324-byte WaveHeader5 with nested
pointers and platform-dependent padding) and igor2 ships no writer, so
synthesising a fake .pxp in pure Python is fragile and adds no real
coverage. Instead this suite focuses on:

  * Parser unit tests for the wave-note → metadata mapping
    (:func:`parse_wave_note`) — pure-Python, no .pxp file needed.
  * Folder-name classification table behaviour.
  * Integration tests that run against ``temp/test.pxp`` when it exists.
    These are skipped when the file is absent so a clean clone still
    passes CI.

Maintainers: the developer fixture ``temp/test.pxp`` is a real legacy
APS USAXS experiment kept out of git (16 MB binary). When you add new
recognised wave-name patterns or note-section markers, exercise them
locally with this fixture before merging.
"""

from __future__ import annotations

from pathlib import Path

import h5py
import pytest

pytest.importorskip("igor2", reason="igor2 is required for the importer to run")


# ---------------------------------------------------------------------------
# Pure-Python unit tests (no .pxp file required)
# ---------------------------------------------------------------------------

def test_parse_wave_note_full_round_trip():
    """The note parser splits sentinel-bracketed sections correctly and
    routes bare key=value pairs to ``notes.*`` without polluting
    ``sample.*`` / ``instrument.*``.
    """
    from pyirena.io.pxp_to_nexus import parse_wave_note

    note = (
        "DATAFILE=foo.dat;DATE=2026-05-25;COMMENT=hello;"
        "Nexus_attributesStartHere;"
        "NXSampleStart;name=Bar;thickness=4;temperature=22.5;NXSampleEnd;"
        "NXInstrumentStart;name=APS;wavelength=0.6;NXInstrumentEnd;"
        "Nexus_attributesEndHere;"
        "extra_key=extra_value;Wname=DSM_Int;Units=cm2/cm3;"
    )
    meta = parse_wave_note(note)

    # Section content lands in sample.*/instrument.*
    assert meta["sample.name"] == "Bar"
    assert meta["sample.thickness"] == 4
    assert meta["sample.temperature"] == pytest.approx(22.5)
    assert meta["instrument.name"] == "APS"
    assert meta["instrument.wavelength"] == pytest.approx(0.6)

    # Well-known top-level keys mapped to canonical names
    assert meta["title"] == "foo.dat"        # DATAFILE
    assert meta["start_time"] == "2026-05-25"

    # Unknown top-level keys go into notes.*
    assert meta["notes.extra_key"] == "extra_value"

    # Wave-level annotations are dropped
    assert "notes.Wname" not in meta
    assert "notes.Units" not in meta


def test_parse_wave_note_handles_empty_and_plain():
    from pyirena.io.pxp_to_nexus import parse_wave_note

    # Empty / plain notes never crash.
    assert parse_wave_note("") == {}
    assert parse_wave_note("no equals here") == {}

    # A note with only wave-level fields should yield an empty dict.
    bare = "units=arb;long_name=Intensity;"
    assert parse_wave_note(bare) == {}


def test_parse_wave_note_section_without_end_marker():
    """A start marker with no matching end marker should not crash —
    section is ignored and the rest of the note is still parsed."""
    from pyirena.io.pxp_to_nexus import parse_wave_note

    note = "NXSampleStart;name=Lost;DATAFILE=ok.dat;"
    meta = parse_wave_note(note)
    # The malformed section is skipped, but DATAFILE is still picked up.
    assert meta["title"] == "ok.dat"
    assert "sample.name" not in meta


def test_technique_folders_table_is_data_driven():
    """The folder-name classification table is the single source of truth.
    A user can add a new mapping here without touching extractor code."""
    from pyirena.io.pxp_to_nexus import _classify_folder

    assert _classify_folder("USAXS") == "USAXS"
    assert _classify_folder("SAXS") == "SAXS"
    assert _classify_folder("WAXS") == "WAXS"
    assert _classify_folder("Imported SAXS") == "SAXS"
    assert _classify_folder("unknown_folder") is None

    # Case-insensitive fallback
    assert _classify_folder("usaxs") == "USAXS"
    assert _classify_folder("Saxs") == "SAXS"


def _pxp_record(record_type: int, payload: bytes, order: str) -> bytes:
    """One packed-file record: 8-byte header (type, version, size) + payload."""
    import struct
    return struct.pack(order + "Hhl", record_type, 0, len(payload)) + payload


def _long_folder_pxp(tmp_path: Path, order: str) -> tuple[Path, list[str]]:
    """Write a folder-only .pxp the way Igor Pro 8 lays out long names:
    type-9 start for short names, type-101 (uint32 length + bytes) for
    names over 31 chars, each followed by a type-33 folder id, closed by
    type 10. Names are taken from the experiment in AFM-analysis/igor2#28."""
    import struct
    long_names = [
        "AMC_25_crystalline_04_07_26_10min_0036",
        "X6_This is very_long_object_name_which_is_allowed_"
        "AMC_25_crystalline_04_02_26_0030",
    ]
    short_name = "AMC_04_07_26_3min_0033"
    parts = [_pxp_record(9, b"USAXS".ljust(32, b"\x00"), order),
             _pxp_record(33, struct.pack(order + "I", 1), order)]
    for name in long_names:
        raw = name.encode()
        parts.append(_pxp_record(101, struct.pack(order + "I", len(raw)) + raw, order))
        parts.append(_pxp_record(33, struct.pack(order + "I", 2), order))
        parts.append(_pxp_record(10, b"", order))
    parts.append(_pxp_record(9, short_name.encode().ljust(32, b"\x00"), order))
    parts.append(_pxp_record(10, b"", order))
    parts.append(_pxp_record(10, b"", order))
    path = tmp_path / f"long_names_{'le' if order == '<' else 'be'}.pxp"
    path.write_bytes(b"".join(parts))
    return path, long_names + [short_name]


@pytest.mark.parametrize("order", ["<", ">"])
def test_pxp_loader_reads_igor8_long_folder_names(tmp_path, order):
    """igor2 >= 0.5.15 decodes type-101 long folder names; our defensive
    loader must place them in the tree. The type-33 records that follow
    every folder start must not disturb it."""
    from pyirena.io.pxp_to_nexus import _load_pxp_filesystem

    path, expected = _long_folder_pxp(tmp_path, order)
    fs, skipped = _load_pxp_filesystem(path)
    assert set(fs["root"]["USAXS"]) == set(expected)
    assert skipped == 0


@pytest.mark.parametrize("version", ["0.5.15", "0.5.13", "0.4.9"])
def test_pxp_loader_refuses_igor2_older_than_0_5_16(tmp_path, monkeypatch, version):
    """Older igor2 silently drops long-named (v7) waves and misplaces
    long-named folders, so loading must stop with an upgrade message."""
    import igor2

    from pyirena.io.pxp_to_nexus import _load_pxp_filesystem

    monkeypatch.setattr(igor2, "__version__", version)
    path, _ = _long_folder_pxp(tmp_path, "<")
    with pytest.raises(ImportError, match=r"igor2>=0\.5\.16"):
        _load_pxp_filesystem(path)


@pytest.mark.parametrize("version", ["0.5.16", "0.5.17.dev3", "0.6.0", "1.0"])
def test_pxp_loader_accepts_igor2_0_5_16_and_later(tmp_path, monkeypatch, version):
    import igor2

    from pyirena.io.pxp_to_nexus import _load_pxp_filesystem

    monkeypatch.setattr(igor2, "__version__", version)
    path, _ = _long_folder_pxp(tmp_path, "<")
    _load_pxp_filesystem(path)


def test_batch_igor_to_nexus_surfaces_old_igor2(tmp_path, monkeypatch):
    """The batch wrapper turns most failures into None; a too-old igor2
    must reach the caller (GUI) with its upgrade message instead."""
    import igor2

    from pyirena.batch import igor_to_nexus

    monkeypatch.setattr(igor2, "__version__", "0.5.15")
    path, _ = _long_folder_pxp(tmp_path, "<")
    with pytest.raises(ImportError, match="Upgrade"):
        igor_to_nexus(str(path), output_folder=str(tmp_path / "out"))


class _FakeWave:
    """Minimal stand-in for an igor2 WaveRecord: data plus wave note."""

    def __init__(self, data, note: str = ""):
        import numpy as np
        self.wave = {"wave": {"wData": np.asarray(data, dtype=float),
                              "note": note.encode()}}


def _extract_tree(tmp_path, root: dict, pickers: str = "pxp"):
    from pyirena.io.pxp_to_nexus import (
        WAVE_PICKERS,
        WAVE_PICKERS_H5XP,
        _extract_filesystem_to_nexus,
    )
    table = WAVE_PICKERS if pickers == "pxp" else WAVE_PICKERS_H5XP
    return _extract_filesystem_to_nexus(
        source_path=tmp_path / "fake.pxp", filesystem={"root": root},
        output_root=tmp_path / "out", keep_techniques=None, overwrite=True,
        pickers_table=table)


def test_pick_waves_matches_qrs_names_case_insensitively(tmp_path):
    """Nika writes lowercase q_/r_/s_ QRS names; Igor names are
    case-insensitive, so they must match the uppercase picker entry."""
    folder = "Q_2p6Inlet_0008_270_30"
    q = [0.01, 0.02, 0.03]
    root = {"SAXS": {folder: {
        f"q_{folder}": _FakeWave(q),
        f"r_{folder}": _FakeWave([3.0, 2.0, 1.0]),
        f"s_{folder}": _FakeWave([0.3, 0.2, 0.1]),
    }}}
    res = _extract_tree(tmp_path, root)
    assert [(f.technique, f.status) for f in res.files] == [("SAXS", "ok")]


@pytest.mark.parametrize("pickers", ["pxp", "h5xp"])
def test_smr_only_usaxs_imports_as_slit_smeared(tmp_path, pickers):
    """A USAXS folder with only SMR_* waves (desmearing not run) imports,
    and the file carries dQl from the note's SlitLength so readers smear
    the model instead of fitting it as pinhole data."""
    from pyirena.io.hdf5 import readGenericNXcanSAS

    note = "NXMetadataStart;SlitLength=0.026108;NXMetadataEnd;"
    root = {"USAXS": {"GC_SRM3600_0003": {
        "SMR_Qvec": _FakeWave([1e-3, 2e-3, 3e-3]),
        "SMR_Int": _FakeWave([1e4, 1e3, 1e2], note),
        "SMR_Error": _FakeWave([1e2, 1e1, 1.0]),
    }}}
    res = _extract_tree(tmp_path, root, pickers)
    (f,) = res.files
    assert f.status == "ok"
    data = readGenericNXcanSAS(str(f.output_path.parent), f.output_path.name)
    assert data["is_slit_smeared"] is True
    assert data["slit_length"] == pytest.approx(0.026108)


def test_smr_without_slit_length_is_refused(tmp_path):
    """Slit-smeared data with no SlitLength must not be written as if it
    were pinhole data."""
    root = {"USAXS": {"S1": {
        "SMR_Qvec": _FakeWave([1e-3, 2e-3, 3e-3]),
        "SMR_Int": _FakeWave([1e4, 1e3, 1e2]),
        "SMR_Error": _FakeWave([1e2, 1e1, 1.0]),
    }}}
    res = _extract_tree(tmp_path, root)
    (f,) = res.files
    assert f.status == "error"
    assert f.output_path is None
    assert "SlitLength" in f.message


def test_dsm_preferred_over_smr_and_not_marked_smeared(tmp_path):
    from pyirena.io.hdf5 import readGenericNXcanSAS

    note = "NXMetadataStart;SlitLength=0.026108;NXMetadataEnd;"
    waves = {}
    for pre, scale in (("DSM", 1.0), ("SMR", 2.0)):
        waves[f"{pre}_Qvec"] = _FakeWave([1e-3, 2e-3, 3e-3])
        waves[f"{pre}_Int"] = _FakeWave([scale, scale, scale], note)
        waves[f"{pre}_Error"] = _FakeWave([0.1, 0.1, 0.1])
    res = _extract_tree(tmp_path, {"USAXS": {"S1": waves}})
    (f,) = res.files
    data = readGenericNXcanSAS(str(f.output_path.parent), f.output_path.name)
    assert data["Intensity"][0] == pytest.approx(1.0)
    assert not data["is_slit_smeared"]


@pytest.mark.parametrize("pickers", ["pxp", "h5xp"])
def test_dsm_and_smr_folder_writes_matilda_style_smr_twin(tmp_path, pickers):
    """With both DSM_* and SMR_* waves the file carries the desmeared
    default plus an ``<sample>_SMR`` twin, as Matilda writes, so the
    slit-smeared selector (prefer_slit_smeared) can load either."""
    import h5py

    from pyirena.io.hdf5 import readGenericNXcanSAS

    note = "NXMetadataStart;SlitLength=0.026108;NXMetadataEnd;"
    waves = {"SMR_dQ": _FakeWave([1e-4, 1e-4, 1e-4])}
    for pre, scale in (("DSM", 1.0), ("SMR", 2.0)):
        waves[f"{pre}_Qvec"] = _FakeWave([1e-3, 2e-3, 3e-3])
        waves[f"{pre}_Int"] = _FakeWave([scale, scale, scale], note)
        waves[f"{pre}_Error"] = _FakeWave([0.1, 0.1, 0.1])
    res = _extract_tree(tmp_path, {"USAXS": {"S1": waves}}, pickers)
    (f,) = res.files
    with h5py.File(f.output_path) as h:
        assert h["entry"].attrs["default"] == "S1"
        sas = h["entry/S1_SMR/sasdata"]
        assert sas["Q"].attrs["resolutions"] == "dQw,dQl"
        assert sas["dQw"].shape == (3,)
    folder, name = str(f.output_path.parent), f.output_path.name
    desm = readGenericNXcanSAS(folder, name)
    smr = readGenericNXcanSAS(folder, name, prefer_slit_smeared=True)
    assert desm["Intensity"][0] == pytest.approx(1.0) and not desm["is_slit_smeared"]
    assert smr["Intensity"][0] == pytest.approx(2.0) and smr["is_slit_smeared"]
    assert smr["slit_length"] == pytest.approx(0.026108)


def test_dsm_and_smr_without_slit_length_writes_no_twin(tmp_path):
    import h5py

    waves = {}
    for pre in ("DSM", "SMR"):
        waves[f"{pre}_Qvec"] = _FakeWave([1e-3, 2e-3, 3e-3])
        waves[f"{pre}_Int"] = _FakeWave([1.0, 1.0, 1.0])
        waves[f"{pre}_Error"] = _FakeWave([0.1, 0.1, 0.1])
    (f,) = _extract_tree(tmp_path, {"USAXS": {"S1": waves}}).files
    assert f.status == "ok"
    with h5py.File(f.output_path) as h:
        assert "S1_SMR" not in h["entry"]


def test_nika_w_wave_imported_as_q_resolution(tmp_path):
    """Nika stores per-point Q resolution in w_<folder>; it is written
    as the sample's Q resolution (Qdev)."""
    from pyirena.io.hdf5 import readGenericNXcanSAS

    folder = "Sample_270_30"
    root = {"SAXS": {folder: {
        f"q_{folder}": _FakeWave([0.01, 0.02, 0.03]),
        f"r_{folder}": _FakeWave([3.0, 2.0, 1.0]),
        f"s_{folder}": _FakeWave([0.3, 0.2, 0.1]),
        f"w_{folder}": _FakeWave([0.001, 0.002, 0.003]),
    }}}
    (f,) = _extract_tree(tmp_path, root).files
    data = readGenericNXcanSAS(str(f.output_path.parent), f.output_path.name)
    assert data["dQ"] == pytest.approx([0.001, 0.002, 0.003])


def test_wave_pickers_table_lists_dsm_for_usaxs():
    """USAXS picker is DSM_* (per design choice)."""
    from pyirena.io.pxp_to_nexus import WAVE_PICKERS

    usaxs = WAVE_PICKERS["USAXS"]
    assert any(t[:3] == ("DSM_Qvec", "DSM_Int", "DSM_Error") for t in usaxs)

    saxs = WAVE_PICKERS["SAXS"]
    assert any(t[:3] == ("R_Qvec", "R_Int", "R_Error") for t in saxs)


def test_walk_tree_yields_root_level_wave_folders_as_multi_technique():
    """Folders that sit at root level (no USAXS/SAXS/WAXS parent) and
    contain at least one wave should be yielded with all three
    techniques as candidates. This is the 'per-sample' layout used by
    in-situ time-series experiments where one folder per timepoint
    holds USAXS, SAXS, and WAXS waves side by side.
    """
    from pyirena.io.pxp_to_nexus import _walk_tree

    # Simulate a filesystem dict: one sample-by-sample folder at root,
    # one technique-organised structure alongside it.
    fake_wave = object()  # any non-dict value counts as a "wave"
    fs = {
        "AMC_3min_0033": {
            "DSM_Int": fake_wave,
            "R_Int": fake_wave,
        },
        "USAXS": {
            "Sample01": {
                "DSM_Int": fake_wave,
            },
        },
        "Packages": {"junk": {}},   # should be skipped as infrastructure
    }
    yielded = list(_walk_tree(fs, None, ()))
    paths = {rp: techs for rp, techs, _ in yielded}

    assert ("AMC_3min_0033",) in paths
    assert paths[("AMC_3min_0033",)] == ["USAXS", "SAXS", "WAXS"]
    assert ("USAXS", "Sample01") in paths
    assert paths[("USAXS", "Sample01")] == ["USAXS"]
    # Packages/* must not appear at all
    assert not any("Packages" in rp for rp in paths)


def test_infrastructure_folders_skipped():
    """Igor's own Packages/, SavedSampleSets/ etc. must be ignored at
    root level so they don't pollute the summary with skip rows."""
    from pyirena.io.pxp_to_nexus import INFRASTRUCTURE_FOLDERS, _walk_tree

    assert "Packages" in INFRASTRUCTURE_FOLDERS
    assert "SavedSampleSets" in INFRASTRUCTURE_FOLDERS

    fake_wave = object()
    fs = {"Packages": {"AMC_sample": {"R_Int": fake_wave}}}
    yielded = list(_walk_tree(fs, None, ()))
    assert yielded == []


def test_safe_filename_strips_illegal_chars():
    from pyirena.io.pxp_to_nexus import _safe_filename

    assert _safe_filename("normal_name") == "normal_name"
    assert "/" not in _safe_filename("a/b/c")
    assert ":" not in _safe_filename("a:b")
    assert _safe_filename("") == "unnamed"
    assert _safe_filename("   ...") == "unnamed"


def test_unique_path_appends_suffix(tmp_path):
    from pyirena.io.pxp_to_nexus import _unique_path

    base = tmp_path / "thing.h5"
    base.write_text("x", encoding="utf-8")
    next_path = _unique_path(base)
    assert next_path.name == "thing_2.h5"


# ---------------------------------------------------------------------------
# Integration tests against a real legacy .pxp
# ---------------------------------------------------------------------------

_REAL_PXP = Path(__file__).resolve().parents[2] / "temp" / "test.pxp"


_LONG_NAMES_PXP = _REAL_PXP.with_name("TestLongNames_liberalAll.pxp")


@pytest.mark.skipif(not _LONG_NAMES_PXP.is_file(),
                    reason="temp/TestLongNames_liberalAll.pxp not present "
                           "(developer-only fixture)")
def test_real_igor8_experiment_with_long_folder_and_wave_names(tmp_path):
    """Igor Pro 10 experiment with long and liberal (quoted, with spaces)
    folder and wave names, saved by Jan for this fix. igor2 0.5.16 must
    give every long-named wave its own full name, and all seven samples
    must import with the expected point counts."""
    from pyirena.io.pxp_to_nexus import _load_pxp_filesystem, extract_pxp_to_nexus

    fs, skipped = _load_pxp_filesystem(_LONG_NAMES_PXP)
    assert skipped == 0
    long_saxs = "ReactedQ_1p73Outlet_this_is_a_very_long_folder_ 0007_270_30"
    folder = fs["root"]["SAXS"][long_saxs]
    for prefix in ("q_", "r_", "s_", "w_"):
        assert prefix + long_saxs in folder
    assert not any("too long" in k for k in folder)

    res = extract_pxp_to_nexus(_LONG_NAMES_PXP, tmp_path)
    written = {f.source_path: f.n_points for f in res.files if f.status == "ok"}
    assert written == {
        "root:SAXS:ReactedQ_1p73Inlet_0006_270_30": 119,
        f"root:SAXS:{long_saxs}": 119,
        "root:SAXS:ReactedQ_2p6Inlet_0008_270_30": 119,
        "root:USAXS:GC_SRM3600_0003": 500,
        "root:USAXS:ReactedQ_1p73Outlet_this_is a very long folder name 0007": 500,
        "root:USAXS:ReactedQ_2p6Inlet_0008": 500,
        "root:USAXS:ReactedQ_2p6Outlet_0009": 500,
    }
    assert res.n_unparseable_records == 0


@pytest.mark.skipif(not _REAL_PXP.is_file(),
                    reason="temp/test.pxp not present (developer-only fixture)")
def test_real_legacy_pxp_writes_expected_count(tmp_path):
    """The bundled APS USAXS experiment should yield ~50+ NeXus files across
    USAXS, SAXS, and WAXS without any extraction errors."""
    from pyirena.io.pxp_to_nexus import extract_pxp_to_nexus

    result = extract_pxp_to_nexus(_REAL_PXP, output_root=tmp_path / "out")

    assert result.n_errors == 0
    assert result.n_written >= 50

    # All three techniques produced something.
    techs = {f.technique for f in result.files if f.status == "ok"}
    assert techs == {"USAXS", "SAXS", "WAXS"}


@pytest.mark.skipif(not _REAL_PXP.is_file(),
                    reason="temp/test.pxp not present (developer-only fixture)")
def test_real_legacy_pxp_usaxs_has_qdev_and_metadata(tmp_path):
    """A USAXS file written from the bundled experiment must carry the dQ
    resolution column AND the rich NeXus sample/instrument metadata that
    the APS USAXS pipeline embeds in wave notes."""
    from pyirena.io.pxp_to_nexus import extract_pxp_to_nexus

    result = extract_pxp_to_nexus(_REAL_PXP, output_root=tmp_path / "out",
                                  techniques=["USAXS"])
    assert result.n_written > 0

    # Verify the structural invariants on every USAXS file: Q + I + Idev
    # must be present, and Qdev must be present on at least some files
    # (DSM_dQ is optional per wave; APS USAXS does emit it). At least one
    # file must carry the NeXus-style sample/instrument metadata.
    n_with_qdev = 0
    n_with_sample_name = 0
    n_with_wavelength = 0
    for fr in result.files:
        if fr.status != "ok" or fr.output_path is None:
            continue
        with h5py.File(fr.output_path, "r") as f:
            sas_paths = [p for p in _walk_paths(f)
                         if p.endswith("/sasdata") and isinstance(f[p], h5py.Group)]
            assert sas_paths, f"no sasdata group in {fr.output_path}"
            sas = f[sas_paths[0]]
            assert "Q" in sas and "I" in sas and "Idev" in sas
            if "Qdev" in sas:
                n_with_qdev += 1
            if "entry/sample/name" in f:
                n_with_sample_name += 1
            if "entry/instrument/beam/incident_wavelength" in f:
                n_with_wavelength += 1

    assert n_with_qdev > 0, "no USAXS files had a Qdev column"
    assert n_with_sample_name > 0, "no USAXS files had sample.name metadata"
    assert n_with_wavelength > 0, "no USAXS files had wavelength metadata"


@pytest.mark.skipif(not _REAL_PXP.is_file(),
                    reason="temp/test.pxp not present (developer-only fixture)")
def test_real_legacy_pxp_round_trips_through_pyirena_reader(tmp_path):
    """Every NeXus file the importer writes must be discoverable by
    pyirena's own NXcanSAS group locator. If this breaks, the imported
    files won't appear in the Data Selector even though they're valid
    HDF5."""
    from pyirena.io.hdf5 import find_matching_groups
    from pyirena.io.pxp_to_nexus import extract_pxp_to_nexus

    result = extract_pxp_to_nexus(_REAL_PXP, output_root=tmp_path / "out")
    n_checked = 0
    for fr in result.files:
        if fr.status != "ok" or fr.output_path is None:
            continue
        with h5py.File(fr.output_path, "r") as f:
            paths = find_matching_groups(f, {"canSAS_class": "SASdata"}, {})
            assert paths, f"no SASdata group in {fr.output_path}"
            sas = f[paths[0]]
            assert sas["I"].shape == (fr.n_points,)
        n_checked += 1
    assert n_checked >= 1


@pytest.mark.skipif(not _REAL_PXP.is_file(),
                    reason="temp/test.pxp not present (developer-only fixture)")
def test_real_legacy_pxp_technique_filter(tmp_path):
    """`techniques=['SAXS']` must skip USAXS and WAXS entirely."""
    from pyirena.io.pxp_to_nexus import extract_pxp_to_nexus

    result = extract_pxp_to_nexus(_REAL_PXP, output_root=tmp_path / "out",
                                  techniques=["SAXS"])
    techs = {f.technique for f in result.files if f.status == "ok"}
    assert techs == {"SAXS"}
    assert not (tmp_path / "out" / "USAXS").exists()
    assert not (tmp_path / "out" / "WAXS").exists()


# ---------------------------------------------------------------------------
# h5xp tests — use pyirena's own writer to synthesise fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def synthetic_h5xp(tmp_path: Path) -> Path:
    """Create a small .h5xp containing one USAXS, one SAXS, one WAXS sample.

    Uses :func:`pyirena.io.h5xp_writer.create_h5xp` + ``write_iq_data``,
    which is the same code path Igor users hit via the Data Explorer's
    "Export to Igor" — guarantees the fixture matches what the importer
    will encounter in the wild.
    """
    import numpy as np

    from pyirena.io.h5xp_writer import create_h5xp, write_iq_data

    pxp_path = tmp_path / "synthetic.h5xp"
    n = 32
    q = np.logspace(-3, 0, n)
    intensity = 1.0 / (1.0 + (q * 30) ** 2)
    err = 0.1 * intensity
    dq = 0.05 * q

    with create_h5xp(pxp_path) as f:
        write_iq_data(f, "Sample_U", q, intensity, error=err, dq=dq,
                      wave_note={"SampleName": "Sample_U", "Title": "USAXS sample",
                                 "Wavelength_A": 0.5904},
                      category="USAXS")
        write_iq_data(f, "Sample_S", q, intensity, error=err,
                      wave_note={"SampleName": "Sample_S"},
                      category="SAXS")
        write_iq_data(f, "Sample_W", q, intensity, error=err,
                      wave_note={"SampleName": "Sample_W"},
                      category="WAXS")
    return pxp_path


def test_h5xp_extract_writes_one_file_per_sample(synthetic_h5xp, tmp_path):
    from pyirena.io.pxp_to_nexus import extract_h5xp_to_nexus

    out_dir = tmp_path / "out"
    result = extract_h5xp_to_nexus(synthetic_h5xp, output_root=out_dir)

    assert result.n_written == 3
    assert result.n_errors == 0
    assert result.n_unparseable_records == 0
    written = {f.output_path for f in result.files if f.status == "ok"}
    assert (out_dir / "USAXS" / "Sample_U.h5") in written
    assert (out_dir / "SAXS"  / "Sample_S.h5") in written
    assert (out_dir / "WAXS"  / "Sample_W.h5") in written


def test_h5xp_extracted_file_has_q_i_idev_qdev(synthetic_h5xp, tmp_path):
    """USAXS sample written with `dq=` should produce a Qdev column."""
    from pyirena.io.pxp_to_nexus import extract_h5xp_to_nexus

    out_dir = tmp_path / "out"
    extract_h5xp_to_nexus(synthetic_h5xp, output_root=out_dir)

    with h5py.File(out_dir / "USAXS" / "Sample_U.h5", "r") as f:
        sas = f["entry/Sample_U/sasdata"]
        assert "Q" in sas and "I" in sas and "Idev" in sas
        assert "Qdev" in sas

    # SAXS sample had no dq → no Qdev
    with h5py.File(out_dir / "SAXS" / "Sample_S.h5", "r") as f:
        sas = f["entry/Sample_S/sasdata"]
        assert "Qdev" not in sas


def test_h5xp_wave_note_parsed_into_metadata(synthetic_h5xp, tmp_path):
    """SampleName from the h5xp wave note must reach entry/notes/SampleName
    (h5xp writer uses the literal ``SampleName`` key, not the prefixed
    ``sample.name`` form)."""
    from pyirena.io.pxp_to_nexus import extract_h5xp_to_nexus

    out_dir = tmp_path / "out"
    extract_h5xp_to_nexus(synthetic_h5xp, output_root=out_dir)

    with h5py.File(out_dir / "USAXS" / "Sample_U.h5", "r") as f:
        # h5xp writer puts SampleName as a top-level note key (not in a
        # section block), so it lands under entry/notes/SampleName.
        assert "entry/notes" in f
        # The wavelength key from the wave-note dict
        notes = f["entry/notes"]
        note_keys = list(notes.keys())
        # SampleName should be present somewhere in notes
        assert any("SampleName" in k for k in note_keys)


def test_h5xp_round_trips_through_pyirena_reader(synthetic_h5xp, tmp_path):
    """Every NeXus file the importer writes must be discoverable by
    pyirena's own NXcanSAS locator."""
    from pyirena.io.hdf5 import find_matching_groups
    from pyirena.io.pxp_to_nexus import extract_h5xp_to_nexus

    out_dir = tmp_path / "out"
    result = extract_h5xp_to_nexus(synthetic_h5xp, output_root=out_dir)
    n_checked = 0
    for fr in result.files:
        if fr.status != "ok" or fr.output_path is None:
            continue
        with h5py.File(fr.output_path, "r") as f:
            paths = find_matching_groups(f, {"canSAS_class": "SASdata"}, {})
            assert paths, f"no SASdata group in {fr.output_path}"
            sas = f[paths[0]]
            assert sas["I"].shape == (32,)
        n_checked += 1
    assert n_checked == 3


def test_h5xp_technique_filter(synthetic_h5xp, tmp_path):
    from pyirena.io.pxp_to_nexus import extract_h5xp_to_nexus

    out_dir = tmp_path / "out"
    result = extract_h5xp_to_nexus(synthetic_h5xp, output_root=out_dir,
                                   techniques=["USAXS"])
    assert result.n_written == 1
    assert (out_dir / "USAXS").exists()
    assert not (out_dir / "SAXS").exists()
    assert not (out_dir / "WAXS").exists()


def test_h5xp_overwrite_false_appends_suffix(synthetic_h5xp, tmp_path):
    from pyirena.io.pxp_to_nexus import extract_h5xp_to_nexus

    out_dir = tmp_path / "out"
    extract_h5xp_to_nexus(synthetic_h5xp, output_root=out_dir)
    extract_h5xp_to_nexus(synthetic_h5xp, output_root=out_dir, overwrite=False)

    usaxs_files = sorted((out_dir / "USAXS").glob("*.h5"))
    assert len(usaxs_files) == 2
    assert usaxs_files[1].name == "Sample_U_2.h5"


def test_extract_igor_experiment_dispatches_on_extension(synthetic_h5xp, tmp_path):
    """The dispatcher must route .h5xp to the h5xp reader without the
    caller specifying which one to use."""
    from pyirena.io.pxp_to_nexus import extract_igor_experiment

    out_dir = tmp_path / "out"
    result = extract_igor_experiment(synthetic_h5xp, output_root=out_dir)
    assert result.n_written == 3
    assert (out_dir / "USAXS" / "Sample_U.h5").is_file()


def test_extract_igor_experiment_rejects_unknown_extension(tmp_path):
    """A random .csv or .txt file must raise ValueError, not silently
    do nothing."""
    from pyirena.io.pxp_to_nexus import extract_igor_experiment

    junk = tmp_path / "not_an_igor_file.csv"
    junk.write_text("not igor\n", encoding="utf-8")
    with pytest.raises(ValueError, match="unsupported Igor experiment format"):
        extract_igor_experiment(junk)


def test_parse_wave_note_handles_h5xp_colon_format():
    """h5xp notes use ``key:value;`` instead of pxp's ``key=value;``.
    The parser must detect the separator and produce the same output."""
    from pyirena.io.pxp_to_nexus import parse_wave_note

    h5xp_note = (
        "DATAFILE:foo.dat;DATE:2026-05-25;"
        "NXSampleStart;name:Bar;thickness:4;NXSampleEnd;"
        "NXInstrumentStart;wavelength:0.6;NXInstrumentEnd;"
        "IgorFolder:my_folder;Wname:DSM_Int;Units:cm2/cm3;"
    )
    meta = parse_wave_note(h5xp_note)
    assert meta["sample.name"] == "Bar"
    assert meta["sample.thickness"] == 4
    assert meta["instrument.wavelength"] == pytest.approx(0.6)
    assert meta["title"] == "foo.dat"
    assert meta["start_time"] == "2026-05-25"
    # IgorFolder / Wname / Units are wave-level, should not appear
    assert "notes.IgorFolder" not in meta
    assert "notes.Wname" not in meta
    assert "notes.Units" not in meta


# ---------------------------------------------------------------------------
# h5xp — Irena/Nika "Use QRS Names" convention (uppercase Q_/R_/S_ + suffix)
# ---------------------------------------------------------------------------

@pytest.fixture
def qrs_names_h5xp(tmp_path: Path) -> Path:
    """Create an h5xp mirroring the Irena/Nika "Use QRS Names" layout.

    Real desktop SAXS instruments (e.g. Xenocs) reduced through Irena/Nika
    store each dataset as ``Packed Data/SAS/ImportedData/<folder>/`` holding
    waves ``Q_<folder>`` / ``R_<folder>`` / ``S_<folder>`` (uppercase prefix,
    data-folder name as suffix; R = intensity, S = error). This is written
    with raw h5py — not pyirena's writer — so it faithfully reproduces the
    on-disk structure the importer meets in the wild rather than pyirena's
    own ``q_<folder>`` output.
    """
    import numpy as np

    path = tmp_path / "qrs_names.h5xp"
    n = 40
    q = np.logspace(-2.7, -0.5, n)
    intensity = 1.0 / (1.0 + (q * 25) ** 2)
    err = 0.1 * intensity

    samples = ["20250512_18740_0_00003_sub", "20250512_18740_0_00025_sub"]
    with h5py.File(path, "w") as f:
        imported = f.create_group("Packed Data/SAS/ImportedData")
        for s in samples:
            g = imported.create_group(s)
            g.create_dataset(f"Q_{s}", data=q)
            g.create_dataset(f"R_{s}", data=intensity)
            g.create_dataset(f"S_{s}", data=err)
    return path


def test_wave_pickers_h5xp_include_uppercase_qrs_convention():
    """The QRS-names triple must be registered for every technique in the
    h5xp picker table (regression: previously only lowercase ``q_<folder>``
    was present, so Irena-reduced desktop SAXS files imported as 0 samples)."""
    from pyirena.io.pxp_to_nexus import _QRS_PICKER, WAVE_PICKERS, WAVE_PICKERS_H5XP

    qrs = _QRS_PICKER
    assert qrs[:3] == ("Q_<folder>", "R_<folder>", "S_<folder>")
    for tech in ("USAXS", "SAXS", "WAXS"):
        assert qrs in WAVE_PICKERS_H5XP[tech], f"missing from h5xp {tech}"
        assert qrs in WAVE_PICKERS[tech], f"missing from pxp {tech}"


def test_h5xp_qrs_names_convention_imports_all_samples(qrs_names_h5xp, tmp_path):
    """Regression: an Irena "Use QRS Names" h5xp (Q_/R_/S_ + folder suffix)
    must import every sample, not skip them all."""
    from pyirena.io.pxp_to_nexus import extract_h5xp_to_nexus

    out_dir = tmp_path / "out"
    result = extract_h5xp_to_nexus(qrs_names_h5xp, output_root=out_dir)

    assert result.n_written == 2
    assert result.n_skipped == 0
    assert result.n_errors == 0
    # SAS folder maps to the SAXS technique
    written = sorted(p.name for p in (out_dir / "SAXS").rglob("*.h5"))
    assert written == [
        "20250512_18740_0_00003_sub.h5",
        "20250512_18740_0_00025_sub.h5",
    ]


def test_h5xp_qrs_names_data_is_intact(qrs_names_h5xp, tmp_path):
    """The R_ wave must land as intensity and S_ as error (Idev)."""
    from pyirena.io.hdf5 import find_matching_groups
    from pyirena.io.pxp_to_nexus import extract_h5xp_to_nexus

    out_dir = tmp_path / "out"
    extract_h5xp_to_nexus(qrs_names_h5xp, output_root=out_dir)

    f0 = next((out_dir / "SAXS").rglob("*.h5"))
    with h5py.File(f0, "r") as f:
        paths = find_matching_groups(f, {"canSAS_class": "SASdata"}, {})
        assert paths
        sas = f[paths[0]]
        assert sas["Q"].shape == (40,)
        assert sas["I"].shape == (40,)
        assert "Idev" in sas               # S_ wave → error column
        assert "Qdev" not in sas           # no dQ_ wave in the fixture


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _walk_paths(h5root):
    out: list[str] = []
    h5root.visit(out.append)
    return out
