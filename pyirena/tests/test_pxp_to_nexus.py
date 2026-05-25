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
    from pyirena.io.pxp_to_nexus import TECHNIQUE_FOLDERS, _classify_folder

    assert _classify_folder("USAXS") == "USAXS"
    assert _classify_folder("SAXS") == "SAXS"
    assert _classify_folder("WAXS") == "WAXS"
    assert _classify_folder("Imported SAXS") == "SAXS"
    assert _classify_folder("unknown_folder") is None

    # Case-insensitive fallback
    assert _classify_folder("usaxs") == "USAXS"
    assert _classify_folder("Saxs") == "SAXS"


def test_wave_pickers_table_lists_dsm_for_usaxs():
    """USAXS picker is DSM_* (per design choice)."""
    from pyirena.io.pxp_to_nexus import WAVE_PICKERS

    usaxs = WAVE_PICKERS["USAXS"]
    assert any(t[:3] == ("DSM_Qvec", "DSM_Int", "DSM_Error") for t in usaxs)

    saxs = WAVE_PICKERS["SAXS"]
    assert any(t[:3] == ("R_Qvec", "R_Int", "R_Error") for t in saxs)


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
    base.write_text("x")
    next_path = _unique_path(base)
    assert next_path.name == "thing_2.h5"


# ---------------------------------------------------------------------------
# Integration tests against a real legacy .pxp
# ---------------------------------------------------------------------------

_REAL_PXP = Path(__file__).resolve().parents[2] / "temp" / "test.pxp"


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
    from pyirena.io.pxp_to_nexus import extract_pxp_to_nexus
    from pyirena.io.hdf5 import find_matching_groups

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
# Helpers
# ---------------------------------------------------------------------------

def _walk_paths(h5root):
    out: list[str] = []
    h5root.visit(out.append)
    return out
