"""Metadata must survive a merge — from both inputs, not just DS1.

A merged file is seeded by copying DS1, so DS1's metadata came across for
free and DS2's was dropped entirely: merging USAXS + SAXS kept everything
recorded about the USAXS measurement and nothing about the SAXS one.  These
tests pin the fix (GitHub issue #21), including the two things that make it
awkward in practice — beamline files disagree on the capitalisation of
``metadata``, and a SAXS/WAXS ``instrument`` group carries the raw 2-D
detector image, which has no business in a merged 1-D curve file.
"""

import h5py
import numpy as np
import pytest

from pyirena.io._nxcansas_common import BULK_DATASET_BYTES, detect_technique
from pyirena.io.nxcansas_data_merge import save_merged_data
from pyirena.io.nxcansas_unified import create_nxcansas_file


@pytest.fixture
def q_i_err():
    q = np.logspace(-3, -0.5, 120)
    intensity = 1e-4 * q**-4 + 0.01
    return q, intensity, 0.02 * intensity


def _add_common(entry, sample_name="sample"):
    """The sample/instrument groups every beamline file carries."""
    s = entry.create_group("sample")
    s.attrs["NX_class"] = "NXsample"
    s.create_dataset("name", data=sample_name)
    s.create_dataset("thickness", data=1.0)
    inst = entry.create_group("instrument")
    inst.attrs["NX_class"] = "NXinstrument"
    mono = inst.create_group("monochromator")
    mono.create_dataset("wavelength", data=0.5904)


def _make_usaxs(path, q, intensity, error, sample_name="usaxs_curve"):
    """A USAXS reduction: lowercase ``metadata``, a ``flyScan`` group.

    *sample_name* becomes the NXcanSAS subentry name, which is separate from
    the ``entry/sample`` NXsample group added below — create_nxcansas_file()
    names the subentry after the sample, so the two collide if you let them.
    """
    create_nxcansas_file(path, q, intensity, error=error, sample_name=sample_name)
    with h5py.File(path, "a") as f:
        entry = f["entry"]
        md = entry.create_group("metadata")
        md.attrs["NX_class"] = "NXcollection"
        md.create_dataset("DCM_energy", data=21.0)
        md.create_dataset("title", data="usaxs scan 42")
        fly = entry.create_group("flyScan")
        fly.attrs["NX_class"] = "NXdata"
        fly.create_dataset("mca3", data=np.arange(10.0))
        _add_common(entry, sample_name="usaxs_sample")
    return path


def _make_swaxs(path, kind, detector_bytes=0):
    """A SAXS or WAXS file: capital ``Metadata``, tilt key picks the technique.

    *detector_bytes* > 0 adds a raw 2-D image under ``instrument/detector/data``
    of roughly that size, the way a real Pilatus/Eiger file does.
    """
    tilt = "pin_ccd_tilt_x" if kind == "saxs" else "waxs_ccd_tilt_x"
    with h5py.File(path, "w") as f:
        entry = f.create_group("entry")
        entry.attrs["NX_class"] = "NXentry"
        md = entry.create_group("Metadata")          # capital M, as the raw files have
        md.attrs["NX_class"] = "NXcollection"
        md.create_dataset(tilt, data=0.25)
        md.create_dataset("StartTime", data="2026-01-02 03:04:05.000000")
        md.create_dataset("I_scaling", data=2.5)
        _add_common(entry, sample_name=f"{kind}_sample")
        if detector_bytes:
            det = f["entry/instrument"].create_group("detector")
            det.attrs["NX_class"] = "NXdetector"
            n = int(detector_bytes // 8)
            det.create_dataset("data", data=np.zeros(n, dtype=np.float64))
            det.create_dataset("distance", data=350.0)
    return path


def _merge(tmp_path, ds1, ds2, q, intensity, error, ds1_is_nxcansas=True,
           suffix="_merged"):
    return save_merged_data(
        output_folder=tmp_path / "out",
        ds1_path=ds1,
        ds1_is_nxcansas=ds1_is_nxcansas,
        q=q, I=intensity, dI=error, dQ=None,
        merge_result_dict={"scale": 1.0},
        ds2_path=ds2,
        output_stem_suffix=suffix,
    )


# ---------------------------------------------------------------------------
# Technique detection
# ---------------------------------------------------------------------------

class TestDetectTechnique:
    def test_saxs_and_waxs_from_tilt_key(self, tmp_path):
        assert detect_technique(_make_swaxs(tmp_path / "s.hdf", "saxs")) == "saxs"
        assert detect_technique(_make_swaxs(tmp_path / "w.hdf", "waxs")) == "waxs"

    def test_usaxs_from_flyscan(self, tmp_path, q_i_err):
        q, i, e = q_i_err
        assert detect_technique(_make_usaxs(tmp_path / "u.h5", q, i, e)) == "usaxs"

    def test_plain_nxcansas_is_unknown(self, tmp_path, q_i_err):
        """No markers at all — say so rather than guess."""
        q, i, e = q_i_err
        p = tmp_path / "plain.h5"
        create_nxcansas_file(p, q, i, error=e, sample_name="plain")
        assert detect_technique(p) == "unknown"

    def test_filename_does_not_decide(self, tmp_path):
        """A SAXS file named 'waxs' is still SAXS — content wins, always."""
        assert detect_technique(_make_swaxs(tmp_path / "waxs_scan.hdf", "saxs")) == "saxs"

    def test_non_hdf5_and_missing_are_not_errors(self, tmp_path):
        text = tmp_path / "curve.dat"
        text.write_text("0.01 100.0\n0.02 50.0\n", encoding="utf-8")
        assert detect_technique(text) == "unknown"
        assert detect_technique(tmp_path / "nope.h5") == "unknown"


# ---------------------------------------------------------------------------
# What lands in the merged file
# ---------------------------------------------------------------------------

class TestMergedMetadata:
    def test_ds2_groups_arrive_under_the_technique_suffix(self, tmp_path, q_i_err):
        q, i, e = q_i_err
        ds1 = _make_usaxs(tmp_path / "u.h5", q, i, e)
        ds2 = _make_swaxs(tmp_path / "s.hdf", "saxs")
        out = _merge(tmp_path, ds1, ds2, q, i, e)

        with h5py.File(out, "r") as f:
            entry = f["entry"]
            # DS1 keeps its plain names …
            assert "metadata" in entry
            assert float(entry["metadata/DCM_energy"][()]) == pytest.approx(21.0)
            # … and DS2 arrives suffixed, from a group spelled 'Metadata'.
            assert "metadata_saxs" in entry
            assert "instrument_saxs" in entry
            assert "sample_saxs" in entry
            assert float(entry["metadata_saxs/pin_ccd_tilt_x"][()]) == pytest.approx(0.25)
            assert float(
                entry["instrument_saxs/monochromator/wavelength"][()]
            ) == pytest.approx(0.5904)
            assert entry["metadata_saxs"].attrs["pyirena_source_file"] == str(ds2)
            # records the source's own spelling, capital M and all
            assert entry["metadata_saxs"].attrs["pyirena_source_group"] == "/entry/Metadata"

    def test_ds1_metadata_is_not_clobbered(self, tmp_path, q_i_err):
        """Carrying DS2 across must not touch what DS1 brought."""
        q, i, e = q_i_err
        ds1 = _make_usaxs(tmp_path / "u.h5", q, i, e)
        ds2 = _make_swaxs(tmp_path / "s.hdf", "saxs")
        out = _merge(tmp_path, ds1, ds2, q, i, e)
        with h5py.File(out, "r") as f:
            assert f["entry/metadata/title"][()].decode() == "usaxs scan 42"
            assert "pin_ccd_tilt_x" not in f["entry/metadata"]
            assert f["entry/sample/name"][()].decode() == "usaxs_sample"

    def test_unknown_technique_falls_back_to_ds2(self, tmp_path, q_i_err):
        q, i, e = q_i_err
        ds1 = _make_usaxs(tmp_path / "u.h5", q, i, e)
        ds2 = tmp_path / "mystery.h5"
        create_nxcansas_file(ds2, q, i, error=e, sample_name="mystery")
        with h5py.File(ds2, "a") as f:
            f["entry"].create_group("metadata").create_dataset("note", data="?")
        out = _merge(tmp_path, ds1, ds2, q, i, e)
        with h5py.File(out, "r") as f:
            assert "metadata_ds2" in f["entry"]
            assert "metadata_unknown" not in f["entry"]

    def test_raw_detector_image_is_left_behind(self, tmp_path, q_i_err):
        """~2 MB of detector image must not ride along into a 1-D curve file."""
        q, i, e = q_i_err
        ds1 = _make_usaxs(tmp_path / "u.h5", q, i, e)
        ds2 = _make_swaxs(tmp_path / "s.hdf", "saxs",
                          detector_bytes=2 * BULK_DATASET_BYTES)
        out = _merge(tmp_path, ds1, ds2, q, i, e)

        with h5py.File(out, "r") as f:
            inst = f["entry/instrument_saxs"]
            assert "detector" in inst
            # the small sibling survives, the bulk array does not
            assert float(inst["detector/distance"][()]) == pytest.approx(350.0)
            assert "data" not in inst["detector"]
            skipped = [s.decode() if isinstance(s, bytes) else str(s)
                       for s in inst.attrs["pyirena_skipped"]]
            assert any(s.endswith("/detector/data") for s in skipped)
        assert out.stat().st_size < BULK_DATASET_BYTES

    def test_text_ds2_is_silently_skipped(self, tmp_path, q_i_err):
        q, i, e = q_i_err
        ds1 = _make_usaxs(tmp_path / "u.h5", q, i, e)
        ds2 = tmp_path / "curve.dat"
        ds2.write_text("0.01 100.0\n", encoding="utf-8")
        out = _merge(tmp_path, ds1, ds2, q, i, e)   # must not raise
        with h5py.File(out, "r") as f:
            assert "metadata" in f["entry"]
            assert not [k for k in f["entry"] if k.startswith("metadata_")]

    def test_fresh_output_still_gets_both_sides(self, tmp_path, q_i_err):
        """DS1 not NXcanSAS: the output is built from scratch and starts empty."""
        q, i, e = q_i_err
        ds1 = _make_usaxs(tmp_path / "u.h5", q, i, e)
        ds2 = _make_swaxs(tmp_path / "w.hdf", "waxs")
        out = _merge(tmp_path, ds1, ds2, q, i, e, ds1_is_nxcansas=False)
        with h5py.File(out, "r") as f:
            entry = f["entry"]
            assert "metadata" in entry          # DS1, copied explicitly
            assert float(entry["metadata/DCM_energy"][()]) == pytest.approx(21.0)
            assert "metadata_waxs" in entry     # DS2

    def test_resaving_does_not_duplicate(self, tmp_path, q_i_err):
        q, i, e = q_i_err
        ds1 = _make_usaxs(tmp_path / "u.h5", q, i, e)
        ds2 = _make_swaxs(tmp_path / "s.hdf", "saxs")
        _merge(tmp_path, ds1, ds2, q, i, e)
        out = _merge(tmp_path, ds1, ds2, q, i, e)   # same output path again
        with h5py.File(out, "r") as f:
            assert float(
                f["entry/metadata_saxs/pin_ccd_tilt_x"][()]
            ) == pytest.approx(0.25)


class TestChainedMerge:
    """USAXS + SAXS, then that + WAXS — the case the issue is actually about."""

    def test_three_techniques_end_up_side_by_side(self, tmp_path, q_i_err):
        q, i, e = q_i_err
        usaxs = _make_usaxs(tmp_path / "u.h5", q, i, e)
        saxs = _make_swaxs(tmp_path / "s.hdf", "saxs")
        waxs = _make_swaxs(tmp_path / "w.hdf", "waxs")

        first = _merge(tmp_path, usaxs, saxs, q, i, e, suffix="_us")
        second = _merge(tmp_path, first, waxs, q, i, e, suffix="_usw")

        with h5py.File(second, "r") as f:
            entry = f["entry"]
            assert "metadata" in entry           # USAXS, carried twice over
            assert "metadata_saxs" in entry      # survived the second copy
            assert "metadata_waxs" in entry      # added by the second merge
            assert float(entry["metadata/DCM_energy"][()]) == pytest.approx(21.0)
            assert float(entry["metadata_saxs/pin_ccd_tilt_x"][()]) == pytest.approx(0.25)
            assert float(entry["metadata_waxs/waxs_ccd_tilt_x"][()]) == pytest.approx(0.25)


class TestNameCollisions:
    """The propagated name can land on top of the merged data itself.

    ``create_nxcansas_file()`` names the NXcanSAS subentry after the sample,
    so a DS1 called ``sample.h5`` produces a fresh output whose merged curve
    lives at ``entry/sample`` — precisely where DS1's NXsample group would be
    written.  Replacing it would delete the merged curve.
    """

    def test_merged_data_is_never_overwritten(self, tmp_path, q_i_err):
        q, i, e = q_i_err
        # stem "sample" ⇒ the fresh output's subentry is entry/sample,
        # and DS1 also carries an entry/sample NXsample group.
        ds1 = _make_usaxs(tmp_path / "sample.h5", q, i, e)
        ds2 = _make_swaxs(tmp_path / "s.hdf", "saxs")
        out = _merge(tmp_path, ds1, ds2, q, i, e, ds1_is_nxcansas=False)

        with h5py.File(out, "r") as f:
            # the merged curve is still there and still right
            np.testing.assert_allclose(f["entry/sample/sasdata/I"][:], i)
            # DS1's NXsample went somewhere harmless instead
            assert f["entry/sample_src/name"][()].decode() == "usaxs_sample"
            # and the rest propagated normally
            assert "metadata" in f["entry"]
            assert "metadata_saxs" in f["entry"]
