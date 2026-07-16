"""
Tests for the Invariant calculation model in Simple Fits.

Analytic ground truth: for dilute monodisperse spheres on absolute scale,

    I(q) = phi * (drho)^2 * V * P(q)     [cm^-1, q in A^-1, V in A^3]

the Porod invariant is exact and morphology-independent:

    Q* = integral q^2 I(q) dq = 2*pi^2 * (drho)^2 * phi*(1-phi)

With contrast drho^2 expressed in 10^20 cm^-4 (Igor Irena convention) the
implementation must recover phi from synthetic data to ~1% (truncation).
"""

import numpy as np
import pytest

from pyirena.core.simple_fits import (
    SimpleFitModel,
    MODEL_REGISTRY,
    calculate_invariant,
)

# ── Synthetic dilute-sphere data (phi = 0.01, contrast = 100e20 cm^-4) ──────
R = 50.0                       # sphere radius [A]
PHI = 0.01                     # volume fraction
CONTRAST_1E20 = 100.0          # drho^2 in 10^20 cm^-4
DRHO2 = CONTRAST_1E20 * 1e20   # cm^-4
V = 4.0 / 3.0 * np.pi * R**3   # A^3

Q = np.logspace(-4, 0.0, 4000)
_u = Q * R
_P = (3.0 * (np.sin(_u) - _u * np.cos(_u)) / _u**3) ** 2
I_SPHERE = PHI * DRHO2 * V * 1e-24 * _P      # cm^-1

EXPECTED_QSTAR = 2.0 * np.pi**2 * DRHO2 * PHI * (1.0 - PHI)   # cm^-4


class TestCalculateInvariant:
    def test_registry_entry(self):
        entry = MODEL_REGISTRY["Invariant"]
        assert entry["calculation"] is True
        assert entry["complex_bg"] is True
        assert entry["linearization"] is None
        names = [p[0] for p in entry["params"]]
        assert names == ["Contrast"]

    def test_recovers_invariant_and_phi(self):
        res = calculate_invariant(Q, I_SPHERE, contrast=CONTRAST_1E20)
        assert res["success"]
        assert res["invariant"] == pytest.approx(EXPECTED_QSTAR, rel=0.02)
        assert res["volume_fraction"] == pytest.approx(PHI, rel=0.02)
        assert res["Q_integral"][0] == 0.0          # extended to Q=0
        assert np.all(np.diff(res["running_integral"]) >= -1e-30)

    def test_porod_tail_improves_truncation(self):
        # Truncate at moderate q so the tail matters more
        m = Q <= 0.3
        res_no = calculate_invariant(Q[m], I_SPHERE[m], contrast=CONTRAST_1E20)
        res_tail = calculate_invariant(Q[m], I_SPHERE[m],
                                       contrast=CONTRAST_1E20, porod_tail=True)
        err_no = abs(res_no["volume_fraction"] - PHI)
        err_tail = abs(res_tail["volume_fraction"] - PHI)
        assert res_tail["porod_tail"] > 0
        assert err_tail < err_no

    def test_flat_background_subtraction(self):
        flat = 0.5 * float(I_SPHERE[-1])
        res = calculate_invariant(Q, I_SPHERE + flat,
                                  contrast=CONTRAST_1E20, bg_flat=flat)
        assert res["volume_fraction"] == pytest.approx(PHI, rel=0.03)

    def test_unphysical_contrast_gives_nan_phi(self):
        # Contrast far too small -> phi*(1-phi) > 0.249 -> NaN + warning
        res = calculate_invariant(Q, I_SPHERE, contrast=0.01)
        assert np.isnan(res["volume_fraction"])
        assert "0.249" in res["warning"]

    def test_too_few_points_fails(self):
        res = calculate_invariant(Q[:3], I_SPHERE[:3], contrast=CONTRAST_1E20)
        assert not res["success"]


class TestModelRouting:
    def test_fit_delegates_to_calculation(self):
        m = SimpleFitModel()
        m.set_model("Invariant")
        assert m.is_calculation
        m.params["Contrast"] = CONTRAST_1E20
        result = m.fit(Q, I_SPHERE)
        assert result["success"]
        assert result["chi2"] is None
        assert result["residuals"] is None
        d = result["derived"]
        assert d["VolumeFraction"] == pytest.approx(PHI, rel=0.02)
        assert d["Invariant"] == pytest.approx(EXPECTED_QSTAR, rel=0.02)
        assert d["QmaxUsed"] > 0
        extra = result["extra_arrays"]
        assert extra["Q_integral"] is not None
        assert len(extra["running_integral"]) == len(extra["Q_integral"])

    def test_complex_bg_used_in_calculation(self):
        flat = 0.5 * float(I_SPHERE[-1])
        m = SimpleFitModel()
        m.set_model("Invariant")
        m.use_complex_bg = True
        m._reset_to_defaults()
        m.params["Contrast"] = CONTRAST_1E20
        m.params["BG_flat"] = flat
        m.params["BG_B"] = 0.0
        result = m.fit(Q, I_SPHERE + flat)
        assert result["derived"]["VolumeFraction"] == pytest.approx(PHI, rel=0.03)
        # I_model is the background curve
        assert result["I_model"] == pytest.approx(
            np.full_like(Q, flat), rel=1e-12)

    def test_serialization_round_trip_with_porod_tail(self):
        m = SimpleFitModel()
        m.set_model("Invariant")
        m.invariant_porod_tail = True
        d = m.to_dict()
        m2 = SimpleFitModel.from_dict(d)
        assert m2.model == "Invariant"
        assert m2.invariant_porod_tail is True

    def test_bg_prefit_serialization_round_trip(self):
        m = SimpleFitModel()
        m.set_model("Invariant")
        m.bg_prefit = {
            "enabled": True,
            "power_law": {"use": True, "q_min": 1e-4, "q_max": 5e-4,
                          "fit_P": True},
            "flat": {"use": True, "q_min": 0.5, "q_max": 1.0},
        }
        m2 = SimpleFitModel.from_dict(m.to_dict())
        assert m2.bg_prefit == m.bg_prefit
        # Old dicts without the key must load with the default
        d = m.to_dict()
        del d["bg_prefit"]
        m3 = SimpleFitModel.from_dict(d)
        assert m3.bg_prefit == {}


class TestBackgroundPrefitReplay:
    """Prefit replay: refit BG terms from saved Q windows before calculating."""

    # Flat background well above the sphere's high-Q signal (median sphere
    # intensity in the window is ~1e-4 cm^-1), so the flat window is truly
    # background-dominated — as it must be for a valid flat prefit.
    FLAT = 0.05
    FLAT_WIN = (0.5, 1.0)

    def _model(self, enabled=True):
        m = SimpleFitModel()
        m.set_model("Invariant")
        m.use_complex_bg = True
        m._reset_to_defaults()
        m.params["Contrast"] = CONTRAST_1E20
        # Deliberately wrong exported values — replay must fix them
        m.params["BG_flat"] = 123.0
        m.params["BG_B"] = 0.0
        m.bg_prefit = {
            "enabled": enabled,
            "flat": {"use": True, "q_min": self.FLAT_WIN[0],
                     "q_max": self.FLAT_WIN[1]},
        }
        return m

    def test_replay_refits_flat_and_recovers_phi(self):
        m = self._model(enabled=True)
        I_bg = I_SPHERE + self.FLAT
        applied = m.prefit_background(Q, I_bg)
        assert applied["BG_flat"] == pytest.approx(self.FLAT, rel=0.05)
        result = m.fit(Q, I_bg)
        assert result["derived"]["VolumeFraction"] == pytest.approx(PHI, rel=0.03)

    def test_replay_disabled_returns_empty(self):
        m = self._model(enabled=False)
        assert m.prefit_background(Q, I_SPHERE + self.FLAT) == {}
        # and the wrong exported flat stays in place
        assert m.params["BG_flat"] == 123.0

    def test_replay_requires_complex_bg(self):
        m = self._model(enabled=True)
        m.use_complex_bg = False
        assert m.prefit_background(Q, I_SPHERE + self.FLAT) == {}

    def test_replay_warns_on_empty_window(self):
        m = self._model(enabled=True)
        m.bg_prefit["flat"] = {"use": True, "q_min": 5.0, "q_max": 10.0}
        applied = m.prefit_background(Q, I_SPHERE + self.FLAT)
        assert "warning" in applied
        assert "BG_flat" not in applied

    def test_power_law_replay(self):
        # Add a q^-4 power-law background dominating the window (the sphere's
        # own Porod tail acts like an effective B of ~4e-5, so B_true must be
        # much larger for a clean single-component window) and refit B at
        # fixed P.
        B_true = 1e-2
        I_bg = I_SPHERE + B_true * Q**-4.0
        m = SimpleFitModel()
        m.set_model("Invariant")
        m.use_complex_bg = True
        m._reset_to_defaults()
        m.params["Contrast"] = CONTRAST_1E20
        m.params["BG_P"] = 4.0
        m.bg_prefit = {
            "enabled": True,
            "power_law": {"use": True, "q_min": 0.5, "q_max": 1.0,
                          "fit_P": False},
        }
        applied = m.prefit_background(Q, I_bg)
        assert applied["BG_B"] == pytest.approx(B_true, rel=0.1)

    def test_batch_replay_end_to_end(self, tmp_path):
        import h5py
        from pyirena.batch.simple import fit_simple

        fp = tmp_path / "sphere_bg.h5"
        with h5py.File(fp, "w") as f:
            e = f.create_group("entry")
            e.attrs["NX_class"] = "NXentry"
            s = e.create_group("sasdata")
            s.attrs["NX_class"] = "NXdata"
            s.attrs["signal"] = "I"
            s.attrs["I_axes"] = "Q"
            s.create_dataset("Q", data=Q)
            s.create_dataset("I", data=I_SPHERE + self.FLAT)

        config = {
            "model": "Invariant",
            "use_complex_bg": True,
            # Wrong exported flat — the per-file replay must correct it
            "params": {"Contrast": CONTRAST_1E20, "BG_B": 0.0,
                       "BG_P": 4.0, "BG_flat": 123.0},
            "bg_prefit": {
                "enabled": True,
                "flat": {"use": True, "q_min": self.FLAT_WIN[0],
                         "q_max": self.FLAT_WIN[1]},
            },
        }
        result = fit_simple(fp, config, verbose=False)
        assert result is not None and result["success"]
        assert result["derived"]["VolumeFraction"] == pytest.approx(PHI, rel=0.03)
        assert result["params"]["BG_flat"] == pytest.approx(self.FLAT, rel=0.05)


class TestHdf5RoundTrip:
    def test_save_and_load(self, tmp_path):
        import h5py
        from pyirena.io.nxcansas_simple_fits import (
            save_simple_fit_results, load_simple_fit_results,
        )

        m = SimpleFitModel()
        m.set_model("Invariant")
        m.params["Contrast"] = CONTRAST_1E20
        m.invariant_porod_tail = False
        result = m.fit(Q, I_SPHERE)

        fp = tmp_path / "invariant_test.h5"
        with h5py.File(fp, "w") as f:
            f.create_group("entry")

        save_simple_fit_results(
            filepath=fp, result=result, model_obj=m,
            intensity_data=I_SPHERE,
        )
        loaded = load_simple_fit_results(fp)
        assert loaded["model"] == "Invariant"
        assert loaded["derived"]["VolumeFraction"] == pytest.approx(PHI, rel=0.02)
        assert loaded["derived"]["Invariant"] == pytest.approx(
            EXPECTED_QSTAR, rel=0.02)
        assert loaded["Q_integral"] is not None
        assert loaded["running_integral"] is not None
        assert len(loaded["running_integral"]) == len(loaded["Q_integral"])
        # Calculation models must not write misleading fit-quality metrics
        assert loaded.get("fit_quality") in (None, {})


class TestBatchConfig:
    def test_fit_simple_with_invariant_config(self, tmp_path):
        import h5py
        from pyirena.batch.simple import fit_simple
        from pyirena.io.nxcansas_simple_fits import load_simple_fit_results

        # Minimal NXcanSAS file with the synthetic data
        fp = tmp_path / "sphere_data.h5"
        with h5py.File(fp, "w") as f:
            entry = f.create_group("entry")
            entry.attrs["NX_class"] = "NXentry"
            sas = entry.create_group("sasdata")
            sas.attrs["NX_class"] = "NXdata"
            sas.attrs["signal"] = "I"
            sas.attrs["I_axes"] = "Q"
            sas.create_dataset("Q", data=Q)
            sas.create_dataset("I", data=I_SPHERE)
            sas.create_dataset("Idev", data=0.01 * I_SPHERE)

        config = {
            "model": "Invariant",
            "params": {"Contrast": CONTRAST_1E20},
            "invariant_porod_tail": False,
        }
        result = fit_simple(fp, config, verbose=False)
        assert result is not None and result["success"]
        assert result["derived"]["VolumeFraction"] == pytest.approx(PHI, rel=0.02)

        loaded = load_simple_fit_results(fp)
        assert loaded["model"] == "Invariant"
        assert loaded["derived"]["VolumeFraction"] == pytest.approx(PHI, rel=0.02)
