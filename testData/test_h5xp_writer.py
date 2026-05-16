"""
testData/test_h5xp_writer.py — Generate a dummy h5xp file to verify the
structure before working with real data.

Run from the repo root:
    python testData/test_h5xp_writer.py

The output file (dummy_output.h5xp) can be opened directly in Igor Pro 9+ to
confirm waves load correctly.  Expected Igor folder layout:

    root:
    └── Packed Data:
        ├── SAXS:
        │   ├── sample_001:
        │   │   Q, R, S, dQ          (I(Q) data waves)
        │   │   UnifiedFitIntensity, UnifiedFitIntensity_X  (model waves)
        │   └── sample_002:
        │       Q, R, S              (no dQ for this sample)
        └── Results:
            └── Unified Fit:
                SampleNames, Rg_L1, G_L1, chi_squared
"""

import sys
from pathlib import Path

# Allow running from testData/ or from repo root
sys.path.insert(0, str(Path(__file__).parent.parent))

import numpy as np

from pyirena.io.h5xp_writer import (
    create_h5xp,
    make_wave_note,
    write_iq_data,
    write_result_wave,
    write_results_table,
)
from pyirena.io.igor_names import TOOL_CROSS_REF, SIMPLE_FIT_MODEL_WAVE

OUT = Path(__file__).parent / "dummy_output.h5xp"


def make_dummy_iq(n: int = 150, q_min: float = 0.001, q_max: float = 0.3,
                  rg: float = 50.0, G: float = 1.0) -> tuple:
    """Generate Guinier-law I(Q) with Porod tail + 5 % Gaussian noise."""
    q = np.geomspace(q_min, q_max, n)
    I = G * np.exp(-(q * rg) ** 2 / 3.0) + 1e-5 * q ** -4
    noise = 0.05 * I * np.abs(np.random.default_rng(42).standard_normal(n))
    dI = 0.05 * I
    dQ = 0.02 * q
    return q, I + noise, dI, dQ


def make_unified_model(q: np.ndarray, rg: float, G: float) -> np.ndarray:
    """Simple Guinier model used as a stand-in for a unified fit."""
    return G * np.exp(-(q * rg) ** 2 / 3.0)


def main() -> None:
    print(f"Writing dummy h5xp -> {OUT}")

    with create_h5xp(OUT, overwrite=True) as f:

        # ── Sample 1 ──────────────────────────────────────────────────────
        rg1, G1 = 45.0, 0.85
        q1, I1, dI1, dQ1 = make_dummy_iq(rg=rg1, G=G1)

        meta1 = {
            "SampleName": "sample_001",
            "Description": "Dummy silica nanoparticles in buffer",
            "StartTime": "2026-05-14T10:00:00",
        }
        write_iq_data(f, "sample_001", q1, I1, error=dI1, dq=dQ1,
                      wave_note=meta1, category="SAXS")

        model_q1 = np.geomspace(q1[0], q1[-1], 200)
        model_I1 = make_unified_model(model_q1, rg1, G1)
        write_result_wave(
            f, "sample_001", "UnifiedFitIntensity",
            model_q1, model_I1,
            params={
                "Rg_L1":      rg1,
                "G_L1":       G1,
                "B_L1":       1e-5,
                "P_L1":       4.0,
                "chi_squared": 0.0031,
                "background":  0.0002,
            },
            category="SAXS",
        )

        # ── Sample 2 (no dQ, different Rg) ────────────────────────────────
        rg2, G2 = 62.0, 1.10
        q2, I2, dI2, _ = make_dummy_iq(rg=rg2, G=G2)

        meta2 = {
            "SampleName": "sample_002",
            "Description": "Dummy sample at higher concentration",
            "StartTime": "2026-05-14T11:30:00",
        }
        write_iq_data(f, "sample_002", q2, I2, error=dI2,
                      wave_note=meta2, category="SAXS")

        model_q2 = np.geomspace(q2[0], q2[-1], 200)
        model_I2 = make_unified_model(model_q2, rg2, G2)
        write_result_wave(
            f, "sample_002", "UnifiedFitIntensity",
            model_q2, model_I2,
            params={
                "Rg_L1":      rg2,
                "G_L1":       G2,
                "B_L1":       1.2e-5,
                "P_L1":       3.9,
                "chi_squared": 0.0045,
                "background":  0.0003,
            },
            category="SAXS",
        )

        # ── Results table (Rg and G across both samples) ───────────────────
        samples = ["sample_001", "sample_002"]
        write_results_table(f, "Unified Fit", "Rg_L1",
                             [rg1, rg2], samples, units="Angstrom")
        write_results_table(f, "Unified Fit", "G_L1",
                             [G1, G2], samples, units="")
        write_results_table(f, "Unified Fit", "chi_squared",
                             [0.0031, 0.0045], samples)

    # ── Verify structure ──────────────────────────────────────────────────
    print("\nFile structure:")
    import h5py
    with h5py.File(OUT, "r") as f:
        def _show(name, obj):
            indent = "  " * name.count("/")
            attrs = {k: obj.attrs[k] for k in obj.attrs}
            if isinstance(obj, h5py.Dataset):
                note = attrs.get("IGORWaveNote", b"")
                if isinstance(note, (bytes, np.bytes_)):
                    note = note.decode("utf-8")
                note_preview = note[:60] + "…" if len(note) > 60 else note
                print(f"  {indent}{name}: shape={obj.shape} dtype={obj.dtype}")
                if note_preview:
                    print(f"  {indent}  note: {note_preview}")
            else:
                print(f"  {indent}{name}/")
        f.visititems(_show)

    print(f"\nDone.  Open {OUT.name} in Igor Pro 9+ to verify.")
    print("\nIgor wave name cross-reference (unified_fit):")
    for w in TOOL_CROSS_REF["unified_fit"]["waves"]:
        print(f"  {w['igor_name']:30s}  <- {w['pyirena_x']} / {w['pyirena_y']}")
    print("\nSimple fits model -> Igor wave name:")
    for model, igor in SIMPLE_FIT_MODEL_WAVE.items():
        print(f"  {model:25s} -> {igor}")


if __name__ == "__main__":
    main()
