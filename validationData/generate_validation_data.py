#!/usr/bin/env python3
"""Generate the pyIrena / Irena synthetic validation dataset.

    python3 validationData/generate_validation_data.py

Writes, for every dataset defined in :mod:`_spec`:

  <folder>/<name>.dat        3-column ASCII  Q, I, dI   (noisy "measurement")
  <folder>/<name>_ideal.dat  2-column ASCII  Q, I       (exact, noise-free model)
  <folder>/<name>.h5         NXcanSAS HDF5 of the noisy data, with the exact
                             curve and the full ground truth stored alongside
                             under entry/ground_truth/

plus ``ground_truth.json`` (machine readable) and ``README.md`` (the manifest).

Everything is deterministic: the pseudo-random noise seed for each dataset is
derived from its name via SHA-256, so regenerating reproduces the files
byte-for-byte and adding a new dataset never changes an existing one.

The scattering models live in :mod:`_models` and are written from the published
literature **without importing pyIrena**, so recovering the parameters below
tests pyIrena's mathematics rather than its self-consistency.
"""

from __future__ import annotations

import hashlib
import json
import sys
from datetime import date
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import _spec as S  # noqa: E402

# ---------------------------------------------------------------------------
# Registry: (builder function, subfolder)
# ---------------------------------------------------------------------------

DATASETS = [
    (S.sizes_sphere_lognormal,          "sizes"),
    (S.sizes_sphere_bimodal,            "sizes"),
    (S.sizes_sphere_flat_background,    "sizes"),
    (S.sizes_spheroid_lognormal,        "sizes"),

    (S.unified_one_level,               "unified"),
    (S.unified_two_level,               "unified"),
    (S.unified_one_level_slitsmeared,   "unified"),

    (S.modeling_sizedist_plus_unified,  "modeling"),
    (S.modeling_sizedist_plus_peak,     "modeling"),
    (S.modeling_hardsphere,             "modeling"),
    (S.modeling_mass_fractal,           "modeling"),
    (S.modeling_surface_fractal,        "modeling"),
    (S.modeling_guinier_porod,          "modeling"),

    (S.simple_guinier,                  "simple_fits"),
    (S.simple_guinier_rod,              "simple_fits"),
    (S.simple_guinier_sheet,            "simple_fits"),
    (S.simple_porod,                    "simple_fits"),
    (S.simple_power_law,                "simple_fits"),
    (S.simple_sphere,                   "simple_fits"),
    (S.simple_debye_chain,              "simple_fits"),
    (S.simple_debye_bueche,             "simple_fits"),

    (S.waxs_three_gauss,                "waxs"),
    (S.waxs_pseudovoigt,                "waxs"),

    (S.merge_usaxs,                     "merge"),
    (S.merge_saxs,                      "merge"),

    (S.manipulate_reference,            "manipulation"),
    (S.manipulate_input,                "manipulation"),
]

# ---------------------------------------------------------------------------
# Noise
# ---------------------------------------------------------------------------

REL_SYSTEMATIC = 0.010      # 1 % point-to-point systematic term
TARGET_REL_AT_P5 = 0.060    # counting term gives 6 % at the 5th-percentile I


def seed_for(name: str) -> int:
    """Deterministic per-dataset seed derived from the dataset name."""
    return int(hashlib.sha256(name.encode()).hexdigest()[:8], 16)


def add_noise(name, I_ideal):
    """Counting-statistics + systematic noise.

        sigma(Q) = sqrt( I/kappa + (rel*I)^2 )
        I_obs    = I_ideal + sigma * N(0,1)

    ``kappa`` (an effective counts-per-unit-intensity monitor level) is chosen
    per dataset so that the counting term equals TARGET_REL_AT_P5 at the 5th
    percentile of the ideal intensity.  Draws that would make I_obs <= 0 are
    redrawn (their count is reported in the metadata); this only ever happens
    at the very bottom of a sphere form-factor minimum.
    """
    rng = np.random.default_rng(seed_for(name))
    i_ref = float(np.percentile(I_ideal, 5))
    kappa = 1.0 / (TARGET_REL_AT_P5 ** 2 * i_ref)
    sigma = np.sqrt(I_ideal / kappa + (REL_SYSTEMATIC * I_ideal) ** 2)

    I_obs = I_ideal + sigma * rng.standard_normal(I_ideal.size)
    n_redrawn = 0
    for _ in range(200):
        bad = I_obs <= 0
        if not bad.any():
            break
        n_redrawn += int(bad.sum())
        I_obs[bad] = I_ideal[bad] + sigma[bad] * rng.standard_normal(int(bad.sum()))
    I_obs = np.maximum(I_obs, 1e-12)

    info = dict(
        model="sigma = sqrt(I/kappa + (rel*I)^2); I_obs = I + sigma*N(0,1)",
        rel_systematic=REL_SYSTEMATIC,
        target_rel_at_5th_percentile=TARGET_REL_AT_P5,
        kappa=kappa,
        seed=seed_for(name),
        median_relative_error=float(np.median(sigma / I_ideal)),
        max_relative_error=float(np.max(sigma / I_ideal)),
        n_points_redrawn=n_redrawn,
    )
    return I_obs, sigma, info


# ---------------------------------------------------------------------------
# Writers
# ---------------------------------------------------------------------------

def _header_lines(meta):
    out = [
        f"pyIrena / Irena validation dataset: {meta['name']}",
        f"tool          : {meta['tool']}",
        f"model         : {meta['model']}",
        f"equation      : {meta['equation']}",
        "",
        "GROUND-TRUTH PARAMETERS",
    ]
    for p in meta["parameters"]:
        v = p["value"]
        vs = f"{v:.6g}" if isinstance(v, (int, float)) else str(v)
        out.append(f"  {p['name']:<26s} = {vs:<14s} [{p['unit']}]"
                   + (f"   {p['note']}" if p["note"] else ""))
    if meta["derived"]:
        out.append("")
        out.append("DERIVED QUANTITIES (not fitted; for cross-checking)")
        for p in meta["derived"]:
            v = p["value"]
            vs = f"{v:.6g}" if isinstance(v, (int, float)) else str(v)
            out.append(f"  {p['name']:<26s} = {vs:<14s} [{p['unit']}]"
                       + (f"   {p['note']}" if p["note"] else ""))
    n = meta["noise"]
    out += [
        "",
        "NOISE",
        f"  {n['model']}",
        f"  rel systematic = {n['rel_systematic']:.4g}, kappa = {n['kappa']:.6g},"
        f" seed = {n['seed']}",
        f"  median rel. error = {n['median_relative_error']:.4g},"
        f" max = {n['max_relative_error']:.4g}",
        "",
        f"Q range       : {meta['q_min']:.6g} to {meta['q_max']:.6g} 1/A,"
        f" {meta['n_points']} points, {meta['q_spacing']}-spaced",
        "Units         : Q [1/A], I [1/cm], dI [1/cm]",
        f"Generated     : {meta['generated']} by generate_validation_data.py",
    ]
    if meta.get("slit_length"):
        out.append(f"Slit smearing : dQl = {meta['slit_length']:.6g} 1/A"
                   " (infinite-slit-length / Lake smearing)")
    return out


def write_ascii(path, meta, q, I, dI=None):
    lines = ["# " + ln for ln in _header_lines(meta)]
    lines.append("# " + ("-" * 70))
    if dI is None:
        lines.append("#   Q[1/A]                I[1/cm]")
        body = "\n".join(f"{a:.8e}  {b:.8e}" for a, b in zip(q, I))
    else:
        lines.append("#   Q[1/A]                I[1/cm]               dI[1/cm]")
        body = "\n".join(f"{a:.8e}  {b:.8e}  {c:.8e}" for a, b, c in zip(q, I, dI))
    path.write_text("\n".join(lines) + "\n" + body + "\n")


def write_nxcansas(path, meta, q, I, dI, I_ideal):
    import h5py
    with h5py.File(path, "w") as f:
        f.attrs["default"] = "entry"
        f.attrs["producer"] = "pyIrena validationData/generate_validation_data.py"

        e = f.create_group("entry")
        e.attrs["NX_class"] = "NXentry"
        e.attrs["canSAS_class"] = "SASentry"
        e.attrs["default"] = "sasdata"
        e.create_dataset("definition", data="NXcanSAS")
        e.create_dataset("title", data=meta["name"])
        e.create_dataset("run", data=meta["name"])

        d = e.create_group("sasdata")
        d.attrs["NX_class"] = "NXdata"
        d.attrs["canSAS_class"] = "SASdata"
        d.attrs["signal"] = "I"
        d.attrs["I_axes"] = "Q"
        d.attrs["Q_indices"] = np.array([0], dtype="i4")

        ds = d.create_dataset("Q", data=q)
        ds.attrs["units"] = "1/angstrom"
        if meta.get("slit_length"):
            ds.attrs["resolutions"] = "dQl"
        ds = d.create_dataset("I", data=I)
        ds.attrs["units"] = "1/cm"
        ds.attrs["uncertainties"] = "Idev"
        ds = d.create_dataset("Idev", data=dI)
        ds.attrs["units"] = "1/cm"
        if meta.get("slit_length"):
            ds = d.create_dataset("dQl", data=float(meta["slit_length"]))
            ds.attrs["units"] = "1/angstrom"

        s = e.create_group("sample")
        s.attrs["NX_class"] = "NXsample"
        s.attrs["canSAS_class"] = "SASsample"
        s.create_dataset("name", data=meta["name"])
        s.create_dataset("thickness", data=1.0)

        # Non-standard, deliberately NOT an NXdata group so that pyIrena's
        # SASdata scan never offers it as a curve: the exact model and the
        # complete ground truth travel with the file.
        g = e.create_group("ground_truth")
        g.create_dataset("I_ideal", data=I_ideal)
        g.create_dataset("parameters_json", data=json.dumps(meta, indent=2))
        g.create_dataset("tool", data=meta["tool"])
        g.create_dataset("model", data=meta["model"])


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------

def build_all():
    today = date.today().isoformat()
    manifest = []
    produced: dict = {}
    for fn, folder in DATASETS:
        name = fn.__name__
        q, I_ideal, meta = fn()
        derive = meta.pop("derive_from", None)
        if derive is not None:
            # This dataset is an exact arithmetic transform of another one's
            # *measured* (noisy) curve, so that inverting the transform must
            # reproduce that curve to machine precision.  No noise is added
            # on top; the parent's uncertainties are scaled with the data.
            src, a, b = derive
            q_src, I_src, dI_src = produced[src]
            assert np.allclose(q, q_src), f"{name}: Q grid differs from {src}"
            I_ideal = a * I_src + b
            I_obs, dI = I_ideal, abs(a) * dI_src
            noise = dict(
                model=f"exact transform of {src}: I = {a:g} * I_{src} + {b:g}",
                rel_systematic=0.0, target_rel_at_5th_percentile=0.0,
                kappa=float("nan"), seed=seed_for(src),
                median_relative_error=float(np.median(dI / I_obs)),
                max_relative_error=float(np.max(dI / I_obs)),
                n_points_redrawn=0)
        else:
            I_obs, dI, noise = add_noise(name, I_ideal)

        spacing = "log" if np.std(np.diff(np.log(q))) < 1e-9 else "linear"
        meta.update(
            name=name, folder=folder, generated=today,
            q_min=float(q[0]), q_max=float(q[-1]), n_points=int(q.size),
            q_spacing=spacing, noise=noise,
            files=dict(ascii=f"{folder}/{name}.dat",
                       ascii_ideal=f"{folder}/{name}_ideal.dat",
                       nxcansas=f"{folder}/{name}.h5"),
        )

        out = HERE / folder
        out.mkdir(parents=True, exist_ok=True)
        write_ascii(out / f"{name}.dat", meta, q, I_obs, dI)
        write_ascii(out / f"{name}_ideal.dat", meta, q, I_ideal)
        write_nxcansas(out / f"{name}.h5", meta, q, I_obs, dI, I_ideal)

        produced[name] = (q, I_obs, dI)
        manifest.append(meta)
        print(f"  {folder}/{name}: {q.size} points, "
              f"I = {I_ideal.min():.4g} .. {I_ideal.max():.4g} 1/cm")
    return manifest


def _fmt(v):
    return f"{v:.6g}" if isinstance(v, (int, float)) and not isinstance(v, bool) else str(v)


README_PREAMBLE = """\
# Synthetic validation data for pyIrena and Irena

This folder holds **synthetic small-angle scattering data with exactly known
parameters**.  Its purpose is to make the correctness of pyIrena checkable, and
to let the same files be analysed with the Igor Pro *Irena* package so that the
two implementations can be compared parameter by parameter on identical input.

Every curve was produced by `generate_validation_data.py`, which computes the
scattering from independent implementations of the published model equations
(`_models.py`) and **does not import pyIrena**.  Recovering the tabulated
parameters is therefore a test of the analysis mathematics, not a
self-consistency check of a package inverting its own forward model.

## What is in each dataset

| file | content |
|---|---|
| `<name>.dat` | 3-column ASCII `Q  I  dI` with the noisy "measurement". The complete ground truth is repeated in the `#` header. |
| `<name>_ideal.dat` | 2-column ASCII `Q  I` — the exact, noise-free model. Use this to compare *forward models* between packages without any fitting. |
| `<name>.h5` | NXcanSAS HDF5 of the noisy data (`entry/sasdata`), with `entry/ground_truth/I_ideal` and `entry/ground_truth/parameters_json` alongside. |

`ground_truth.json` carries the same information machine-readably; it is what
`run_validation_report.py` reads when it fits every file and tabulates the
recovered parameters.

## Conventions

* `Q` in Å⁻¹, `I` in cm⁻¹ (absolute), radii and `Rg` in Å.
* `contrast` is (Δρ)² in units of 10²⁰ cm⁻⁴, matching Irena and pyIrena.
* Size distributions are quoted as **volume** distributions `P_V(r)`, normalised
  so that `∫P_V(r)dr` is the stated volume fraction.
* A log-normal `P_V(r)` is written with its **median** `R_median` and log-space
  standard deviation `sigma_log`; this is exactly Irena's `LogNormal` with
  `min_size = 0`, `mean_size = R_median`, `sdeviation = sigma_log`.

## Noise

`σ(Q) = sqrt(I/κ + (0.01·I)²)`, then `I_obs = I + σ·N(0,1)`, with `κ` chosen per
dataset so the counting term is 6 % at the 5th percentile of the intensity.
The random seed is `SHA-256(dataset name)[:8]`, so the files regenerate
byte-for-byte and adding a dataset never perturbs the others.

## Reproducing

```bash
python3 validationData/generate_validation_data.py    # regenerate the data
python3 validationData/run_validation_report.py       # fit it all, write the tables
```

## Tolerances

Recovered parameters are expected to agree with the ground truth to the
tolerances listed in `run_validation_report.py`.  Two systematic effects are
known and deliberately not hidden:

1. **Grid truncation.** Irena and pyIrena integrate a modelled size
   distribution over a finite radius grid that omits 1 % of the cumulative
   distribution at each tail, so a recovered volume fraction or `scale` may sit
   ~1–2 % below the generated value. The generator integrates over ±6σ.
2. **Regularised inversions.** The Size Distribution tool solves an
   ill-conditioned inverse problem; the *shape* it returns is smoothed relative
   to the generating distribution. Integral quantities (volume fraction, mean
   and RMS radius) are the meaningful comparison, not bin-by-bin values.

---
"""


def write_readme(manifest):
    lines = [README_PREAMBLE, "## Datasets\n"]
    lines.append("| # | dataset | tool | model | Q range (1/Å) | points |")
    lines.append("|---|---|---|---|---|---|")
    for i, m in enumerate(manifest, 1):
        lines.append(f"| {i} | `{m['folder']}/{m['name']}` | {m['tool']} | "
                     f"{m['model']} | {m['q_min']:.4g} – {m['q_max']:.4g} | "
                     f"{m['n_points']} |")
    lines.append("")
    lines.append("## Ground truth, dataset by dataset\n")
    for i, m in enumerate(manifest, 1):
        lines.append(f"### {i}. `{m['name']}`\n")
        lines.append(f"**Tool:** {m['tool']}  ")
        lines.append(f"**Model:** {m['model']}  ")
        lines.append(f"**Equation:** `{m['equation']}`  ")
        lines.append(f"**Q:** {m['q_min']:.6g} – {m['q_max']:.6g} Å⁻¹, "
                     f"{m['n_points']} points, {m['q_spacing']}-spaced  ")
        if m.get("slit_length"):
            lines.append(f"**Slit length dQl:** {m['slit_length']:.6g} Å⁻¹  ")
        lines.append("")
        lines.append("| parameter | value | unit | note |")
        lines.append("|---|---|---|---|")
        for p in m["parameters"]:
            lines.append(f"| `{p['name']}` | {_fmt(p['value'])} | {p['unit']} | {p['note']} |")
        if m["derived"]:
            lines.append("")
            lines.append("Derived (not fitted):\n")
            lines.append("| quantity | value | unit | note |")
            lines.append("|---|---|---|---|")
            for p in m["derived"]:
                lines.append(f"| `{p['name']}` | {_fmt(p['value'])} | {p['unit']} | {p['note']} |")
        n = m["noise"]
        lines.append("")
        lines.append(f"Noise: κ = {n['kappa']:.6g}, seed = {n['seed']}, "
                     f"median relative error {100*n['median_relative_error']:.2f} %, "
                     f"max {100*n['max_relative_error']:.2f} %"
                     + (f", {n['n_points_redrawn']} point(s) redrawn"
                        if n["n_points_redrawn"] else "") + ".")
        lines.append("")
    (HERE / "README.md").write_text("\n".join(lines) + "\n")


def main():
    print("Generating pyIrena/Irena validation data ...")
    manifest = build_all()
    (HERE / "ground_truth.json").write_text(json.dumps(manifest, indent=2) + "\n")
    write_readme(manifest)
    print(f"\n{len(manifest)} datasets written to {HERE}")
    print("  ground_truth.json, README.md updated")


if __name__ == "__main__":
    main()
