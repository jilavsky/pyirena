# pyIrena

**Coded by Claude from SAXS_IgorCode Irena. Planned, defined, debugged, and validated by Jan Ilavsky.**

Python tools for small-angle scattering (SAS) data analysis. A port of the Igor Pro
[Irena](https://usaxs.xray.aps.anl.gov/software/irena) package. Includes interactive
GUI tools for fitting, modeling, data merging, and visualization of SAXS/SANS/USAXS data.

**Current release: v1.2.1**

[![PyPI version](https://img.shields.io/pypi/v/pyirena.svg)](https://pypi.org/project/pyirena/)
[![Python Version](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

---

## Installation

**From PyPI (recommended):**

```bash
pip install pyirena[gui]
```
on macos with zsh:
```zsh
pip install "pyirena[gui]"
```
This installs pyirena with all GUI dependencies (PySide6, pyqtgraph, etc.). For the
core library only (no GUI), use `pip install pyirena`.

**Trying a pre-release (beta) version:**

Beta releases (e.g. `1.2.0b2`) are published to PyPI ahead of a stable release
for early testing. `pip` does not install pre-releases by default, so pass
`--pre` explicitly:

```bash
pip install --pre "pyirena[gui]"
```

To install a specific version rather than whatever is newest:

```bash
pip install "pyirena[gui]==1.2.1"
```

Pre-releases may contain incomplete or breaking changes — see
[CHANGELOG.md](CHANGELOG.md) for what's new, and please report issues at
https://github.com/jilavsky/pyirena/issues.

**From source (for development):**

```bash
git clone https://github.com/jilavsky/pyirena.git
cd pyirena
pip install -e ".[gui]"
```

**With conda:**

```bash
git clone https://github.com/jilavsky/pyirena.git
cd pyirena
conda env create -f environment.yml
conda activate pyirena
```

`environment.yml` uses conda for Python + scientific stack and pip (inside
conda) for the entire Qt6/GUI stack — this split is deliberate. Do **not**
run `conda install pyside6` in this env, and use Python 3.10 – 3.13
(Python 3.14 wheels are still flaky for compiled bindings).

See [docs/getting-started/installation.md](docs/getting-started/installation.md) for full details, troubleshooting, and platform notes.

---

## Running the GUI

```bash
pyirena-gui
```

This launches the Data Selector, the main entry point for all analysis tools.
See [docs/getting-started/quickstart.md](docs/getting-started/quickstart.md) for a walkthrough.

Individual tools can also be launched directly:

| Command | Tool |
|---------|------|
| `pyirena-gui` | Data Selector (main entry point) |
| `pyirena-viewer` | Data Explorer |
| `pyirena-modeling` | Modeling tool (standalone) |
| `pyirena-datamerge` | Data Merge tool (standalone) |
| `pyirena-contrast` | Scattering Contrast Calculator |
| `pyirena-datamanip` | Data Manipulation tool (standalone) |
| `pyirena-saxsmorph` | SAXS Morph tool (standalone) |
| `pyirena-mcp` | MCP server for AI assistants (stdio) — see [AI integration](#ai-integration-mcp) |
| `pyirena-zmq` | ZMQ service for remote/beamline orchestration ([docs](docs/reference/zmq-service.md)) |
| `pyirena-doctor` | Install/startup diagnostics |

`[project.scripts]` in `pyproject.toml` is the authoritative list.

---

## Analysis Tools

### Modeling
Parametric forward-modeling of SAS data. Combine up to 10 populations, each of which
can be a Size Distribution, Unified Fit Level, or Diffraction Peak. Fits the combined
model to experimental data using least-squares optimization.

- **5 distribution functions**: Gaussian, LogNormal, LSW, Schulz-Zimm, Ardell
- **10 form factors**: Sphere, Spheroid, Cylinder (aspect ratio / fixed length),
  Core-Shell Sphere (by core R / shell t / total R), Core-Shell-Shell Sphere,
  Core-Shell Spheroid (by core R / total R)
- **2 structure factors**: Born-Green interferences, Hard Sphere (Percus-Yevick)
- Monte Carlo uncertainty estimation
- [Modeling GUI guide](docs/tools/modeling.md)

### Unified Fit
Beaucage hierarchical model (1995, 1996) with 1-5 structural levels, each combining
Guinier and power-law contributions. Optional Born-Green correlation function.

- [Unified Fit GUI guide](docs/tools/unified-fit.md)

### Size Distribution
Indirect Fourier transform to recover particle size distributions from SAS data.
Four inversion methods: MaxEnt, Regularization, TNNLS, and Monte Carlo.

- [Size Distribution tool](docs/tools/size-distribution.md)
- [Size Distribution methods](docs/methods/size-distribution-methods.md)

### Simple Fits
15 direct analytical models: Guinier (sphere, rod, sheet), Porod, Power Law,
Debye polymer chain, Sphere, Spheroid, Debye-Bueche, Teubner-Strey,
Benedetti-Ciccariello, Hermans, Hybrid Hermans, Unified Born-Green, and the
Invariant calculation. Each fit has linearization plots and Monte Carlo
uncertainty estimation.

- [Simple Fits GUI guide](docs/tools/simple-fits.md)

### WAXS Peak Fit
Fit diffraction peaks (Gaussian, Lorentzian, Pseudo-Voigt, Log-Normal) on linear
I vs Q scale. Auto-detect peaks via Savitzky-Golay + `scipy.signal.find_peaks`.
Simultaneous background fitting (constant, linear, cubic, or 5th-order polynomial).

- [WAXS Peak Fit GUI guide](docs/tools/waxs-peak-fit.md)

### Data Merge
Merge two SAS datasets (e.g., SAXS + WAXS) onto a common Q scale. Optimizes scale
factor, flat background, and optional Q-shift using Nelder-Mead.

- [Data Merge GUI guide](docs/tools/data-merge.md)

### Data Explorer
Browse NXcanSAS HDF5 files, inspect raw data and analysis results, extract and plot
datasets, and export to Igor Pro h5xp format. Supports all pyIrena result types
(Unified Fit, Sizes, Simple Fits, WAXS Peaks, Modeling).

- [Data Explorer guide](docs/tools/data-explorer.md)

### Import Igor Experiment
Open an Igor Pro packed experiment (`.pxp` or `.h5xp`) and export each
reduced USAXS/SAXS/WAXS sample as a stand-alone NXcanSAS `.h5` file.
Brings legacy and modern Igor data into pyIrena's analysis tools
without re-reducing from raw detector files.

- [Import Igor Experiment guide](docs/guides/igor-pxp-import.md)

### Scattering Contrast Calculator
Look up X-ray and neutron scattering length densities for materials by chemical
formula. Computes contrast (Delta-rho-squared) between two materials.

- [Contrast Calculator guide](docs/tools/scattering-contrast.md)

### Data Selector
Central GUI panel for managing data files and launching analysis tools. Load HDF5
files, select datasets, tabulate results across files, and generate reports.

---

## Batch Scripting API

All analysis tools can be run headlessly from Python scripts or JSON configuration files:

```python
from pyirena.batch import fit_pyirena

results = fit_pyirena(
    data_file='sample.h5',
    config_file='pyirena_config.json',
)
```

Individual functions: `fit_unified`, `fit_sizes`, `fit_simple`, `fit_waxs`,
`fit_modeling`, `merge_data`.

See [docs/guides/batch-api.md](docs/guides/batch-api.md) for the full scripting guide.

---

## NXcanSAS I/O

All data and results are stored in HDF5 files using the
[NXcanSAS](https://www.nexusformat.org/NXcanSAS.html) format. Fit results are saved
alongside raw data, making files self-contained and shareable.

- [NXcanSAS format details](docs/reference/hdf5-nxcansas-structure.md)

---

## AI integration (MCP)

pyirena ships a [Model Context Protocol](https://modelcontextprotocol.io/) server
that exposes its HDF5 readers, parameter aggregation, and headless plotting to any
MCP client — Claude Desktop, Claude Code, custom AI agents.

```bash
pip install pyirena[mcp]
```

Add to `claude_desktop_config.json` (or your client's MCP server config):

```json
{
  "mcpServers": {
    "pyirena": {
      "command": "pyirena-mcp",
      "env": {
        "PYIRENA_DATA_ROOT": "/path/to/your/data",
        "PYIRENA_PLOT_CACHE": "/tmp/pyirena-plots",
        "PYIRENA_MAX_ARRAY_POINTS": "500"
      }
    }
  }
}
```

The agent can then call tools like `pyirena_summarize_folder`,
`pyirena_tabulate_parameter`, `pyirena_plot_iq`, `pyirena_read_modeling` etc.
to answer questions about analysis results in plain language. The underlying
`pyirena.api` module is also usable as a regular Python library, with unprefixed
names (`api.summarize_folder`, `api.plot_iq`, …) — no MCP needed.

Full instructions, including configs for **Claude Desktop**, **Claude Code**,
and **AnythingLLM**:

- [docs/guides/ai-integration.md](docs/guides/ai-integration.md) — install & configure
- [docs/reference/ai-tools-reference.md](docs/reference/ai-tools-reference.md) — tool reference for AI agents (system-prompt material)

---

## Documentation

The documentation is published at
<https://jilavsky-docs.readthedocs.io/projects/pyirena/>; the Help buttons in the
GUI open the matching page there. The same pages, as Markdown, are in `docs/`:

| Topic | File |
|-------|------|
| Installation | [docs/getting-started/installation.md](docs/getting-started/installation.md) |
| Quick start (GUI) | [docs/getting-started/quickstart.md](docs/getting-started/quickstart.md) |
| Modeling GUI guide | [docs/tools/modeling.md](docs/tools/modeling.md) |
| Unified Fit GUI guide | [docs/tools/unified-fit.md](docs/tools/unified-fit.md) |
| Size Distribution tool | [docs/tools/size-distribution.md](docs/tools/size-distribution.md) |
| Size Distribution methods | [docs/methods/size-distribution-methods.md](docs/methods/size-distribution-methods.md) |
| Simple Fits GUI guide | [docs/tools/simple-fits.md](docs/tools/simple-fits.md) |
| WAXS Peak Fit GUI guide | [docs/tools/waxs-peak-fit.md](docs/tools/waxs-peak-fit.md) |
| Data Merge GUI guide | [docs/tools/data-merge.md](docs/tools/data-merge.md) |
| Data Explorer guide | [docs/tools/data-explorer.md](docs/tools/data-explorer.md) |
| Import Igor Experiment (.pxp / .h5xp → NeXus) | [docs/guides/igor-pxp-import.md](docs/guides/igor-pxp-import.md) |
| Contrast Calculator guide | [docs/tools/scattering-contrast.md](docs/tools/scattering-contrast.md) |
| NXcanSAS file format | [docs/reference/hdf5-nxcansas-structure.md](docs/reference/hdf5-nxcansas-structure.md) |
| Batch fitting API | [docs/guides/batch-api.md](docs/guides/batch-api.md) |
| Usage guide (scripting) | [docs/guides/scripting-unified-fit.md](docs/guides/scripting-unified-fit.md) |
| AI integration (MCP) — install & configure | [docs/guides/ai-integration.md](docs/guides/ai-integration.md) |
| AI integration (MCP) — tool reference for agents | [docs/reference/ai-tools-reference.md](docs/reference/ai-tools-reference.md) |
| Developer: adding features (master checklist) | [docs/developer/adding-features.md](docs/developer/adding-features.md) |
| Developer: adding form factors | [docs/developer/adding-form-factors.md](docs/developer/adding-form-factors.md) |
| Developer: adding structure factors | [docs/developer/adding-structure-factors.md](docs/developer/adding-structure-factors.md) |
| Testing | [docs/developer/testing.md](docs/developer/testing.md) |
| Distribution / packaging | [docs/developer/distribution.md](docs/developer/distribution.md) |
| Contributing | [CONTRIBUTING.md](CONTRIBUTING.md) |
| Changelog | [CHANGELOG.md](CHANGELOG.md) |

---

## License

MIT -- see [LICENSE](LICENSE).

## Contact

- **Author**: Jan Ilavsky -- ilavsky@aps.anl.gov
- **Affiliation**: X-ray Science Division, Argonne National Laboratory
- **Issues**: https://github.com/jilavsky/pyirena/issues
