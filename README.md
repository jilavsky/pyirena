# pyIrena

> **Development version — not ready for production use.**
> Coded by Claude from SAXS_IgorCode Irena.

Python tools for small-angle scattering data analysis, centered on the **Unified Fit model** (Beaucage method) with an interactive GUI.

[![Python Version](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

---

## Documentation

| Topic | File |
|-------|------|
| Installation | [docs/installation.md](docs/installation.md) |
| Quick start (GUI) | [docs/gui_quickstart.md](docs/gui_quickstart.md) |
| Unified Fit GUI guide | [docs/unified_fit_gui.md](docs/unified_fit_gui.md) |
| Unified Fit features & parameters | [docs/unified_fit_features.md](docs/unified_fit_features.md) |
| Simple Fits GUI guide | [docs/simple_fits_gui.md](docs/simple_fits_gui.md) |
| WAXS Peak Fit GUI guide | [docs/waxs_peakfit_gui.md](docs/waxs_peakfit_gui.md) |
| NXcanSAS file format | [docs/NXcanSAS_UnifiedFit_Format.md](docs/NXcanSAS_UnifiedFit_Format.md) |
| Batch fitting API (scripting & automation) | [docs/batch_api.md](docs/batch_api.md) |
| Usage guide (scripting) | [docs/usage_guide.md](docs/usage_guide.md) |
| Testing | [docs/testing.md](docs/testing.md) |
| Distribution / packaging | [docs/distribution.md](docs/distribution.md) |
| Contributing | [CONTRIBUTING.md](CONTRIBUTING.md) |
| Changelog | [CHANGELOG.md](CHANGELOG.md) |

---

## Installation

```bash
git clone https://github.com/jilavsky/pyirena.git
cd pyirena
pip install -e ".[gui]"
```

This installs pyirena with GUI support (PySide6 + pyqtgraph).

---

## Running the GUI

```bash
pyirena-gui
```

See [docs/gui_quickstart.md](docs/gui_quickstart.md) for a full walkthrough.

---

## Overview

`pyIrena` is a Python port of the Igor Pro **Irena** package for SAS data analysis. Current capabilities:

- **Unified Fit Model** — Beaucage hierarchical fit, 1–5 structural levels, Born-Green correlation function
- **Size Distribution** — indirect Fourier transform with four inversion methods (MaxEnt, Regularization, TNNLS, Monte Carlo)
- **Simple Fits** — 13 direct analytical models (Guinier family, Porod, Sphere, Spheroid, Debye-Bueche, Treubner-Strey, and more) with linearization plots and Monte Carlo uncertainty
- **WAXS Peak Fit** — fit Gaussian, Lorentzian, Pseudo-Voigt, or log-normal diffraction peaks on a linear/linear intensity vs Q scale; auto-detect peaks via Savitzky-Golay + `scipy.signal.find_peaks`; fit background simultaneously (constant, linear, cubic, or 5th-order polynomial); interactive overlay of each peak + background on the main plot; results saved to `entry/waxs_peakfit_results` in NXcanSAS HDF5
- **Interactive GUI** — load data, adjust parameters, fit, store results; all tools share a common Data Selector
- **NXcanSAS I/O** — read/write HDF5 files in NXcanSAS format; all fit results stored alongside raw data
- **Batch scripting API** — `fit_unified`, `fit_sizes`, `fit_simple`, `fit_waxs`, `fit_pyirena`; headless fitting from JSON config files

---

## Mathematical Model

The Unified fit model (Beaucage 1995, 1996) combines Guinier and power-law contributions for each structural level:

```
I_i(q) = G_i · exp(-q²Rg_i²/3) + exp(-q²Rg_(i-1)²/3) · B_i · [erf(k·q·Rg_i/√6)³/q]^P_i
```

Optional correlation function (Born-Green approximation):

```
I_i(q) = I_i(q) / (1 + PACK · f(q, ETA))
```

References:
- Beaucage, G. (1995). *J. Appl. Cryst.* **28**, 717–728
- Beaucage, G. (1996). *J. Appl. Cryst.* **29**, 134–146

---

## License

MIT — see [LICENSE](LICENSE).

## Contact

- **Author**: Jan Ilavsky — ilavsky@aps.anl.gov
- **Issues**: https://github.com/jilavsky/pyirena/issues
