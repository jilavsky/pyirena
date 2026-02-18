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
| NXcanSAS file format | [docs/NXcanSAS_UnifiedFit_Format.md](docs/NXcanSAS_UnifiedFit_Format.md) |
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

- **Unified Fit Model** — Beaucage hierarchical fit, 1–5 structural levels
- **Interactive GUI** — load data, adjust parameters, fit, store results
- **NXcanSAS I/O** — read/write HDF5 files in NXcanSAS format
- **Local fits** — Guinier and Porod region fits with cursor selection

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
