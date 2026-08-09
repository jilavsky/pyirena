"""
Source-level guard for the panel state-method convention.

Setup save/restore is the feature every panel must have and the one most easily
broken by a rename: the same two operations had grown six names
(``save_state`` / ``_save_state`` / ``load_state`` / ``_load_state`` /
``_restore_state`` / ``_apply_state``), so a developer reading one panel learned
nothing about the next, and a shared helper could not call them.

The convention, also stated in ``docs/developer_adding_features.md``:

* **panels** — public ``save_state()`` / ``load_state()``; the private halves
  ``_collect_state()`` (returns a dict) and ``_apply_state(state)`` where the
  panel is large enough to want them;
* **embedded components** driven by a parent panel (the Diffraction Lines tab
  inside WAXS) — public ``collect_state()`` / ``apply_state()``, because the
  parent, not the component, owns the persistence lifecycle.

Source text only: no Qt, no display, so this runs in every CI job.
"""

from __future__ import annotations

import re
from pathlib import Path

GUI_DIR = Path(__file__).resolve().parents[1] / "gui"
REPO = GUI_DIR.parent.parent

#: Top-level panels that own a StateManager section.
PANELS = [
    "pyirena/gui/unified_fit.py",
    "pyirena/gui/sizes_panel.py",
    "pyirena/gui/simple_fits_panel.py",
    "pyirena/gui/modeling_panel.py",
    "pyirena/gui/waxs_peakfit_panel.py",
    "pyirena/gui/fractals_panel.py",
    "pyirena/gui/saxs_morph_panel.py",
    "pyirena/gui/contrast_panel.py",
    "pyirena/gui/data_manipulation_panel.py",
    "pyirena/gui/data_merge_panel.py",
    "pyirena/gui/hdf5viewer/main_window.py",
]

#: Components whose parent drives their state.
COMPONENTS = ["pyirena/gui/diffraction_lines_panel.py"]

#: Names this convention replaced.  None may come back.
RETIRED = ("_save_state", "_load_state", "_restore_state", "_get_current_state")


def _defs(path: str) -> set[str]:
    text = (REPO / path).read_text(encoding="utf-8")
    return set(re.findall(r"^\s+def (\w+)\(", text, re.MULTILINE))


def test_every_panel_exposes_the_public_pair():
    missing = []
    for panel in PANELS:
        names = _defs(panel)
        for required in ("save_state", "load_state"):
            if required not in names:
                missing.append(f"{panel}: {required}()")
    assert not missing, (
        "CLAUDE.md requires every panel to expose save_state()/load_state() so "
        "setup save/restore works:\n  " + "\n  ".join(missing)
    )


def test_components_expose_the_collect_apply_pair():
    missing = []
    for component in COMPONENTS:
        names = _defs(component)
        for required in ("collect_state", "apply_state"):
            if required not in names:
                missing.append(f"{component}: {required}()")
    assert not missing, (
        "An embedded component must let its parent collect and apply its "
        "state:\n  " + "\n  ".join(missing)
    )


def test_retired_state_method_names_do_not_come_back():
    offenders = []
    for path in sorted(GUI_DIR.rglob("*.py")):
        text = path.read_text(encoding="utf-8")
        for name in RETIRED:
            if re.search(rf"^\s+def {name}\(", text, re.MULTILINE):
                offenders.append(f"{path.relative_to(REPO)}: {name}()")
    assert not offenders, (
        "These are the old names for save_state/load_state/_collect_state; "
        "use the current convention:\n  " + "\n  ".join(offenders)
    )


def test_unified_fit_keeps_its_documented_public_aliases():
    """``get_current_state`` / ``apply_state`` are quoted as the setup-state
    shape in pyirena/io/setup_config.py and the api/control layer, and agent
    scripts call them, so the rename kept them as aliases."""
    names = _defs("pyirena/gui/unified_fit.py")
    assert {"get_current_state", "apply_state"} <= names
    assert {"_collect_state", "_apply_state"} <= names
