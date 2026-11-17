"""Per-package part of the documentation screenshot machinery.

Copied from jilavsky-docs/templates/package-docs; **this is the only file in
scripts/doc_shots/ that differs between packages** (like the PER-PACKAGE
block of docs/conf.py). pyIrena adds ``VALIDATION``/``TEST_DATA`` and
:func:`wait_jobs` for its background fits.

``isolate_settings`` runs before anything imports Qt or the package, and
again before every shot, so each shot starts from the coded defaults and
nothing is written to the user's real settings.
"""

from __future__ import annotations

import os
from pathlib import Path

from ._common import REPO

#: Tracked validation datasets (the default for Context.dataset()).
VALIDATION: Path = REPO / "validationData"
#: Real measured data, for tools that need it (``ctx.dataset(..., root=TEST_DATA)``).
TEST_DATA: Path = REPO / "testData"

#: Default folder for Context.dataset(): tracked input files for the shots.
DATA_ROOT: Path = VALIDATION

#: Callout circles: an orange that reads on white and on plots.
CALLOUT_FILL = "#c2410c"
CALLOUT_RADIUS = 12


def isolate_settings(scratch: Path) -> None:
    """Point every settings/state location of the package into *scratch*.

    Find every place the package reads or writes user state: environment
    variables it honours, ``Path.home()`` / ``~`` (set ``HOME`` and
    ``USERPROFILE``), ``QSettings``, ``QStandardPaths`` (on macOS Qt ignores
    ``HOME`` for these: replace the method that asks), and the keychain.
    A shot must never read or write the developer's real ones; check
    modification times before and after a run. Examples: pyIrena
    ``PYIRENA_STATE_DIR``; AIDA ``AIDA_HOME`` and ``HOME``; Bernardyn ``HOME``
    plus a replaced template-folder method.
    """
    os.environ["PYIRENA_STATE_DIR"] = str(scratch / "state")


def make_app():
    """Create the offscreen QApplication as the real app sets it up.

    Same organisation/application names, light theme, window-geometry
    persistence off, and the app's default font applied once here (otherwise
    a shot that changes the font size changes every later image). A package
    without a theme sets ``app.setStyle("Fusion")`` so images look the same on
    every platform. ``refuse_message_boxes()`` makes a stray modal box fail
    the shot instead of hanging the run.
    """
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    os.environ.pop("PYIRENA_NATIVE_THEME", None)
    from pyirena.gui import window_state

    window_state.GEOMETRY_PERSISTENCE_ENABLED = False
    from pyirena.gui._qt import QApplication
    from pyirena.gui.theme import apply_theme

    from ._common import refuse_message_boxes

    app = QApplication.instance() or QApplication(["make_doc_screenshots"])
    apply_theme(app)
    refuse_message_boxes()
    return app


def wait_jobs(ctx, panel, timeout: float = 300.0) -> None:
    """Wait for a panel's background fit/MC job (``FitJobRunner``) to finish."""
    jobs = getattr(panel, "_jobs", None)
    ctx.settle(5)
    if jobs is not None:
        ctx.wait_until(lambda: not jobs.busy, timeout)
