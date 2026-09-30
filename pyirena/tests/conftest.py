"""Test-wide isolation from the developer's own pyIrena settings.

Constructing a GUI panel constructs a :class:`StateManager`, and closing one
saves its window geometry. Both go to ``~/.pyirena`` unless told otherwise, so
running the suite used to rewrite the developer's real ``state.json`` and
``window_geometry.json``.

That was not merely untidy. A panel built in a test is never *shown*, so its
splitter reports layout hints rather than real pane widths — the right pane
collapses to its minimum. Saving those and rescaling them on the next launch
faithfully reproduced the ratio, so after running the tests the Unified Fit,
Size Distribution, Modeling and Simple Fits control panels opened at about
two-thirds of the window width. The panels were fine; the saved state was not.

Pointing ``$PYIRENA_STATE_DIR`` at a temporary directory for the session fixes
it at the root: the tests cannot reach the real settings at all. The guard in
``gui/window_state.save_window_state`` — which refuses to record geometry for
a window that was never shown — is the second half, and the half that protects
anyone embedding pyIrena rather than testing it.
"""

from __future__ import annotations

import os
import tempfile
from pathlib import Path

import pytest

from pyirena.state.state_manager import STATE_DIR_ENV_VAR


@pytest.fixture(scope="session", autouse=True)
def isolated_pyirena_state_dir():
    """Redirect every per-user pyIrena file into a throwaway directory.

    Session-scoped and autouse: it has to be in place before the first
    ``StateManager`` is built, and no test should have to remember it.
    """
    previous = os.environ.get(STATE_DIR_ENV_VAR)
    with tempfile.TemporaryDirectory(prefix="pyirena-test-state-") as tmp:
        os.environ[STATE_DIR_ENV_VAR] = tmp
        try:
            yield Path(tmp)
        finally:
            if previous is None:
                os.environ.pop(STATE_DIR_ENV_VAR, None)
            else:
                os.environ[STATE_DIR_ENV_VAR] = previous
