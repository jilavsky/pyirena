"""Documentation screenshots, one module per docs page.

Copied from jilavsky-docs/templates/package-docs. Each ``<page>.py`` module
(page name with ``-`` written as ``_``) registers its shots with
``@shot("<page>-<nn>-<what>")`` from :mod:`._common`. :func:`load_all`
imports every module not starting with ``_``, filling ``REGISTRY``.
See ``_example.py`` for the pattern.
"""

import importlib
import pkgutil

from ._common import REGISTRY  # noqa: F401  (re-exported)


def load_all() -> None:
    """Import every page module in this package so its shots register."""
    for info in sorted(pkgutil.iter_modules(__path__), key=lambda m: m.name):
        if not info.name.startswith("_"):
            importlib.import_module(f"{__name__}.{info.name}")
