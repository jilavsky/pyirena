"""Sphinx configuration for <PACKAGE> (an RTD subproject of jilavsky-docs).

Copied from jilavsky-docs/templates/package-docs. Change only the block
marked PER-PACKAGE; keep the rest identical across the family so the
copies can be diffed against the template.
"""

import os
import sys
from importlib import metadata
from pathlib import Path

# ---- PER-PACKAGE ----------------------------------------------------------
PACKAGE = "pyirena"  # import name
DIST = "pyirena"  # PyPI distribution name (AIDA: "aida-workbench")
project = "pyIrena"  # display name
GITHUB_REPO = "pyirena"  # github.com/jilavsky/<GITHUB_REPO>
SRC_DIR = ".."  # where PACKAGE lives relative to docs/ ("../src" for src layout)
# Heavy or GUI-only imports that autodoc must not try to import on RTD.
autodoc_mock_imports = ["PySide6", "pyqtgraph", "matplotlib", "OpenGL"]
# Old GitHub-linked files kept as one-line "moved" stubs for released Help
# buttons (see PLAN.md, "Transition"). Excluded from the Sphinx build.
LEGACY_GITHUB_STUBS: list[str] = [
    "unified_fit_gui.md",
    "sizes_methods.md",
    "modeling_gui.md",
    "simple_fits_gui.md",
    "fractals_gui.md",
    "waxs_peakfit_gui.md",
    "carbon_fit_gui.md",
    "saxs_morph_gui.md",
    "data_merge_gui.md",
    "data_manipulation_gui.md",
    "scattering_contrast_gui.md",
    "hdf5_viewer_gui.md",
    "GUI_README.md",
]
# ---------------------------------------------------------------------------

sys.path.insert(0, str((Path(__file__).parent / SRC_DIR).resolve()))
author = "Jan Ilavsky"
copyright = "2026, Jan Ilavsky"  # noqa: A001  (Sphinx name)
# Never bind the name `version` at module level: Sphinx reads it as a setting.
try:
    release = metadata.version(DIST)
except metadata.PackageNotFoundError:
    release = "unknown"
version = release

extensions = [
    "myst_parser",
    "sphinx_design",
    "sphinx_copybutton",
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.intersphinx",
]

myst_enable_extensions = ["colon_fence", "deflist", "dollarmath", "substitution"]
myst_heading_anchors = 3
napoleon_google_docstring = True
napoleon_numpy_docstring = True
autodoc_typehints = "description"

exclude_patterns = ["_build", "Thumbs.db", ".DS_Store", "**/_template*.md", *LEGACY_GITHUB_STUBS]

# ---- Family navbar: identical block in every family conf.py ---------------
# Version-less addresses: RTD redirects them to each project's default version.
# The home page has no releases, so it never has a "stable".
_RTD = "https://jilavsky-docs.readthedocs.io"
_IGOR = "https://saxs-igorcodedocs.readthedocs.io/en/stable/"
FAMILY = [
    {"name": "Home", "url": f"{_RTD}/"},
    {"name": "pyIrena", "url": f"{_RTD}/projects/pyirena/"},
    {"name": "Bernardyn", "url": f"{_RTD}/projects/bernardyn/"},
    {"name": "AIDA", "url": f"{_RTD}/projects/aida/"},
    {"name": "Igor: Irena/Nika/Indra", "url": _IGOR},
]
# ---------------------------------------------------------------------------

# Add a family member only once its RTD project is live (the build uses -W).
intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
    "igor": (_IGOR, None),
    # "family": (f"{_RTD}/en/latest/", None),  # enable once jilavsky-docs is live
}
intersphinx_disabled_reftypes = ["std:doc"]  # cross-project pages need an explicit prefix

html_theme = "pydata_sphinx_theme"
html_title = f"{project} {release}"
html_static_path = ["_static"]
html_css_files = ["family.css"]
html_baseurl = os.environ.get("READTHEDOCS_CANONICAL_URL", "")
html_theme_options = {
    "external_links": [link for link in FAMILY if link["name"] != project],
    "header_links_before_dropdown": 6,
    "navbar_align": "left",
    "show_toc_level": 2,
    "use_edit_page_button": True,
    "icon_links": [
        {
            "name": "GitHub",
            "url": f"https://github.com/jilavsky/{GITHUB_REPO}",
            "icon": "fa-brands fa-github",
        },
    ],
}
html_context = {
    "github_user": "jilavsky",
    "github_repo": GITHUB_REPO,
    "github_version": "main",
    "doc_path": "docs",
}
