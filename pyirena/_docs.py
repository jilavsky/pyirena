"""Links from the application to its online documentation.

Copied from jilavsky-docs/templates/package-docs/pkg/_docs.py; place it at
``<package>/_docs.py``. No Qt imports here: the GUI calls
``QDesktopServices.openUrl(QUrl(docs_url("unified_fit")))`` itself.

The page names in HELP_TARGETS are a contract with released copies of the
program. Never rename a page they point to; add a redirect instead.
tests/test_help_targets.py checks that every target exists.
"""

from importlib.metadata import PackageNotFoundError, version

DIST = "pyirena"  # PyPI distribution name
SUBPROJECT = "pyirena"  # RTD subproject alias under jilavsky-docs
_BASE = f"https://jilavsky-docs.readthedocs.io/projects/{SUBPROJECT}/en/{{ver}}/"

# GUI id -> page path relative to docs/, without suffix ("" = home page).
HELP_TARGETS: dict[str, str] = {
    "home": "",
    "unified_fit": "tools/unified-fit",
    "size_distribution": "tools/size-distribution",
    "modeling": "tools/modeling",
    "simple_fits": "tools/simple-fits",
    "fractals": "tools/fractals",
    "waxs_peakfit": "tools/waxs-peak-fit",
    "carbon_fit": "tools/carbon-fit",
    "saxs_morph": "tools/saxs-morph",
    "data_merge": "tools/data-merge",
    "data_manipulation": "tools/data-manipulation",
    "scattering_contrast": "tools/scattering-contrast",
    "data_explorer": "tools/data-explorer",
    "data_selector": "tools/data-selector",
}

# Explicit MyST labels, "(label)=", that the code links to with an anchor.
# GUI id -> (page path, label)
HELP_ANCHORS: dict[str, tuple[str, str]] = {
    # "unified_fit_levels": ("tools/unified-fit", "unified-fit-levels"),
}


def docs_version() -> str:
    """Return the RTD version to open: "stable" for releases, "latest" otherwise.

    Development, local and pre-release installs ("dev", "+", "a"/"b"/"rc")
    open "latest", which is built from the main branch.
    """
    try:
        v = version(DIST)
    except PackageNotFoundError:
        return "latest"
    is_release = v.replace(".", "").isdigit()
    return "stable" if is_release else "latest"


def docs_url(target: str = "home", anchor: str = "") -> str:
    """Return the documentation URL for a GUI component.

    Args:
        target: Key in HELP_TARGETS or HELP_ANCHORS.
        anchor: Optional explicit label on the page (ignored for HELP_ANCHORS keys).

    Returns:
        Absolute URL on Read the Docs.
    """
    if target in HELP_ANCHORS:
        page, anchor = HELP_ANCHORS[target]
    else:
        page = HELP_TARGETS[target]
    url = _BASE.format(ver=docs_version())
    if page:
        url += page + ".html"
    return url + (f"#{anchor}" if anchor else "")
