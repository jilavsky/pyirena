"""Every Help-button target in the code must exist in docs/.

Copied from jilavsky-docs/templates/package-docs. Catches "renamed the page,
broke the Help button" at PR time. Needs no Sphinx and no network.
"""

import re
from pathlib import Path

import pytest

from pyirena._docs import HELP_ANCHORS, HELP_TARGETS

DOCS = Path(__file__).resolve().parents[2] / "docs"


def _page_file(page: str) -> Path:
    return DOCS / ("index.md" if page == "" else f"{page}.md")


@pytest.mark.parametrize("target,page", sorted(HELP_TARGETS.items()))
def test_help_target_page_exists(target, page):
    assert _page_file(page).is_file(), f"HELP_TARGETS[{target!r}] -> missing docs/{page}.md"


@pytest.mark.parametrize("target,spec", sorted(HELP_ANCHORS.items()))
def test_help_anchor_label_exists(target, spec):
    page, label = spec
    path = _page_file(page)
    assert path.is_file(), f"HELP_ANCHORS[{target!r}] -> missing docs/{page}.md"
    labels = re.findall(r"^\(([^)]+)\)=\s*$", path.read_text(encoding="utf-8"), flags=re.M)
    assert label in labels, f"HELP_ANCHORS[{target!r}] -> no '({label})=' in docs/{page}.md"
