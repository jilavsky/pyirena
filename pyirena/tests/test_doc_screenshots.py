"""The documentation screenshots, the shot list and the pages agree.

Screenshots in ``docs/tools/`` have an ID ``<page>-<nn>-<what>`` and live at
``docs/tools/img/<page>/<id>.png``. A page either shows the image (a MyST
``figure``) or holds an invisible placeholder comment
``% SCREENSHOT <id>: <what to shoot>``. Every ID has one row in
``planning/docs-screenshots/SHOTLIST.md``; automated shots have a
``@shot("<id>")`` function in ``scripts/doc_shots/``.

No Qt, no network: this reads text files only. Sphinx ``-W`` already fails
on a figure whose image is missing.
"""

from __future__ import annotations

import re
from collections import Counter
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
DOCS = REPO / "docs"
IMG = DOCS / "tools" / "img"
SHOTLIST = REPO / "planning" / "docs-screenshots" / "SHOTLIST.md"
SHOTS_DIR = REPO / "scripts" / "doc_shots"

ID = r"[a-z0-9]+(?:-[a-z0-9]+)*-\d\d-[a-z0-9]+(?:-[a-z0-9]+)*"
PLACEHOLDER_RE = re.compile(rf"^\s*%\s*SCREENSHOT\s+({ID})\s*:", re.M)
IMAGE_REF_RE = re.compile(rf"img/([a-z0-9-]+)/({ID})\.png")
ROW_RE = re.compile(rf"^\|\s*({ID})\s*\|(.*)\|\s*$", re.M)
SHOT_DECO_RE = re.compile(rf"@shot\(\s*[\"']({ID})[\"']\s*\)")

pytestmark = pytest.mark.skipif(not DOCS.is_dir(), reason="docs/ not present (installed package)")


def _pages() -> dict[Path, str]:
    return {p: p.read_text(encoding="utf-8") for p in DOCS.rglob("*.md")
            if "_build" not in p.parts}


def _placeholders() -> list[tuple[str, Path]]:
    return [(m.group(1), p) for p, text in _pages().items()
            for m in PLACEHOLDER_RE.finditer(text)]


def _image_refs() -> list[tuple[str, str, Path]]:
    """(page folder, id, referencing page) for every docs/tools/img reference."""
    return [(m.group(1), m.group(2), p) for p, text in _pages().items()
            for m in IMAGE_REF_RE.finditer(text)]


def _shotlist() -> dict[str, dict[str, str]]:
    rows = {}
    cols = ["page", "dataset", "state", "must_show", "automated", "status"]
    text = SHOTLIST.read_text(encoding="utf-8")
    for m in ROW_RE.finditer(text):
        cells = [c.strip() for c in m.group(2).split("|")]
        assert len(cells) == len(cols), f"SHOTLIST row {m.group(1)}: expected 7 columns"
        assert m.group(1) not in rows, f"duplicate ID in SHOTLIST: {m.group(1)}"
        rows[m.group(1)] = dict(zip(cols, cells, strict=True))
    return rows


def _page_of(shot_id: str) -> str:
    return re.split(r"-\d\d-", shot_id, maxsplit=1)[0]


def test_shotlist_ids_are_unique_and_well_formed():
    text = SHOTLIST.read_text(encoding="utf-8")
    ids = [m.group(1) for m in ROW_RE.finditer(text)]
    dupes = [k for k, n in Counter(ids).items() if n > 1]
    assert not dupes, f"duplicate IDs in SHOTLIST.md: {dupes}"
    for shot_id, row in _shotlist().items():
        assert row["page"] == _page_of(shot_id), f"{shot_id}: page column is {row['page']!r}"
        assert row["automated"] in ("yes", "no"), f"{shot_id}: automated? must be yes/no"
        assert row["status"] in ("done", "todo"), f"{shot_id}: status must be done/todo"


def test_every_doc_screenshot_is_in_the_shotlist_and_vice_versa():
    in_docs = {i for i, _ in _placeholders()} | {i for _, i, _ in _image_refs()}
    listed = set(_shotlist())
    assert not in_docs - listed, f"in docs/ but not in SHOTLIST.md: {sorted(in_docs - listed)}"
    assert not listed - in_docs, f"in SHOTLIST.md but not in docs/: {sorted(listed - in_docs)}"


def test_placeholder_ids_are_unique_and_not_also_shown():
    ids = [i for i, _ in _placeholders()]
    dupes = [k for k, n in Counter(ids).items() if n > 1]
    assert not dupes, f"placeholder used twice: {dupes}"
    shown = {i for _, i, _ in _image_refs()}
    both = set(ids) & shown
    assert not both, f"image exists but placeholder comment left in: {sorted(both)}"


def test_image_paths_follow_the_convention():
    for folder, shot_id, page in _image_refs():
        assert folder == _page_of(shot_id), f"{page.name}: img/{folder}/{shot_id}.png"
        assert page.parent == DOCS / "tools", f"{shot_id} referenced outside docs/tools/"


def test_every_image_is_referenced_by_a_page():
    referenced = {(f, i) for f, i, _ in _image_refs()}
    on_disk = {(p.parent.name, p.stem) for p in IMG.rglob("*.png")} if IMG.is_dir() else set()
    orphans = on_disk - referenced
    assert not orphans, f"images no page references: {sorted(orphans)}"


def test_status_matches_the_pages():
    shown = {i for _, i, _ in _image_refs()}
    for shot_id, row in _shotlist().items():
        if row["status"] == "done":
            assert shot_id in shown, f"{shot_id}: status done but no page shows the image"
        else:
            assert shot_id not in shown, f"{shot_id}: image is shown; set status to done"


def test_automated_shots_have_a_script_function():
    registered = set()
    for py in SHOTS_DIR.glob("*.py"):
        registered |= set(SHOT_DECO_RE.findall(py.read_text(encoding="utf-8")))
    automated = {k for k, r in _shotlist().items() if r["automated"] == "yes"}
    assert not automated - registered, (
        f"automated: yes but no @shot in scripts/doc_shots/: {sorted(automated - registered)}")
    assert not registered - automated, (
        f"@shot in scripts/doc_shots/ but not automated: yes in SHOTLIST.md: "
        f"{sorted(registered - automated)}")
