"""The documentation screenshots, the shot list and the pages agree.

Copied from jilavsky-docs/templates/package-docs (derived from pyIrena).
A screenshot has an ID ``<page>-<nn>-<what>``, where ``<page>`` is the file
name of a page under ``docs/``, and lives next to it at
``<page dir>/img/<page>/<id>.png``. A page either shows the image (a MyST
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


def _find_repo() -> Path:
    for parent in Path(__file__).resolve().parents:
        if (parent / "pyproject.toml").is_file():
            return parent
    return Path(__file__).resolve().parent


REPO = _find_repo()
DOCS = REPO / "docs"
SHOTLIST = REPO / "planning" / "docs-screenshots" / "SHOTLIST.md"
SHOTS_DIR = REPO / "scripts" / "doc_shots"

ID = r"[a-z0-9]+(?:-[a-z0-9]+)*-\d\d-[a-z0-9]+(?:-[a-z0-9]+)*"
PLACEHOLDER_RE = re.compile(rf"^\s*%\s*SCREENSHOT\s+({ID})\s*:", re.M)
IMAGE_REF_RE = re.compile(rf"img/([a-z0-9-]+)/({ID})\.png")
ROW_RE = re.compile(rf"^\|\s*({ID})\s*\|(.*)\|\s*$", re.M)
SHOT_DECO_RE = re.compile(rf"@shot\(\s*[\"']({ID})[\"']\s*\)")
COLS = ["page", "dataset", "state", "must_show", "automated", "status"]

pytestmark = pytest.mark.skipif(
    not (DOCS.is_dir() and SHOTLIST.is_file()), reason="docs/ or SHOTLIST.md not present"
)


def _pages() -> dict[Path, str]:
    return {
        p: p.read_text(encoding="utf-8")
        for p in DOCS.rglob("*.md")
        if "_build" not in p.parts and not p.name.startswith("_")
    }


def _placeholders() -> list[tuple[str, Path]]:
    return [(m.group(1), p) for p, text in _pages().items() for m in PLACEHOLDER_RE.finditer(text)]


def _image_refs() -> list[tuple[str, str, Path]]:
    """(img folder, id, referencing page) for every img/<page>/<id>.png reference."""
    return [
        (m.group(1), m.group(2), p)
        for p, text in _pages().items()
        for m in IMAGE_REF_RE.finditer(text)
    ]


def _shotlist() -> dict[str, dict[str, str]]:
    rows = {}
    for m in ROW_RE.finditer(SHOTLIST.read_text(encoding="utf-8")):
        cells = [c.strip() for c in m.group(2).split("|")]
        assert len(cells) == len(COLS), f"SHOTLIST row {m.group(1)}: expected 7 columns"
        assert m.group(1) not in rows, f"duplicate ID in SHOTLIST: {m.group(1)}"
        rows[m.group(1)] = dict(zip(COLS, cells, strict=True))
    return rows


def _page_of(shot_id: str) -> str:
    return re.split(r"-\d\d-", shot_id, maxsplit=1)[0]


def test_shotlist_rows_are_well_formed():
    ids = [m.group(1) for m in ROW_RE.finditer(SHOTLIST.read_text(encoding="utf-8"))]
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
    both = set(ids) & {i for _, i, _ in _image_refs()}
    assert not both, f"image exists but placeholder comment left in: {sorted(both)}"


def test_images_sit_next_to_their_page():
    for folder, shot_id, page in _image_refs():
        assert folder == _page_of(shot_id), f"{page.name}: img/{folder}/{shot_id}.png"
        assert folder == page.stem, f"{shot_id} is shown on {page.name}, not on {folder}.md"
        img = page.parent / "img" / folder / f"{shot_id}.png"
        assert img.is_file(), f"{page.name}: missing {img.relative_to(REPO)}"


def test_every_image_is_referenced_by_a_page():
    referenced = {(page.parent / "img" / f / f"{i}.png").resolve() for f, i, page in _image_refs()}
    on_disk = {
        p.resolve()
        for p in DOCS.rglob("img/*/*.png")
        if "_build" not in p.parts and re.fullmatch(ID, p.stem)
    }
    orphans = sorted(str(p.relative_to(REPO)) for p in on_disk - referenced)
    assert not orphans, f"images no page references: {orphans}"


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
        if not py.name.startswith("_"):
            registered |= set(SHOT_DECO_RE.findall(py.read_text(encoding="utf-8")))
    automated = {k for k, r in _shotlist().items() if r["automated"] == "yes"}
    assert not automated - registered, (
        f"automated: yes but no @shot in scripts/doc_shots/: {sorted(automated - registered)}"
    )
    assert not registered - automated, (
        f"@shot in scripts/doc_shots/ but not automated: yes in SHOTLIST.md: "
        f"{sorted(registered - automated)}"
    )
