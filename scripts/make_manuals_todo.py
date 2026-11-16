#!/usr/bin/env python
"""Write planning/manuals-todo.md: everything in the manuals that needs Jan's eyes.

Collects every ``% CHECK(Jan): …`` comment in ``docs/`` (a MyST comment, so it
never appears on the published site) with its file and line, and every
``automated: no`` row of ``planning/docs-screenshots/SHOTLIST.md``, into one
checklist for the documentation walk-through. Re-run after resolving items:

    python scripts/make_manuals_todo.py

The comments stay the source of truth: resolve an item by fixing the text and
deleting its comment, then re-run this script.
"""

from __future__ import annotations

import datetime as _dt
import re
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
DOCS = REPO / "docs"
SHOTLIST = REPO / "planning" / "docs-screenshots" / "SHOTLIST.md"
OUT = REPO / "planning" / "manuals-todo.md"

#: Tool pages first, in Help-button order; everything else after, by path.
PAGE_ORDER = [
    "unified-fit", "size-distribution", "modeling", "simple-fits", "fractals",
    "waxs-peak-fit", "carbon-fit", "saxs-morph", "data-merge", "data-manipulation",
    "scattering-contrast", "data-explorer", "data-selector", "feature-identifier",
]

CHECK_RE = re.compile(r"^\s*%\s*CHECK\(Jan\):\s*(.*)$")
CONT_RE = re.compile(r"^\s*%\s+(?!CHECK|SCREENSHOT)(.*)$")
ROW_RE = re.compile(r"^\|\s*([a-z0-9-]+-\d\d-[a-z0-9-]+)\s*\|(.*)\|\s*$", re.M)


def collect_checks() -> dict[Path, list[tuple[int, str]]]:
    """``{page: [(line, text), …]}`` for every CHECK comment (continuation lines joined)."""
    found: dict[Path, list[tuple[int, str]]] = {}
    for page in sorted(DOCS.rglob("*.md")):
        if "_build" in page.parts:
            continue
        lines = page.read_text(encoding="utf-8").splitlines()
        i = 0
        while i < len(lines):
            m = CHECK_RE.match(lines[i])
            if not m:
                i += 1
                continue
            line_no, text = i + 1, [m.group(1).strip()]
            i += 1
            while i < len(lines) and (c := CONT_RE.match(lines[i])):
                text.append(c.group(1).strip())
                i += 1
            found.setdefault(page, []).append((line_no, " ".join(text)))
    return found


def collect_manual_shots() -> list[list[str]]:
    """``[id, dataset, state, must_show]`` for every ``automated: no`` shot."""
    rows = []
    for m in ROW_RE.finditer(SHOTLIST.read_text(encoding="utf-8")):
        cells = [c.strip() for c in m.group(2).split("|")]
        if len(cells) == 6 and cells[4] == "no":
            rows.append([m.group(1), cells[1], cells[2], cells[3]])
    return rows


def _sort_key(page: Path):
    stem = page.stem
    if page.parent == DOCS / "tools" and stem in PAGE_ORDER:
        return (0, PAGE_ORDER.index(stem), "")
    return (1, 0, str(page))


def main() -> int:
    checks = collect_checks()
    shots = collect_manual_shots()
    n_checks = sum(len(v) for v in checks.values())
    out = [
        "# Manuals to do",
        "",
        "Everything in the manuals that needs Jan's eyes, in one list for the walk-through.",
        f"Generated {_dt.date.today():%d-%m-%Y} by `python scripts/make_manuals_todo.py` from the",
        "`% CHECK(Jan): …` comments in `docs/` (invisible on the site) and the `automated: no` rows",
        "of `planning/docs-screenshots/SHOTLIST.md`. Resolve an item by fixing the page and deleting",
        "its comment (or taking the screenshot), then re-run the script.",
        "",
        f"**{n_checks} text checks, {len(shots)} screenshots.**",
        "",
        "## Screenshots to take",
        "",
        "Save as `docs/tools/img/<page>/<id>.png`, replace the `% SCREENSHOT <id>` comment in the",
        "page with a `figure` (alt text and caption), and set the row in `SHOTLIST.md` to `done`.",
        "",
        "| ID | dataset | state to set up | must show |",
        "|---|---|---|---|",
    ]
    out += [f"| {' | '.join(r)} |" for r in shots]
    out += ["", "## Text to check, by page", ""]
    for page in sorted(checks, key=_sort_key):
        rel = page.relative_to(DOCS).as_posix()
        out += [f"### {rel}", ""]
        out += [f"- [ ] [line {n}](../docs/{rel}#L{n}) — {text}" for n, text in checks[page]]
        out.append("")
    OUT.write_text("\n".join(out).rstrip() + "\n", encoding="utf-8")
    print(f"{OUT.relative_to(REPO)}: {n_checks} text checks, {len(shots)} screenshots")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
