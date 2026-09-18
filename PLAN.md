# pyIrena — Code Quality Plan

Living record of what is still open and, more importantly, **what was
deliberately decided against**. Completed work is not repeated here — it is in
[CHANGELOG.md](CHANGELOG.md) (current series) and
[docs/CHANGELOG_archive.md](docs/CHANGELOG_archive.md) (1.0.1 and earlier),
with the reasoning for cross-cutting choices in
[docs/developer_adding_features.md](docs/developer_adding_features.md).

A decision to *not* do something is the part that gets re-litigated every few
months, so it is the part worth keeping.

---

## Open — no bugs, no user impact

**Large-file comfort refactors.** Revisit if these keep growing or become hard
to navigate; none of them is causing a problem today.

- `gui/data_selector/panel.py` (~2,800 lines) — split by mixin/topic
- `gui/unified_fit.py` (~4,300 lines) — split along panel/widget lines
- `gui/modeling_panel.py` (~3,800 lines)
- `gui/sizes_panel.py` (~3,000 lines)

**Batch path for Data Manipulation.** Wanted, but the tool needs revising
first. Filed as its own GitHub issue.

**MCP 2.x migration.** `mcp` is pinned `<2`. Migrating is its own project;
nothing needs it yet, and it should be decided before anything new is built on
top of the dispatcher.

---

## Decided against — deliberate, not forgotten

### A shared `FileBrowserWidget`

The four browsers now share all the *logic* — filtering, sorting, the
file-type table, folder listing, drag-and-drop. What is still duplicated is
widget assembly with no logic in it, and the four differ in ways a common
widget would have to be configured around anyway: the Data Explorer is a tree
with lazy sub-folder expansion, Data Merge shows two linked instances, Data
Manipulation adds a context menu, and the Data Selector alone offers text files
and the convert-on-load path. **Duplicated assembly with no logic in it does
not drift; duplicated logic does.** Recorded in
`docs/developer_adding_features.md`.

### A shared batch config loader

Data Merge and Data Manipulation each run from their **own** JSON at a specific
point in an instrument reduction pipeline — they produce new *data files* that
a later stage consumes, so their config lives in a different folder and a
different pipeline step from the analysis config. Unifying the loader would
break those pipelines.

### A batch path for Fractals

Fractals is a visualization tool, not an analysis technique. It has no batch
use case and no JSON config by design.

### Exhaustive schema default-value parity tests

The parity test already enforces `properties ⊆ parameters`,
`required == mandatory parameters`, and that every parameter is exposed.
Asserting each JSON `default` equals the Python default is brittle — several
defaults are descriptive (e.g. "scipy default") — for little extra safety.

### A dedicated lowest-supported-dependency CI job

The NumPy floor matches the code (≥ 2.0), so there is no 1.x/2.x split left to
guard. A lowest-pins resolver job would add CI cost for no benefit unless lower
floors return.

### A full install smoke matrix (core / plotting / gui / mcp / all)

CI already runs a core `test` job and a `[gui]` `test-gui` job. A 5-way
extras × OS install matrix is more machinery than a single-maintainer
scientific package needs.

### A heavier release pipeline

The cheap, high-value guard is in place: `publish.yml` verifies the release tag
equals the `pyproject.toml` version. Gating publish on a re-run of the Tests
workflow, and reusing exact previously-validated artifacts, are deferred.

---

## How this file is used

Add an entry when a review or a design discussion produces a decision that will
not be obvious from the code six months later. Remove an "Open" entry when it
ships — the CHANGELOG records what happened; this file records what is left and
what we chose not to do.
