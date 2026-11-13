# pyIrena — Code Quality Plan

Living record of what is still open and, more importantly, **what was
deliberately decided against**. Completed work is not repeated here — it is in
[CHANGELOG.md](CHANGELOG.md) (current series) and
[docs/changelog-archive.md](docs/changelog-archive.md) (1.0.1 and earlier),
with the reasoning for cross-cutting choices in
[docs/developer/adding-features.md](docs/developer/adding-features.md).

A decision to *not* do something is the part that gets re-litigated every few
months, so it is the part worth keeping.

---

## Open — no bugs, no user impact

**Large-file comfort refactors.** Revisit if these keep growing or become hard
to navigate. **The trigger is met** (review of 09-11-2026, issue #38): all four
grew since this list was written. Sequence it after the shared fit worker and
a `FitPanelBase` (review §4.1/§4.2), which give the split its natural seam.

- `gui/data_selector/panel.py` (~2,800 → 3,117 lines) — split by mixin/topic
- `gui/unified_fit.py` (~4,300 → 4,420 lines) — split along panel/widget lines
- `gui/modeling_panel.py` (~3,800 → 3,962 lines)
- `gui/sizes_panel.py` (~3,000 → 3,230 lines)

**Ruff `UP` and `B` passes, and `ruff format`.** House standard, not yet
enabled here. Measured 09-11-2026: `UP` ≈ 1,500 auto-fixes, `B` ≈ 120 manual
(look first at `B905` bare `zip`, `B904` in `core/scattering_contrast.py`,
`B023` in `core/reporting.py`), and `ruff format` would rewrite 251 of 271
files. Deferred deliberately: each is a whole-tree diff that conflicts with
every open branch, so it runs on its own branch once those are merged — `UP`
and `B` as two separate commits. Whether to adopt `ruff format` at all is
still to be decided when that branch is cut; until then, do not reformat
files you are not otherwise changing.

**Shared GUI fit orchestration (`FitPanelBase`).** The batch/API half of the
fit-runner work shipped with #38; the panels still orchestrate their own
fits. Plan and the one GUI bug waiting on it are in
`planning/fit-runner/README.md` ("Remaining").

**One NXcanSAS reader, not two generations** (review §3.8). `io/hdf5.py`
(camelCase `readGenericNXcanSAS` / `readMyNXcanSAS` / `saveNXcanSAS`, no type
hints) sits beside the documented `io/scattering.py` stable API. Make the old
functions thin shims over the new ones so there is one reader to fix. Saved
files must not change.

**Very long functions that hold logic** (review §3.7), split when next
touched: `core/reporting.py:_build_report` (~690 lines, one branch per tool —
a per-tool builder table would match the rest of the design);
`gui/data_selector/panel.py:tabulate_results` (~480 lines — check whether it
duplicates `api/aggregate.py`; if so the panel should call that);
`io/hdf5.py:saveNXcanSAS` / `readMyNXcanSAS` (with the item above). Pure
widget-assembly giants (`sizes_panel._create_control_panel`,
`data_selector.init_ui`) only along group boxes, when touched.

**Tests: a coverage floor and a panel smoke test** (review §3.10). CI runs
`pytest --cov` and discards the number; publish it and set a floor for
`core/` and `io/` only. Add an offscreen test that constructs every
registered panel (the contract tests check source conventions, not
construction).

**Small packaging leftovers.** `pyIrena.icns` (470 KB) and `.kilo/` are
tracked in git; `pyirena/examples/` ships in the wheel. Decide per item.

**`print()` in library code** — mostly resolved (review §3.5). What is left
is deliberate console output: the `print_*_results` helpers in
`io/nxcansas_*.py` and CLI `main()`s. The MCP server logs to stderr; convert
any remaining print only if it becomes reachable from an MCP/ZMQ tool, where
stdout is the protocol stream.

**Batch path for Data Manipulation.** Wanted, but the tool needs revising
first. Filed as its own GitHub issue (#19).

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
`docs/developer/adding-features.md`.

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
