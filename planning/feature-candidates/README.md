# Feature candidates — what users might want next

Internal planning artifact. Collected from the 2026-11 review (§5) on
11-11-2026; ranked by estimated user value, to be re-ranked by what Irena
users actually ask for. Nothing here is started. When one is picked up, give
it its own folder (design + decisions) and remove it from this list.

| # | Feature | Issue | Notes |
|---|---|---|---|
| 1 | **Q-resolution (dQ) smearing in fits** | #30 | Data, sessions and io already carry `dQ`; only slit smearing is applied. For SANS and pinhole users this decides whether a fit is quantitative. `core/smearing.py` already has the matrix machinery. Decide first — see `planning/unfinished/README.md` §1. |
| 2 | **More form and structure factors in Modeling** | — | pyIrena has 10 form factors and 2 structure factors; Igor Irena offers more (tube, core-shell cylinder, fractal aggregate, disk/rod variants; sticky hard sphere, square well, Hayter-Penfold MSA, interprecipitate). `docs/developer_adding_*` make each a contained task; pick by user request. |
| 3 | **Series / in-situ workflow from the GUI** | — | Fit a sequence using the previous result as the start, live trend plot of chosen parameters, stop on a bad fit. The pieces exist (batch, tabulate, Data Explorer trend plots); check how much is already wired before designing. |
| 4 | **Data Manipulation: trim/delete points** | #36 | Requested. |
| 5 | **Data Manipulation: batch path** | #19 | Requested; the tool needs revising first (`PLAN.md`). |
| 6 | **Pair-distance distribution p(r)** | — | No PDDF tool yet; the regularisation code in Size Distribution is most of the maths. Valuable for solution-scattering users. |
| 7 | **Async jobs over ZMQ** | #33 | Only if Modeling timeouts actually appear — `planning/zmq-service/README.md` §3.2. |
| 8 | **Simultaneous multi-dataset fit** (contrast variation, SAXS+SANS) | — | Larger; only on demand. |
| 9 | **Built, searchable documentation site** (mkdocs or Sphinx on GitHub Pages) | — | The 38 Markdown files are browsed on GitHub today. Needs a `docs` extra — a dependency decision. |
| 10 | **Recent files / folders** in the file browsers | — | No such memory exists today (`gui/window_state.py`, state manager). |

**Not recommended now:** reflectivity and desmearing. Model-side smearing,
which pyIrena already does, is the better approach, and neither has been
requested.
