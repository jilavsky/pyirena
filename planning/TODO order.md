## Recommended sequence

**2. Step 0 of the config-dialects plan on a short branch off main → `v1.2.0b2`.** I read the code the fixes touch and Step 0 is honestly S, not S–M: the Unified fix is reading `fit_X`/`X_limits` in the bare-number branch of `flatten_level_config` (four lines); the Modeling fix is merging the `uf`/`peak`/`gp`/`mf`/`sf2` sub-block into `d` before `cls.from_dict(d)` in `population_from_dict` (three lines); the Sizes aspect ratio is reordering two branches; the slit fields are two dict entries. The tests are the real work. This can land within days of b1.

**3. Steps 1–4 on `feature/config-envelope` off main, targeting `1.2.0` final** — with one deviation from the plan I'll explain below.

**Branch rule going forward:** every branch comes off `main`, never off another feature branch. The reason zmq contains Carbon is that you branched from a branch. That worked here only because you never touched Carbon on main. Don't branch from zmq again — there'll be nothing left on it after the fast-forward anyway.

## Is the dialect plan good enough?

As diagnosis, it's excellent — I reproduced all four bugs exactly as written, plus the `{tool, model}` envelope silently building an all-defaults model. As a *structural* answer to "make sure nobody can get this out of order again", it's a list of rules, and rules don't enforce themselves. Five things would turn it into a data structure:

**Name the one writer.** The plan says "three readers, one internal representation, one writer going forward" but never names the writer. It already exists: `ToolSetup` in `core/tool_config.py` is exactly *model + fit Q range + held parameters + slit settings*. Give it `to_dict()`/`from_dict()` and make it the *only* thing that writes a config envelope — *Export Parameters*, the HDF5 `_pyirena_config` attribute and `export_results` all call it. Then "is this config complete?" has one answer and one place.

**Pick one parameter encoding and one legacy reader.** §1.5 lists five encodings of "value + fit flag + bounds" and the plan doesn't choose. Choose the one three tools already use — `X`, `fit_X`, `X_limits` — as canonical, and add one small shared function (say `core/param_spec.py: read_param(d, name, defaults) -> (value, fit, lo, hi)`) that accepts all five legacy shapes. Every `from_dict`/legacy reader calls it. The Unified bug is what happens when each reader hand-rolls shape detection; a single function can only be half-written once.

**Version the envelope explicitly instead of sniffing shape.** The panel state already carries `schema_version` (it's in `SizeDis.json`). Put `schema_version` in the `_pyirena_config` header alongside `tool`, `pyirena_version` and the `_note` the plan proposes, bump it when a dialect changes, and let readers branch on it. Shape-sniffing ("a dict means panel dialect") stays only for files that predate the version field. Sniffing is the mechanism that produced the bug in §2.1.

**Add a writer-contract test, not just a round-trip test.** Step 2.3 (invert every flag, narrow every bound, round-trip) is the single most valuable test in the whole plan — write it *before* Step 0 so it goes red then green; that's how you know Step 0 is complete. But add a sibling in the style of `test_tool_registration.py`: enumerate every config writer in the package and assert `build_setup` accepts its output at flags-and-bounds fidelity. That's what permanently closes the "three writers, one reader" class.

**Make Step 3 optional for 1.2.0.** Steps 0–2 make every envelope *readable and correct*; Step 3 (switching the Unified panel and Modeling `_collect_state()` to write the core dialect) only removes writer duplication — and it's the part that "may break a lot", since it changes what goes into every HDF5 file and StateManager, in the 4,400-line file users exercise most. The plan doesn't say this, but it's a real option: ship 1.2.0 with Steps 0–2, do Step 3 for 1.3 on its own branch with its own release note. Also decide Fractals now (leave it on the legacy reader, documented) rather than leaving the "or".

## Small things worth doing before 1.2.0 final

From the unfinished list, item 6 (shipped "Phase 1" strings in MCP/ZMQ replies) is cheap and user-visible; item 4 (Sizes `unc_n_runs` ignored on replay) fits naturally into Step 0's tests; item 2 (MC uncertainties missing from five `run_*_fit` control functions) is the gap an agent will notice first but is not a blocker. The SAXS Morph decision (item 3) is worth ten minutes to *decide*, whichever way. And per its own opening paragraph, `planning/zmq-service/README.md` should shrink to §11's leftovers after the merge — a 750-line shipped plan is what you said you didn't want.

One documentation point I checked because it matters for USAXS: over ZMQ, slit smearing for Sizes comes from `data.is_slit_smeared`/`data.slit_length`, not from the config (the GUI export doesn't carry it). That is documented in `docs/zmq_service.md`, but the orchestrator team should be told explicitly, because a forgotten flag gives a plausible pinhole fit rather than an error.

Want me to do the merge prep now — version bump, carbon README status fix, changelog check — so you can review it as a single commit before tagging?