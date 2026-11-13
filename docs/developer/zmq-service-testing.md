# ZMQ service — deployment and test checklist

A checklist to work through once, on a test machine, before the service goes
anywhere near a beamtime. Tick the boxes as you go; anything that fails has a
line in [Troubleshooting](#troubleshooting).

Reference documentation is [zmq_service.md](zmq_service.md) — this file is the
procedure, that one is the protocol.

**You will need:** the analysis workstation (RHEL — called *server* below),
the machine the orchestrator runs on (*client*), and about an hour. Parts 1–3
are on the server, Part 4 onwards needs both.

---

## Part 1 — Install on the server

- [X] **1.1** Pick the Python. Anything ≥ 3.10; 3.12 or 3.13 preferred.
      ```bash
      python3 --version
      ```

- [X] **1.2** Make an environment of its own, so the service cannot be broken
      by an unrelated `pip install`.
      ```bash
      python3 -m venv ~/pyirena-zmq-env
      source ~/pyirena-zmq-env/bin/activate
      ```
      (or `conda create -n pyirena-zmq python=3.13 && conda activate pyirena-zmq`)

- [X] **1.3** Install pyIrena with the `zmq` extra. From a checkout of this
      branch:
      ```bash
      cd /path/to/pyirena
      pip install -e ".[zmq]"
      ```
      No Qt and no matplotlib are needed — the service runs headless.

- [X] **1.4** Check the entry point exists.
      ```bash
      pyirena-zmq --help
      ```
      → the usage block, listing `--port`, `--request-budget-s`, and the rest.

- [X] **1.5** Note the versions for your record:
      ```bash
      python -c "import pyirena, zmq, sys; print('pyirena', pyirena.__version__, '| pyzmq', zmq.__version__, '| python', sys.version.split()[0])"
      ```
      pyirena 1.1.1  pyzmq 27.2.0  python 3.13.15

---

## Part 2 — First run, on loopback

Prove it works before exposing it to the network.

- [X] **2.1** Start it bound to loopback only, logging to the terminal:
      ```bash
      pyirena-zmq --bind tcp://127.0.0.1:9865 --log-stderr
      ```
      → `pyirena-zmq listening on tcp://127.0.0.1:9865 (profile=json_only,
      budget=55s, max_message=16MB)`

- [X] **2.2** In a **second terminal on the server**, run the smoke test:
      ```bash
      source ~/pyirena-zmq-env/bin/activate
      python /path/to/pyirena/scripts/zmq_smoke_test.py tcp://127.0.0.1:9865
      ```
      → `19 passed, 0 failed` (17 without `--config`).
      NOTE: reports only 16 passed, 0 failed. Some seemed to be skipped... 

      The fit in step 4 is a plumbing check, not a scientific one — it runs a
      default one-level Unified Fit on a synthetic curve, so ignore the χ².
      Part 5 is where the numbers matter.

- [X] **2.3** Watch the first terminal. One log line per request:
      `id=… op=call tool=run_fit session=… ok 0.061s`

- [X] **2.4** Stop it with Ctrl+C.
      → `pyirena-zmq stopped after N request(s)`, and the prompt returns
      promptly. If Ctrl+C hangs, say so — that is a bug, not your machine.

- [X] **2.5** Check the log file was written:
      ```bash
      ls -l ~/.pyirena/logs/zmq.log && tail -5 ~/.pyirena/logs/zmq.log
      ```

---

## Part 3 — Expose it to the network

- [X] **3.1** Find the server's address, and the client's:
      ```bash
      hostname -f ; ip -4 addr show | grep inet
      ```
      server usaxscontrol.xray.aps.anl.gov  client ANLUF4C4VKRQ09

- [X] **3.2** Start it on all interfaces (the default):
      ```bash
      pyirena-zmq --port 9865 --log-stderr
      ```
      → a second warning line: *bound to every interface with no
      authentication — restrict this port to the orchestrator host at the
      firewall*. That warning is the reason for the next step.

- [ ] **3.3** Open the port **to the client only**, not to the world:
      ```bash
      sudo firewall-cmd --permanent --add-rich-rule \
        'rule family="ipv4" source address="CLIENT.IP.GOES.HERE" port port="9865" protocol="tcp" accept'
      sudo firewall-cmd --reload
      sudo firewall-cmd --list-rich-rules
      ```

- [ ] **3.4** Confirm nothing else already owns 9865 (you mentioned another
      service on `purple`):
      ```bash
      ss -tlnp | grep 9865
      ```
      → exactly one line, and it is this service.

---

## Part 4 — Reach it from the client

- [X] **4.1** On the client, install pyzmq. **Nothing else is needed** — the
      smoke test deliberately has no pyIrena, numpy or h5py dependency.
      ```bash
      pip install pyzmq
      ```

- [X] **4.2** Copy the test script over:
      ```bash
      scp user@server:/path/to/pyirena/scripts/zmq_smoke_test.py .
      ```

- [X] **4.3** Run it across the network:
      ```bash
      python zmq_smoke_test.py tcp://usaxscontrol:9865
      ```
      → `17 passed, 0 failed`.
      NOTE: 16 passes, 0 failed. 
      **This is the step that has never been tested.** Everything up to here
      has been exercised on loopback in CI; a real network path between two
      machines has not. If something is going to be wrong, it is most likely
      here.

- [X] **4.4** Note the timings from step 4 of the output. Compare with the
      loopback run: the difference is your network overhead for a ~30 KB
      request. Expect milliseconds; if it is seconds, investigate before
      building on it.

      Measured 29-09-2026, laptop → `usaxscontrol`: **0.6 ms** median for a
      tiny request (`server_info` ×10), and **9.7 ms** median of
      *wall minus the server's own `elapsed_s`* for a 165 KB request (×5).
      Transport is not a factor at these sizes.

---

## Part 5 — Your real data and configs

This is where you find out whether the answers are right, not just whether
the plumbing works.

- [ ] **5.1** Put one measured curve into the JSON the service takes. On a
      machine with pyIrena:
      ```python
      import json
      import numpy as np
      from pyirena.io.hdf5 import readGenericNXcanSAS

      d = readGenericNXcanSAS("testData/Scripting", "Modeling_PP15.h5")
      json.dump({"q": np.asarray(d["Q"]).tolist(),
                 "intensity": np.asarray(d["Intensity"]).tolist(),
                 "error": np.asarray(d["Error"]).tolist(),
                 "label": "PP15"}, open("curve.json", "w"))
      ```

- [X] **5.2** Replay a real exported config through `analyze`:
      ```bash
      python zmq_smoke_test.py tcp://usaxscontrol:9865 \
             --config modeling.json --data curve.json
      ```
      → section 5 of the output reports the tool, χ², the fitted Q range and
      any notes.

- [X] **5.3** **Compare the numbers against the GUI.** Open the same file and
      config in the pyIrena GUI, fit, and check χ² and the parameters match.
      This is the check that matters; everything else is transport.

      GUI χ² 2455138  service χ² 2325500 = close enough

      **Worth one more look, though.** Checked 29-09-2026: a Modeling fit
      through `analyze` is *bit-deterministic* — four runs of
      `modeling.json` + `Modeling_PP15.h5` gave
      `297.608163707637` every time, and the server gave the same. So a 5 %
      gap is not run-to-run scatter; something genuinely differed between
      the two fits. The likely candidate is visible in the reply's
      `quality.message` for that config:

      > population 2 B = 7.45e-06 is pinned at its upper fit limit (7.45e-06)

      A parameter sitting on its limit is exactly where two runs with
      slightly different limits, Q range or error scaling land in different
      places. Compare `quality.message`, `analyze.fit_q_min/fit_q_max` and
      the fit limits in the config against what the GUI had, and the 5 %
      should explain itself. (Your numbers are ~2.4 M absolute χ², whereas
      `modeling.json` + `Modeling_PP15.h5` gives 151 780 — so you fitted a
      different config or curve, and I could not reproduce your pair here.)

- [X] **5.4** Repeat 5.2 for each tool you will use, and 5.3 (the GUI
      comparison) for the ones whose numbers you care about. Configs and data
      are in `testData/Scripting/`.

      Run 29-09-2026 against the live service on `usaxscontrol:9865`, each
      config replayed through `analyze` over the wire and the identical call
      also run locally, so the two numbers are directly comparable:

      | Tool | config | data | χ² over the wire | vs local |
      |---|---|---|---|---|
      | Modeling | `modeling.json` | `Modeling_PP15.h5` | 297.6081637076363 | 2e-15 |
      | Modeling | `modeling.json` | `ModelingSF_SaD.h5` | 58284.20959373278 | **identical** |
      | Unified Fit | `UnifiedFit` | `UnifiedFit_PP15.h5` | 23.393567477100092 | 1e-15 |
      | Size Distribution | `SizeDis.json` | `SizeDis_PP15.h5` | 208.8052021373908 | 2e-8 |
      | Simple Fits | `SimpleFits_Porod.json` | `SimpleFits_Porod.h5` | 7.7786373968865545 | 3e-16 |
      | WAXS Peak Fit | `WAXS` | `WAXS_Al_7075.hdf` | 0.022413849290221697 | 5e-15 |
      | Carbon model | `carbon_CE_1400.json` | `Carbon_CE_1400.h5` | 42.951768104669846 | 1e-13 |

      The last column is the worst relative difference in χ²; every tool also
      had its whole parameter block compared (87 numbers for Carbon, 25 for
      Sizes), worst disagreement anywhere **3.8e-5**, in a WAXS peak FWHM.
      These are floating-point last-digit differences between two machines,
      not disagreements — nothing here needs chasing.

      **saxsMorph is not one of the tools `analyze` knows** (`BAD_CONFIG`,
      naming the six it does). The config in `testData/Scripting/saxsMorph.json`
      is therefore untestable this way. It was listed as optional; if it is
      ever needed over the service, that is a feature, not a bug fix.

      The server is roughly **3× slower per fit than a laptop** —
      `ModelingSF_SaD` took 13.4 s there against 4.8 s here. Budget
      accordingly: the 55 s ceiling has less headroom than a local timing
      suggests.

- [X] **5.5** Read the `notes` in the reply. A clipped Q range or a failed
      pre-fit is reported there rather than thrown away, and a note you did
      not expect is worth chasing.

      All notes from the 5.4 run, and all of them expected:
      - `ModelingSF_SaD` — *"The config's Q range was clipped to this curve:
        [0.0001385, 0.3]."* The config was exported against a wider curve.
      - `SizeDis_PP15` — *"Power-law pre-fit: B = 0.0004554, P = 3.183"* and
        *"Background pre-fit: 0.4996 cm⁻¹"*. The background matches the
        `background: 0.499628` in the config, which is the check that the
        pre-fit ran on the right Q window.
      - `WAXS_Al_7075` — *"Per-peak Q0 scan done (window=±0.036 1/A,
        50 steps)."*
      - Modeling/PP15, Unified Fit, Simple Fits, Carbon — no notes.

- [ ] **5.6** Time your **slowest realistic** fit and compare it with the
      55 s budget. Modeling with several populations is the one to watch.
      Monte-Carlo uncertainties are off by default for exactly this reason —
      at `n_mc_runs=50` a Modeling fit takes minutes, not seconds.

      Slowest over the wire in the 5.4 run: **13.4 s** (Modeling,
      `ModelingSF_SaD`), against the 55 s budget. That config carries
      `n_mc_runs=50`; uncertainties were not requested, and if they are, the
      same fit runs for minutes and will hit `TIMEOUT`.

      slowest observed with *your* data ________ s

---

## Part 6 — Make it fail on purpose

Five minutes here saves an argument during beamtime about whose fault it is.

- [X] **6.1** **Client timeout.** Ask for an impossible deadline and confirm
      you get a `TIMEOUT` *reply* rather than silence:
      ```bash
      python - <<'EOF'
      import json, zmq
      s = zmq.Context.instance().socket(zmq.REQ); s.connect("tcp://usaxscontrol:9865")
      s.send_string(json.dumps({"op": "call", "tool": "analyze",
                                "args": {"data": json.load(open("curve.json")),
                                         "config": json.load(open("modeling.json"))},
                                "options": {"request_budget_s": 0.5}}))
      print(json.loads(s.recv_string())["error"]["code"])
      EOF
      ```
      → `TIMEOUT`. The server answering instead of going quiet is the whole
      point of the design — if you ever see silence, report it.

- [X] **6.2** **Kill the client mid-request** (Ctrl+C during a fit), then run
      the smoke test again. The service must still be healthy.

      Done 29-09-2026: a client was closed 1.0 s into a ~13 s Modeling fit.
      The service stayed healthy and the next fit returned a bit-identical
      χ². **But the next request waited 12.3 s** — abandoning a client does
      not stop the fit, it only stops anyone hearing the answer, and the
      service handles one request at a time. So a client that times out and
      retries immediately puts the retry behind the original fit, where it
      can time out in turn. Tell the orchestrator team to back off between
      retries rather than hammer.

- [ ] **6.3** **Restart the server** while the client holds a session id. The
      next call returns `NO_SESSION` — sessions do not survive a restart, by
      design. Make sure the orchestrator team knows.

      The restart itself still needs doing on `usaxscontrol` — but the
      symptom was confirmed 29-09-2026 by presenting an id the server never
      issued: `get_session_summary`, `call run_fit` and `close_session` all
      return `NO_SESSION` with *"No session found with id '…'."* A full
      session round trip over the wire (`open_dataset` with arrays →
      `list_sessions` → `get_session_summary` → `close_session` → gone) also
      works, and a session with no model attached refuses to fit with
      `NO_MODEL` rather than inventing one.

- [X] **6.4** **A file tool is refused.** Already covered by the smoke test
      (`NOT_AVAILABLE_OVER_ZMQ`), but confirm the team understands *why*:
      the service shares no filesystem with them, so a tool that returns a
      path would be lying.

      Confirmed 29-09-2026 on the live service, in both directions: the 15
      file- and image-returning tools are refused by `call` *and* by
      `describe_tool` with the same code, and they are absent from
      `list_tools` — so an agent enumerating the service never learns they
      exist. Message: *"Tool 'save_fit' reads or writes a file on the
      server's filesystem and is not available in profile 'json_only'."*

- [X] **6.5** **Bad input of every shape**, all answered rather than dropped:
      `BAD_JSON` (a bare string), `BAD_ENVELOPE` (no `op`), `UNKNOWN_OP`,
      `UNKNOWN_TOOL`, `SHAPE_MISMATCH` (q and I different lengths),
      `EMPTY_DATA`, `BAD_VALUES` (a string in the intensities),
      `TOO_MANY_POINTS` (200 000 against the 100 000 cap, refused in 0.2 s
      without attempting the fit).

- [X] **6.6** **Per-call options narrow but never widen.** Asking for a 600 s
      budget is refused — *"the server's budget is 55 s; a call may shorten
      its own deadline but not extend it"* — while asking for 5 s is granted.
      `allow_images: true` is refused on a `json_only` server; `max_points:
      50` is granted. An unrecognised key is `BAD_OPTION` rather than being
      silently ignored, which is what you want when a typo would otherwise
      look like it worked.

---

## Part 7 — Run it as a service

Only once Parts 1–6 pass.

- [ ] **7.1** Write a config file so the options are not buried in a unit
      file. `/etc/pyirena/zmq.json`:
      ```json
      {
        "bind": "tcp://0.0.0.0:9865",
        "request_budget_s": 55,
        "max_sessions": 16,
        "log_level": "INFO"
      }
      ```

- [ ] **7.2** Install the unit — the template is in
      [zmq_service.md](zmq_service.md#keeping-it-running-systemd). Set `User=`
      to a service account, and `ExecStart=` to the venv's `pyirena-zmq`.

- [ ] **7.3** Start and enable it:
      ```bash
      sudo systemctl daemon-reload
      sudo systemctl enable --now pyirena-zmq
      systemctl status pyirena-zmq
      ```

- [ ] **7.4** Run the smoke test once more, from the client, against the
      service-managed instance.

- [ ] **7.5** Confirm it comes back after a reboot, or at least after
      `sudo systemctl restart pyirena-zmq`.

- [ ] **7.6** Check where the logs go:
      ```bash
      journalctl -u pyirena-zmq -n 30
      tail -20 ~/.pyirena/logs/zmq.log   # as the service user
      ```

---

## Part 8 — Hand over to the orchestrator team

- [ ] **8.1** Give them the address, `tcp://SERVER:9865`, and
      [zmq_service.md](zmq_service.md).

- [ ] **8.2** Tell them the three things that are easy to get wrong:
      1. **Set the client timeout above 55 s.** 60 s is the intended pairing.
      2. **A REQ socket that times out is dead** — close it and make a new
         one, or every later call fails for an unrelated-looking reason.
         `pyirena/zmq/client.py` shows the pattern in ~150 lines.
      3. **Check `ok`, and nothing else.** A failed fit and a malformed
         request both arrive as `ok: false` with a code.

- [ ] **8.3** Point them at `analyze` as the normal path, and the session
      tools only for fits that need iteration.

- [ ] **8.4** Ask them to call `server_info` on startup and log what it
      returns. When something behaves oddly in three months, that one line
      says which deployment and which limits were in force.

---

## Troubleshooting

| Symptom | Likely cause |
|---|---|
| Smoke test: "no reply in 5 s" | Service not running, wrong port, or the firewall. Check `ss -tlnp \| grep 9865` on the server first, then the rich rule. |
| Works on loopback, not from the client | Firewall, or the service is bound to `127.0.0.1`. `--port 9865` binds all interfaces; `--bind tcp://127.0.0.1:9865` does not. |
| `zmq.error.ZMQError: Operation cannot be accomplished in current state` | A REQ socket reused after a timeout. Close it and reconnect — this is the one in 8.2.2. |
| Every call returns `NOT_AVAILABLE_OVER_ZMQ` | Expected for `save_*`, `open_dataset` with a path, data operations and image tools. The service has no shared filesystem. |
| `TIMEOUT` on a fit that is fine in the GUI | Monte-Carlo uncertainties, or a Modeling config with many populations. Check `n_mc_runs`. |
| `SESSION_STALE` | A previous call on that session hit the budget. Close it and start again. |
| `BAD_CONFIG` from `analyze` | The config has no recognisable tool section. It wants the whole exported file, `{"_pyirena_config": …, "<tool>": …}`. |
| Results differ from the GUI | Not a transport problem. Check the notes in the reply — a clipped Q range or a skipped pre-fit — then raise it. |
| `pyirena-zmq: command not found` | Wrong environment, or installed without the extra: `pip install -e ".[zmq]"`. |

---

## What to report back

Whatever happens, these are the useful things to write down:

- the versions from 1.5, and the server's OS
- `17 passed` / `19 passed`, or which checks failed and their output
- the timing from 4.4 (network overhead) and 5.6 (slowest real fit)
- the GUI-vs-service comparison from 5.3, per tool
- anything in `notes` you did not expect
- anything where the service went **silent** rather than answering — that is
  the one failure mode the design is supposed to make impossible
