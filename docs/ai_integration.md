# AI Integration (MCP) — Installation & Configuration

This guide explains how to expose pyirena's analysis results to AI
assistants via the [Model Context Protocol](https://modelcontextprotocol.io/).
With this set up, you can ask an AI assistant questions like:

- "Summarise what's in `/data/run42`."
- "Plot I(Q) for the three latest scans of sample_A."
- "Tabulate Rg from the Unified Fit level 1 across the folder — is it
  trending up?"
- "What was the volume fraction in the size distribution fit of
  `scan_017.h5`?"

The AI calls pyirena's API; pyirena reads the NXcanSAS HDF5 files; the
AI summarises in plain language.

> **Companion document:** [ai_tools_reference.md](ai_tools_reference.md)
> describes each MCP tool and is useful as system-prompt material for AI
> agents.

---

## Architecture in 30 seconds

```
┌──────────────────┐       stdio JSON-RPC      ┌──────────────────┐
│  AI client       │  ◄───────────────────────►│ pyirena-mcp      │
│  (Claude         │                            │ stdio server     │
│   Desktop /      │                            │                  │
│   Claude Code /  │                            │   pyirena.api    │
│   AnythingLLM /  │                            │   ├ discovery    │
│   custom agent)  │                            │   ├ readers      │
│                  │                            │   ├ aggregation  │
│                  │                            │   ├ plotting     │
│                  │                            │   ├ calculators  │
└──────────────────┘                            │   └ data ops     │
                                                 └────────┬─────────┘
                                                          │
                                                          ▼
                                                ┌──────────────────┐
                                                │  NXcanSAS HDF5   │
                                                │  files on disk   │
                                                └──────────────────┘
```

The MCP server is a small process that the AI client spawns on demand. It
exposes four families of tools — see
[ai_tools_reference.md](ai_tools_reference.md):

- **Read-only tools** (`pyirena_` prefix) — discovery, per-tool result
  reading, parameter aggregation across files, headless plotting
  (e.g. `pyirena_summarize_folder`, `pyirena_list_files`).
- **Control tools** — drive fitting interactively and write results back
  to HDF5. **Every fitting tool** is covered: Unified Fit, Size
  Distribution, Simple Fits, Modeling and WAXS Peak Fit. Session
  lifecycle (`pyirena_ctrl_open_dataset` and friends) is its own MCP
  tool; everything else — model selection, parameters, fit execution,
  quality, persistence (`*_save_fit`) — is reached through a small fixed
  dispatcher (`pyirena_list_categories` / `pyirena_list_tools` /
  `pyirena_describe_tool` / `pyirena_call`) instead of one MCP tool per
  function, so the server's registered tool count stays small (~26)
  regardless of how many control functions exist. These tools are
  stateful (session-based) and, unlike the read-only tools, can modify
  files.
- **Calculators** — stateless support calculations that need no dataset
  and open no session: scattering contrast and scattering length
  densities from chemical formulas and densities, anomalous contrast and
  transmission at a given energy, contrast-vs-energy scans for anomalous
  SAXS planning, element lookup, and read-only access to the user's saved
  compound library. Reached through the same dispatcher, as category
  `calculators`, so they add no registered MCP tools at all. They need
  the `pyirena[contrast]` extra (included in `pyirena[mcp]`).
- **Data operations** — average, subtract, divide, scale, trim, rebin and
  merge datasets. Reached through the dispatcher as category `data`, so
  again no extra registered tools. **These create new data files**: output
  goes to a sibling of the source folder (`/data/run42` →
  `/data/run42_manip`, or `_merged`) with a per-operation filename suffix,
  matching what the GUI and `pyirena.batch` already do.

---

## Installation

```bash
pip install pyirena[mcp]
```

This pulls in the MCP SDK + matplotlib (for headless plotting). It adds
one new CLI entry point: `pyirena-mcp`.

Verify:

```bash
which pyirena-mcp            # macOS / Linux
where pyirena-mcp             # Windows PowerShell
```

Note the full path — most macOS GUI clients (Claude Desktop, AnythingLLM)
do **not** inherit your shell's `PATH`, so you'll need this absolute path
in the client config.

---

## Environment variables

pyirena-mcp reads three optional environment variables:

| Variable | Purpose | Default |
|----------|---------|---------|
| `PYIRENA_DATA_ROOT` | Restrict all file access to this directory subtree. Strongly recommended when exposing the server to an AI agent. | none (any absolute path accepted) |
| `PYIRENA_MAX_ARRAY_POINTS` | Decimation cap for arrays returned in tool responses. Lower = less context bloat. | `500` |
| `PYIRENA_PLOT_CACHE` | Where generated plot PNGs are written. | `<tempdir>/pyirena-mcp` |

### Where to set them

For an AI-client-spawned MCP server, put env vars **inside the client's
MCP config JSON `env` block** rather than your shell profile. macOS and
Windows GUI apps don't reliably inherit shell environment — putting them
in the JSON guarantees they take effect. Examples below.

If you also use `pyirena-mcp` from the terminal, you can set them in
`~/.zshrc` / `~/.bashrc` (macOS/Linux) or via `System → Environment
Variables` (Windows). Both work in parallel.

---

## Client configuration

### Claude Desktop

Config file location:

| OS | Path |
|---|---|
| macOS | `~/Library/Application Support/Claude/claude_desktop_config.json` |
| Windows | `%APPDATA%\Claude\claude_desktop_config.json` |

Add `pyirena` to `mcpServers`:

```json
{
  "mcpServers": {
    "pyirena": {
      "command": "/Users/you/miniconda3/envs/pyirena/bin/pyirena-mcp",
      "env": {
        "PYIRENA_DATA_ROOT": "/Users/you/data/saxs",
        "PYIRENA_MAX_ARRAY_POINTS": "500"
      }
    }
  }
}
```

Restart Claude Desktop. Click the 🔌 / tools icon — `pyirena` should
appear with 18 tools.

### Claude Code

In the project directory:

```bash
claude mcp add pyirena \
  /Users/you/miniconda3/envs/pyirena/bin/pyirena-mcp \
  --env PYIRENA_DATA_ROOT=/Users/you/data/saxs
```

Or edit the project's `.mcp.json` directly using the same JSON schema as
Claude Desktop.

### AnythingLLM Desktop (macOS)

Config file: `~/Library/Application Support/anythingllm-desktop/storage/plugins/anythingllm_mcp_servers.json`

This file is auto-created the first time you open the **Agent Skills**
page in AnythingLLM. Then quit the app and edit:

```json
{
  "mcpServers": {
    "pyirena": {
      "command": "/Users/you/miniconda3/envs/pyirena/bin/pyirena-mcp",
      "args": [],
      "env": {
        "PYIRENA_DATA_ROOT": "/Users/you/data/saxs",
        "PYIRENA_PLOT_CACHE": "/Users/you/Library/Caches/pyirena-mcp",
        "PYIRENA_MAX_ARRAY_POINTS": "500"
      },
      "anythingllm": {
        "autoStart": true
      }
    }
  }
}
```

Restart AnythingLLM (or hit **Refresh** in Agent Skills). The pyirena
tools should be listed.

> **Important on AnythingLLM:** the `command` path **must be absolute**.
> AnythingLLM does not inherit your shell `PATH`, so `pyirena-mcp` alone
> will fail with "command not found". Always paste the output of
> `which pyirena-mcp` into the JSON.

### Generic custom agent

Any MCP-capable client that supports stdio servers can use pyirena. The
launch command is just:

```
/abs/path/to/pyirena-mcp
```

with env vars set per the client's convention.

---

## Verification

Three increasing levels of confidence.

### 1. Smoke check — server starts

```bash
pyirena-mcp
```

It will read from stdin and wait. Press Ctrl-C. No traceback = healthy.

### 2. Manual JSON-RPC handshake

Server is up, your stdin sends one JSON message per line:

```bash
pyirena-mcp
```

Paste (one line):

```json
{"jsonrpc":"2.0","id":1,"method":"initialize","params":{"protocolVersion":"2024-11-05","capabilities":{},"clientInfo":{"name":"manual","version":"1.0"}}}
```

You should get a JSON response with server capabilities. Then:

```json
{"jsonrpc":"2.0","method":"notifications/initialized"}
{"jsonrpc":"2.0","id":2,"method":"tools/list"}
```

The `tools/list` response should enumerate 18 pyirena tools.

### 3. Browser inspector (recommended for interactive testing)

```bash
npx @modelcontextprotocol/inspector pyirena-mcp
```

Opens a web UI where you can call any tool with a form-built input and
see the JSON response. Requires Node/npm; the inspector itself is
fetched on the fly.

---

## Troubleshooting

| Symptom | Likely cause | Fix |
|---|---|---|
| Server doesn't appear in the client | `command` not found by the client | Use the absolute path from `which pyirena-mcp` |
| `command not found` in client logs | macOS GUI apps don't inherit shell `PATH` | Same — absolute path |
| Tool calls fail with `PathSecurityError` | File is outside `PYIRENA_DATA_ROOT` | Either widen the root or remove the env var while debugging |
| File-read tools return `{"found": false, ...}` | Group not present in the file, or wrong path | Use `inspect_file()` first to see which analyses are present |
| Plots don't render inline in the AI client | Client doesn't support MCP `ImageContent` blocks | The PNG is on disk under `PYIRENA_PLOT_CACHE` — ask the agent for the path from the text item (content[0]) |
| Agent prints a wall of garbled characters when a plot is requested | Agent/pipeline iterated over tool content items without checking `type`, stringifying the raw base64 PNG data | Instruct the agent: plot tools return **two content items** — `text` (file path) and `image` (base64 PNG). It must branch on `item.type`; never print or forward an `image` item as text. See [ai_tools_reference.md § Plotting](ai_tools_reference.md#plotting-returns-mixed-text--image-content) |
| Wrong python is used | `pyirena-mcp` resolves to a different env | `head -1 $(which pyirena-mcp)` should show the python interpreter; if wrong, prepend the conda env's `bin/` to `PATH` or use a different absolute path |
| AI calls succeed but answers are wrong | Model is the bottleneck (limited tool-use training) | Try a model trained for tool use: Claude, GPT-4, Llama 3.1+ Instruct, Qwen 2.5+, Gemma 2 Instruct |

---

## Security model

- **Not read-only.** Two groups write: the control tools save fit results
  into NXcanSAS files (in place by default), and the `data` operations
  create new data files next to the source. Nothing deletes or overwrites
  an input dataset — a data operation always writes to a new name — but a
  repeated operation does overwrite its own previous output, and saving a
  fit without an explicit `output_path` updates the source file in place.
- File-access boundary: `PYIRENA_DATA_ROOT` is enforced on every public
  call. Without it set, any absolute path is accepted — use this only on
  a fully trusted client (e.g. local CLI), never when exposing the server
  to a remote agent.
- Stdio transport: the MCP process inherits its parent client's
  credentials. There is no separate auth layer in v0.7.
- Array size bounding: large arrays are decimated to
  `PYIRENA_MAX_ARRAY_POINTS` to prevent context-window blowup or
  denial-of-service via huge response payloads.
- Calculators sit outside the `PYIRENA_DATA_ROOT` boundary because they
  take no path: they compute from formulas and densities alone. The one
  file they touch is the user's own compound library in their home
  directory, and only for reading — saving and deleting are not exposed,
  so an agent cannot mutate it.
- Data operations honour `PYIRENA_DATA_ROOT` on both the inputs and the
  derived output folder. Note the default output folder is a *sibling* of
  the source folder, so setting the root to the data folder itself makes
  the default illegal; the tool then returns `PATH_NOT_ALLOWED` and the
  agent must pass an explicit in-root `output_folder`. Set the root one
  level above your data to avoid this.

---

## Programmatic use (no MCP)

The same surface is usable as a regular Python library:

```python
from pyirena import api

api.summarize_folder("/data/run42")
api.read_unified_fit("/data/run42/sample_A_scan_017.h5")
api.tabulate_parameter("/data/run42", tool="unified_fit",
                        parameter="Rg", subgroup_index=1)
api.plot_iq(["/data/run42/sample_A_scan_017.h5",
             "/data/run42/sample_A_scan_018.h5"],
            output_path="/tmp/iq.png")

# Calculators need no data at all
api.calc_contrast("TiO", 4.95, "Ti2O3", 4.49)["xray_contrast"]  # 11.63

# Data operations write a new file and return where it went
api.average_data(["/data/run42/f001.h5", "/data/run42/f002.h5"])
api.subtract_data("/data/run42/sample.h5", "/data/run42/buffer.h5")
api.merge_datasets("/data/usaxs/s_001.h5", "/data/saxs/s_001.h5")
```

See [pyirena/api/README.md](../pyirena/api/README.md) for a complete
function list with example return shapes.
