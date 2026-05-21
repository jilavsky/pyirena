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
└──────────────────┘                            │   └ plotting     │
                                                 └────────┬─────────┘
                                                          │
                                                          ▼
                                                ┌──────────────────┐
                                                │  NXcanSAS HDF5   │
                                                │  files on disk   │
                                                └──────────────────┘
```

The MCP server is a small process that the AI client spawns on demand. It
exposes 18 tools (discovery, per-tool result reading, parameter
aggregation across files, headless plotting), all prefixed `pyirena_`
(e.g. `pyirena_summarize_folder`, `pyirena_list_files`) — see
[ai_tools_reference.md](ai_tools_reference.md).

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
| Plots don't render inline in the AI client | Client doesn't support inline image content from this tool | The PNG is still on disk under `PYIRENA_PLOT_CACHE` — ask the agent to give you the path |
| Wrong python is used | `pyirena-mcp` resolves to a different env | `head -1 $(which pyirena-mcp)` should show the python interpreter; if wrong, prepend the conda env's `bin/` to `PATH` or use a different absolute path |
| AI calls succeed but answers are wrong | Model is the bottleneck (limited tool-use training) | Try a model trained for tool use: Claude, GPT-4, Llama 3.1+ Instruct, Qwen 2.5+, Gemma 2 Instruct |

---

## Security model

- Read-only: the v0.7 API has no write functions. The AI cannot modify or
  delete files.
- File-access boundary: `PYIRENA_DATA_ROOT` is enforced on every public
  call. Without it set, any absolute path is accepted — use this only on
  a fully trusted client (e.g. local CLI), never when exposing the server
  to a remote agent.
- Stdio transport: the MCP process inherits its parent client's
  credentials. There is no separate auth layer in v0.7.
- Array size bounding: large arrays are decimated to
  `PYIRENA_MAX_ARRAY_POINTS` to prevent context-window blowup or
  denial-of-service via huge response payloads.

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
```

See [pyirena/api/README.md](../pyirena/api/README.md) for a complete
function list with example return shapes.
