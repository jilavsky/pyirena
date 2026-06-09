# 03 — In-GUI AI Advisor

**Status:** Draft
**Last updated:** 2026-06-09
**Depends on:** [01 API & MCP extensions](01-api-and-mcp-extensions.md) — Phase 1 only (read parameters + capture image)
**Related:** [02 Standalone AI app](02-standalone-ai-app.md) — shares LLM provider abstraction

---

## Goal

Add a small **"Ask AI advisor"** panel to the existing pyirena GUI. User
clicks a button while manually fitting; the panel grabs the current fit
parameters and a screenshot of the fit, sends them to an LLM, and displays
the AI's evaluation and suggestions.

**Low effort, high value.** User testing already confirmed that even local
models, given a fit screenshot and parameter list, give meaningful advice.
This subproject productizes that experiment.

---

## Why this is separate from Subproject 2

| Aspect | This (Subproject 3) | Standalone app (Subproject 2) |
|--------|---------------------|-------------------------------|
| Lives in | Existing pyirena GUI | Separate package |
| AI control | None — advice only | Full agentic control over fitting |
| User experience | "Press a button, get advice" | "Watch AI fit my data" |
| Effort | Small | Large |
| Dependency on Subproject 1 | Minimal (Phase 1 only) | Full (all phases) |
| When to ship | Quick win, ship early | After full foundation is ready |

The advisor is *advisory*; the standalone app is *autonomous*. They serve
different needs and can coexist.

---

## User experience

### Trigger
A new "Ask AI advisor" button in the existing pyirena fitting panel(s).
**Initial scope: Unified Fit only.** Unified Fit is the test case for the
entire advisor pattern (button, capture, configurator, LLM round-trip,
display). Once the pattern is validated, expand to other supported tools
([00 — scope](00-overall-plan.md#scope-of-pyirena-tools-supported)) in an
order to be decided.

### What happens on click
1. App captures the current fit plot as PNG bytes
2. App reads the current parameter table (names, values, fixed/free, bounds)
3. App reads any quality metrics available (χ², residuals)
4. App opens a small modal/panel showing "asking AI..."
5. Bundle is sent to the configured LLM with a system prompt explaining the
   context ("You are a SAXS analysis expert. The user is fitting a dataset
   with the Unified Fit model. Here are the current parameters and a plot
   of the fit. Suggest what to consider next.")
6. AI response is displayed in the panel (markdown rendered)
7. User reads, decides, applies suggestions *manually*

### What the user does NOT do
- Approve tool calls
- Watch the AI execute anything
- Configure complex strategies

The AI never touches the fit. It only speaks. This keeps the feature simple
and trustworthy — no risk of the AI silently breaking the user's work.

---

## Architecture

Lives inside `pyirena/gui/` as an additional panel. Imports a minimal LLM
client module.

```
pyirena/
├── gui/
│   ├── ai_advisor/
│   │   ├── __init__.py
│   │   ├── panel.py          # Qt panel: "Ask AI" button, response display
│   │   ├── prompt.py         # System prompt + per-model context
│   │   └── client.py         # Thin LLM client (or import from shared module)
│   └── unified_fit.py        # Add "Ask AI advisor" button → opens panel
└── config/
    └── ai_settings.py        # API key (keyring), model selection, custom instructions
```

### LLM client
Two options:
- **A. Inline client** (HTTP calls directly to Anthropic/OpenAI APIs). Small,
  no extra dependency surface beyond `requests`/`httpx`.
- **B. Shared client with Subproject 2** (extract the `llm/` provider layer
  into a small helper package both consume). More code reuse, more
  coordination.

**Recommendation:** start with A (inline client) for v1; refactor to B if
Subproject 2 ships and the duplication becomes painful. The two have
different needs (advisor is one-shot text generation; agent is multi-turn
tool use) so the shared surface is smaller than it looks.

---

## Configuration UI (AI configurator)

The configurator is reachable from the advisor panel (button: "Configure
AI..."). It is the central place where users set up their LLM access and
add domain context. Should be designed carefully — pyirena's user base is
diverse (different institutions, different LLM access).

### Sections

**1. Provider & model**
- Provider (radio): **Anthropic / OpenAI / Local (OpenAI-compatible endpoint)**
- Endpoint URL (only shown for Local; e.g. `http://localhost:1234/v1` for
  LM Studio, `http://localhost:11434/v1` for Ollama)
- Model name (text field with provider-appropriate defaults; e.g.
  `claude-opus-4-7` / `gpt-4o` / `gemma-3-...`)
- API key (password field; stored via OS keyring, never plaintext; not
  required for local)
- Test connection button (sends a trivial prompt; confirms credentials and
  model name work)

**2. Cost preview**
- Show estimated cost per advisor request given current provider/model
  (e.g. "~$0.01 per request" or "Free (local)")
- Show running total for the current session

**3. Instructions (layered, all combined when sending)**
- **Tool prompt** (read-only, app-supplied): "You are advising the user on
  use of the *Unified Fit* model. Here is the current parameter table and
  a plot of the fit..." — added by the app per fitting panel; user cannot
  edit
- **User general instructions** (editable textarea): persistent across all
  advisor invocations. Example: "Always note if χ² < 2; flag any parameter
  close to its bound"
- **Sample / project context** (editable textarea): the user's
  domain-specific framing that narrows the AI's response. Example: "We are
  dealing with polymer samples based on PEG-PCL diblock copolymers, swollen
  in water. Expected feature sizes: 5-50 nm." This is the most important
  field for fit quality — it gives the AI the scientific context the plot
  alone can't convey.

### Why three instruction layers
- App-supplied tool prompt = consistent baseline for each panel
- User general instructions = personal preferences across all fits
- Sample/project context = project-specific knowledge, often changes
  between datasets

Combined at send time as a single system prompt. User sees the combined
result in a preview if they ask for it.

---

## Image format sent to LLM

User testing confirmed JPG works with local Gemma. Need a default that is:
- High-fidelity for thin plot lines (SAXS plots have lots of fine detail)
- Compatible with all three providers (Anthropic, OpenAI, local)
- Reasonable on token cost

### Decision: PNG at 1024×768

**Why PNG over JPEG**: SAXS log-log plots are thin lines on near-uniform
backgrounds. JPEG quantization adds visible artifacts near sharp lines
(can confuse AI about residual scatter); PNG keeps lines crisp.

**Why 1024×768**: matches the typical pyirena plot aspect ratio; both
Anthropic and OpenAI charge vision tokens based on image dimensions, not
file size, so this resolution costs ~750-800 tokens regardless of format.
Larger images cost more without improving advice quality. Smaller (e.g.
768×576) is a config option for cost-sensitive users.

**Why not WebP**: support is patchy in some local model runtimes.

**Configurable** in the AI configurator under an "Advanced" section:
- Image format: PNG (default) / JPEG (smaller bandwidth)
- Image max dimension: 1024 (default) / 768 (cheaper) / 1536 (more detail)

### Cost notes (June 2026 pricing reference)
- Claude Opus 4.7 vision: ~750 tokens for 1024×768 image ≈ $0.011 per
  advisor request including text
- OpenAI gpt-4o vision (high detail): ~765 tokens for 1024×768 image ≈
  $0.005 per request
- Local: free (compute only)

Cost should be displayed in the configurator so users have realistic
expectations. Update if pricing changes.

---

## Phasing

### Phase 1 — MVP (small-medium)
- "Ask AI advisor" button on the **Unified Fit panel only** (test case)
- **Multi-LLM from day one**: Anthropic, OpenAI, and local
  (OpenAI-compatible endpoint, e.g. LM Studio / Ollama). Many users have
  institutional or convenience constraints — supporting all three early
  removes a major adoption blocker for the advisor.
- Configurator with all three instruction layers (tool prompt, user
  general instructions, sample/project context)
- API key storage via OS keyring; env var fallback
- Test connection button
- Cost preview and running total
- PNG image capture at 1024×768
- Markdown response display

Ships when Unified Fit is validated end-to-end. The configurator design
work pays off here even at MVP scope.

### Phase 2 — Roll out to other supported tools (small per tool)
- Add advisor button to each remaining supported tool, in an order TBD:
  Modeling, Size distribution, Simple fits, WAXS
- Each tool gets its own app-supplied tool prompt
- Utility tools (Scattering contrast, Data merge, Data manipulation) may
  not need an advisor button — evaluate per tool

### Phase 3 — Optional polish (small)
- Persistent advice history (save AI responses to HDF5 alongside the fit;
  "Show previous AI advice" button)
- Could feed into Subproject 2's audit trail design
- Multi-turn follow-up (user can ask the AI to elaborate without
  re-sending the image)

---

## What we already validated

User has confirmed (experimentally, with a manually-taken screenshot and
parameter list passed to local model):
- LLMs *do* give useful, specific advice on SAXS fits
- Even a local model is good enough for first-pass advice
- The user-perceived value is real

This subproject just removes the manual screenshot-and-paste friction.

---

## Resolved decisions

| Decision | Rationale |
|----------|-----------|
| **Unified Fit is the first and only Phase 1 target** | Used as the test case for the entire advisor pattern. Other tools added in Phase 2 once the pattern is validated. |
| **Image format: PNG at 1024×768** (configurable) | Best fidelity for thin plot lines; compatible with all providers; ~750 tokens regardless of provider. JPEG fallback available in advanced config. |
| **Multi-LLM in Phase 1** (Anthropic, OpenAI, local) | Many users have institutional or convenience constraints on LLM choice; advisor adoption depends on flexibility from day one. |
| **API key via OS keyring in Phase 1** (env var fallback) | Security baseline; not worth deferring to Phase 2. |
| **Three-layer prompt** (app tool prompt + user general + project context) | Tool prompt is consistent; user general carries personal preferences; project context is the highest-value user input for fit-specific advice. |
| **Cost transparency in Phase 1** | Per-request and per-session running total shown in configurator. |

## Open questions

| Question | Notes |
|----------|-------|
| **Configurator placement** | Standalone "Configure AI..." dialog reachable from the advisor button, or a new tab in existing pyirena preferences? Standalone is simpler; pref-tab integrates better. Decide during UI design. |
| **Order of Phase 2 rollout** | After Unified Fit, which tool next? Modeling has most parameters → highest value but most work. Size distribution might be a safer second. Decide based on Unified Fit experience and user requests. |
| **Local model endpoint defaults** | LM Studio default is `http://localhost:1234/v1`; Ollama is `http://localhost:11434/v1`. Provide both as quick-select buttons in configurator. |
| **Persistent advice history** | Phase 3 if at all — defer; revisit when Subproject 2 audit-trail design is settled, to avoid two incompatible formats. |

---

## Cross-cutting concerns

(From `00-overall-plan.md` — only the ones that apply.)

- **Multi-LLM support** — Phase 1 (Anthropic, OpenAI, local OpenAI-compatible)
- **Audit trail** — only if Phase 3 ships; otherwise out of scope
- **Cost transparency** — Phase 1 (per-request + session total in configurator)
- **Custom instructions per user/lab** — Phase 1 (the three-layer prompt)
- **API key management** — Phase 1 (OS keyring, env var fallback)
- **Human-in-the-loop checkpoints** — not applicable; the AI never acts

---

## Out of scope

- AI executing fits or changing parameters (that's Subproject 2)
- Multi-turn conversation (each click is independent unless user asks for
  follow-up — defer that)
- Comparing fits across files (Subproject 2 territory)
- Strategy save/replay

---

## Success criteria

1. User clicks "Ask AI advisor" while fitting a dataset in the Unified Fit
   panel
2. Within a few seconds, useful, specific suggestions appear in a side panel
3. The user applies one or more suggestions and the fit improves
4. The feature is configurable enough for users with different LLM
   preferences (provider, model, custom instructions)
5. No risk of the AI silently breaking the user's in-progress work
