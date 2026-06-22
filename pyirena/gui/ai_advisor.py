"""AI advisor for pyirena fitting panels.

Provides:
  - AiAdvisorConfigDialog  — provider / model / API key / instructions settings
  - AiAdvisorThread        — background LLM call (Anthropic or local OpenAI-compat)
  - AiAdvisorResultPanel   — non-modal result display window
  - launch_unified_fit_advisor(panel) — entry point called by the button handler

Phase 1 supports:
  - Anthropic (Claude) via the official `anthropic` SDK
  - Local OpenAI-compatible endpoint (LM Studio, Ollama) via httpx

API keys are stored in the OS keyring (never in state.json).
Fallback: ANTHROPIC_API_KEY / OPENAI_API_KEY environment variables.
"""
from __future__ import annotations

import base64
import html as _html_mod
import os
import re
import traceback
from pathlib import Path
from typing import Optional

try:
    from PySide6.QtCore import Qt, QThread, QTimer, Signal
    from PySide6.QtGui import QFont
    from PySide6.QtWidgets import (
        QApplication, QDialog, QDialogButtonBox, QFormLayout, QGroupBox,
        QHBoxLayout, QLabel, QLineEdit, QPlainTextEdit, QPushButton,
        QRadioButton, QSizePolicy, QTextBrowser, QTextEdit, QVBoxLayout,
        QWidget,
    )
    from PySide6.QtCore import QBuffer, QByteArray, QIODevice
except ImportError:
    try:
        from PyQt6.QtCore import Qt, QThread, QTimer, pyqtSignal as Signal
        from PyQt6.QtGui import QFont
        from PyQt6.QtWidgets import (
            QApplication, QDialog, QDialogButtonBox, QFormLayout, QGroupBox,
            QHBoxLayout, QLabel, QLineEdit, QPlainTextEdit, QPushButton,
            QRadioButton, QSizePolicy, QTextBrowser, QTextEdit, QVBoxLayout,
            QWidget,
        )
        from PyQt6.QtCore import QBuffer, QByteArray, QIODevice
    except ImportError:
        from PyQt5.QtCore import Qt, QThread, QTimer, pyqtSignal as Signal
        from PyQt5.QtGui import QFont
        from PyQt5.QtWidgets import (
            QApplication, QDialog, QDialogButtonBox, QFormLayout, QGroupBox,
            QHBoxLayout, QLabel, QLineEdit, QPlainTextEdit, QPushButton,
            QRadioButton, QSizePolicy, QTextBrowser, QTextEdit, QVBoxLayout,
            QWidget,
        )
        from PyQt5.QtCore import QBuffer, QByteArray, QIODevice


# ---------------------------------------------------------------------------
# Response formatting: LaTeX math → Unicode + Markdown → HTML
# ---------------------------------------------------------------------------

_SUP = str.maketrans("0123456789+-nxi", "⁰¹²³⁴⁵⁶⁷⁸⁹⁺⁻ⁿˣⁱ")
_SUB = str.maketrans("0123456789+-", "₀₁₂₃₄₅₆₇₈₉₊₋")

_LATEX_SYMBOLS = [
    # Angstrom (must come before generic \text{} rule)
    (r"\\text\{\\AA\}", "Å"), (r"\\text\{Å\}", "Å"), (r"\\AA\b", "Å"),
    # Remove \text{} wrapper, keep content
    (r"\\text\{([^}]+)\}", r"\1"),
    # Math operators
    (r"\\approx\b", "≈"), (r"\\times\b", "×"), (r"\\cdot\b", "·"),
    (r"\\leq\b", "≤"),    (r"\\geq\b", "≥"),   (r"\\neq\b", "≠"),
    (r"\\pm\b", "±"),     (r"\\infty\b", "∞"),
    # Greek lower
    (r"\\alpha\b", "α"), (r"\\beta\b", "β"),  (r"\\gamma\b", "γ"),
    (r"\\delta\b", "δ"), (r"\\epsilon\b", "ε"),(r"\\zeta\b", "ζ"),
    (r"\\eta\b", "η"),   (r"\\theta\b", "θ"),  (r"\\kappa\b", "κ"),
    (r"\\lambda\b", "λ"),(r"\\mu\b", "μ"),     (r"\\nu\b", "ν"),
    (r"\\xi\b", "ξ"),    (r"\\pi\b", "π"),     (r"\\rho\b", "ρ"),
    (r"\\sigma\b", "σ"), (r"\\tau\b", "τ"),    (r"\\phi\b", "φ"),
    (r"\\psi\b", "ψ"),   (r"\\omega\b", "ω"),  (r"\\chi\b", "χ"),
    # Greek upper
    (r"\\Gamma\b", "Γ"), (r"\\Delta\b", "Δ"),  (r"\\Sigma\b", "Σ"),
    (r"\\Omega\b", "Ω"), (r"\\Lambda\b", "Λ"),
]


def _convert_latex(s: str) -> str:
    """Convert common LaTeX math notation to plain Unicode."""
    for pattern, repl in _LATEX_SYMBOLS:
        s = re.sub(pattern, repl, s)

    # Superscripts: ^{...}  or  ^digit/sign
    def _sup_block(m):
        t = m.group(1).translate(_SUP)
        return t if t != m.group(1) else f"^{m.group(1)}"
    s = re.sub(r"\^\{([^}]+)\}", _sup_block, s)
    s = re.sub(r"\^([0-9+\-nxi])", lambda m: m.group(1).translate(_SUP), s)

    # Subscripts: _{...}  or  _digit/sign
    def _sub_block(m):
        t = m.group(1).translate(_SUB)
        return t if t != m.group(1) else f"_{m.group(1)}"
    s = re.sub(r"_\{([^}]+)\}", _sub_block, s)
    s = re.sub(r"_([0-9+\-])", lambda m: m.group(1).translate(_SUB), s)

    return s


def _format_response(text: str) -> str:
    """Convert LLM response (Markdown + inline LaTeX) to HTML for QTextBrowser.

    Handles:
    - $...$ and $$...$$ inline math  → Unicode italic span
    - **bold** and *italic* markdown
    - Unordered bullet lists (* item or - item)
    - Plain paragraphs
    """
    # Split into alternating [plain, math, plain, math, ...] segments.
    # Match $$...$$ before $...$ so double-dollar is not consumed as two singles.
    segments = re.split(r"(\$\$[^$]+\$\$|\$[^$\n]+\$)", text)

    parts: list[str] = []
    for i, seg in enumerate(segments):
        if i % 2 == 1:  # math segment
            is_display = seg.startswith("$$")
            content = seg[2:-2] if is_display else seg[1:-1]
            content = _convert_latex(content)
            parts.append(f'<i style="color:#1a5276;">{_html_mod.escape(content)}</i>')
        else:
            parts.append(_html_mod.escape(seg))

    processed = "".join(parts)

    # Bold **...**  (may contain already-HTML math spans — use non-greedy)
    processed = re.sub(r"\*\*(.+?)\*\*", r"<b>\1</b>", processed, flags=re.DOTALL)
    # Italic *...*  (avoid matching ** that bold already consumed)
    processed = re.sub(r"(?<!\*)\*(?!\*)(.+?)(?<!\*)\*(?!\*)", r"<i>\1</i>", processed)

    # Build HTML from lines
    lines = processed.split("\n")
    html_parts: list[str] = []
    in_list = False

    for line in lines:
        stripped = line.strip()
        is_bullet = stripped.startswith("* ") or stripped.startswith("- ")

        if is_bullet:
            if not in_list:
                html_parts.append('<ul style="margin:4px 0 4px 0;padding-left:20px;">')
                in_list = True
            html_parts.append(
                f'<li style="margin-bottom:5px;">{stripped[2:]}</li>'
            )
        else:
            if in_list:
                html_parts.append("</ul>")
                in_list = False
            if stripped:
                html_parts.append(f'<p style="margin:4px 0;">{line}</p>')

    if in_list:
        html_parts.append("</ul>")

    return "".join(html_parts)


# ---------------------------------------------------------------------------
# API key helpers (OS keyring ↔ env-var fallback)
# ---------------------------------------------------------------------------

_KEYRING_SERVICE = "pyirena-ai"
_KEYRING_KEYS = {
    "anthropic": "anthropic_api_key",
    "openai":    "openai_api_key",
    "local":     "local_api_key",
}


def _get_api_key(provider: str) -> str:
    """Return API key for *provider*, preferring OS keyring over env vars."""
    try:
        import keyring
        key = keyring.get_password(_KEYRING_SERVICE, _KEYRING_KEYS.get(provider, provider))
        if key:
            return key
    except Exception:
        pass
    # env-var fallback
    if provider == "anthropic":
        return os.environ.get("ANTHROPIC_API_KEY", "")
    if provider == "openai":
        return os.environ.get("OPENAI_API_KEY", "")
    return ""


def _set_api_key(provider: str, key: str) -> None:
    """Store API key in OS keyring."""
    try:
        import keyring
        keyring.set_password(_KEYRING_SERVICE, _KEYRING_KEYS.get(provider, provider), key)
    except Exception:
        pass  # silently fall back to env var usage


# ---------------------------------------------------------------------------
# Skills file loading (expert fitting guidance, per tool)
# ---------------------------------------------------------------------------

def _load_skills(tool_key: str) -> str:
    """Load expert fitting guidance for a tool.

    Search order:
    1. ~/.pyirena/ai_skills/<tool_key>.md  (user override)
    2. pyirena/gui/ai_skills/<tool_key>.md (bundled default)
    """
    user_path = Path.home() / ".pyirena" / "ai_skills" / f"{tool_key}.md"
    if user_path.exists():
        try:
            return user_path.read_text(encoding="utf-8")
        except Exception:
            pass
    pkg_path = Path(__file__).parent / "ai_skills" / f"{tool_key}.md"
    if pkg_path.exists():
        try:
            return pkg_path.read_text(encoding="utf-8")
        except Exception:
            pass
    return ""


# ---------------------------------------------------------------------------
# Background LLM worker
# ---------------------------------------------------------------------------

_TOOL_SYSTEM_PROMPT = """\
You are advising a SAXS/USAXS scientist on their Unified Fit model analysis.

The Unified Fit model (Beaucage formalism) describes hierarchical scattering \
structures. Each structural level contributes a Guinier term (radius of \
gyration Rg, amplitude G) and a power-law term (slope P, prefactor B). \
Multiple levels represent multi-scale structures. Background is a flat \
incoherent offset.

You will receive:
1. A screenshot of the fit panel showing the data, current model, and \
residuals (if a fit has been run).
2. A formatted table of current parameter values, bounds, and fit/fixed status.
3. Fit quality (reduced χ², if available).

Give specific, actionable advice on what to try next. Be concise — 3 to 8 \
bullet points. Reference specific parameter names (e.g. "try freeing P_1", \
"Rg_1 appears too large given the Guinier knee position"). Do not give \
generic fitting advice that ignores the specific values shown."""

# Approximate token costs ($/million tokens): (input, output)
_COST_PER_1M: dict[str, tuple[float, float]] = {
    "claude-opus-4-7":          (15.0,  75.0),
    "claude-sonnet-4-6":        ( 3.0,  15.0),
    "claude-haiku-4-5-20251001":( 0.8,   4.0),
    "gpt-4o":                   ( 2.5,  10.0),
    "gpt-4o-mini":              ( 0.15,  0.60),
    "o3":                       (10.0,  40.0),
    "o4-mini":                  ( 1.1,   4.4),
}


def _estimate_cost(model: str, input_tok: int, output_tok: int) -> Optional[float]:
    key = model.lower()
    for k, (ci, co) in _COST_PER_1M.items():
        if k in key:
            return (input_tok * ci + output_tok * co) / 1_000_000
    return None  # unknown model — don't show estimate


def _build_system_prompt(tool_key: str, user_instructions: str,
                          project_context: str) -> str:
    parts = [_TOOL_SYSTEM_PROMPT]
    skills = _load_skills(tool_key)
    if skills.strip():
        parts.append("\n## Expert fitting guidance for this tool\n" + skills.strip())
    if user_instructions.strip():
        parts.append("\nAdditional user preferences:\n" + user_instructions.strip())
    if project_context.strip():
        parts.append("\nSample/project context:\n" + project_context.strip())
    return "\n".join(parts)


def _format_param_table(levels: list[dict], background: float,
                         num_levels: int, chi_sq: Optional[float],
                         file_path: Optional[str]) -> str:
    lines = ["## Current Unified Fit State"]
    if file_path:
        lines.append(f"File: {Path(file_path).name}")
    lines.append(f"Levels: {num_levels}")
    if chi_sq is not None:
        lines.append(f"Reduced χ²: {chi_sq:.4f}")
    else:
        lines.append("Reduced χ²: no fit run yet")
    lines.append("")

    for i, lv in enumerate(levels, 1):
        lines.append(f"### Level {i}")
        correlated = lv.get("correlated", False)
        lines.append("| Parameter | Value | Fixed | Lo | Hi |")
        lines.append("|-----------|-------|-------|----|-----|")
        # ETA and PACK are only active when the Correlations checkbox is checked.
        # Omit them entirely when unused so stale values don't mislead the AI.
        active_params = ("Rg", "G", "B", "P", "RgCutoff")
        if correlated:
            active_params = ("Rg", "G", "B", "P", "ETA", "PACK", "RgCutoff")
        for param in active_params:
            val = lv.get(param, "—")
            fixed = not lv.get(f"fit_{param}", True)
            lo = lv.get(f"{param}_low", "—")
            hi = lv.get(f"{param}_high", "—")
            fixed_str = "Yes" if fixed else "No"
            val_s  = f"{val:.4g}"  if isinstance(val, float) else str(val)
            lo_s   = f"{lo:.4g}"   if isinstance(lo, float)  else str(lo)
            hi_s   = f"{hi:.4g}"   if isinstance(hi, float)  else str(hi)
            lines.append(f"| {param} | {val_s} | {fixed_str} | {lo_s} | {hi_s} |")
        extras = []
        if correlated:
            extras.append("correlations ON (ETA, PACK active)")
        else:
            extras.append("correlations OFF (ETA/PACK not used in fit)")
        if lv.get("estimate_B"):
            extras.append("estimate_B ON")
        if lv.get("link_rgco"):
            extras.append("link_RgCO ON")
        lines.append(f"Flags: {', '.join(extras)}")
        lines.append("")

    lines.append("### Background")
    lines.append(f"Value: {background:.4g}  |  Fixed: No")
    return "\n".join(lines)


class AiAdvisorThread(QThread):
    """Background thread that calls the LLM and emits the response + usage."""

    # (response_text, usage_dict) where usage_dict has input_tokens, output_tokens
    result_ready   = Signal(str, dict)
    error_occurred = Signal(str)

    def __init__(
        self,
        provider:      str,
        api_key:       str,
        model:         str,
        endpoint:      str,   # base URL for all providers
        system_prompt: str,
        param_text:    str,
        img_b64:       str,
        parent=None,
    ):
        super().__init__(parent)
        self.provider      = provider
        self.api_key       = api_key
        self.model         = model
        self.endpoint      = endpoint.rstrip("/")
        self.system_prompt = system_prompt
        self.param_text    = param_text
        self.img_b64       = img_b64

    def run(self):
        try:
            if self.provider == "anthropic":
                text, usage = self._call_anthropic()
            else:  # "openai" or "local" — both use OpenAI-compatible format
                text, usage = self._call_openai_compat()
            self.result_ready.emit(text, usage)
        except Exception as exc:
            self.error_occurred.emit(f"{type(exc).__name__}: {exc}\n\n{traceback.format_exc()}")

    def _call_anthropic(self) -> tuple[str, dict]:
        import anthropic  # noqa: PLC0415
        kwargs: dict = {"api_key": self.api_key}
        if self.endpoint:
            kwargs["base_url"] = self.endpoint
        client = anthropic.Anthropic(**kwargs)
        content = []
        if self.img_b64:
            content.append({
                "type": "image",
                "source": {"type": "base64", "media_type": "image/png",
                           "data": self.img_b64},
            })
        content.append({"type": "text", "text": self.param_text})
        msg = client.messages.create(
            model=self.model, max_tokens=1024,
            system=self.system_prompt,
            messages=[{"role": "user", "content": content}],
        )
        usage = {
            "input_tokens":  msg.usage.input_tokens,
            "output_tokens": msg.usage.output_tokens,
        }
        return msg.content[0].text, usage

    def _call_openai_compat(self) -> tuple[str, dict]:
        """Call an OpenAI-compatible endpoint (OpenAI, LM Studio, Ollama)."""
        import httpx  # noqa: PLC0415
        if not self.endpoint:
            raise ValueError("No endpoint URL configured for this provider.")
        messages = [{"role": "system", "content": self.system_prompt}]
        user_parts: list = []
        if self.img_b64:
            user_parts.append({
                "type":      "image_url",
                "image_url": {"url": f"data:image/png;base64,{self.img_b64}"},
            })
        user_parts.append({"type": "text", "text": self.param_text})
        messages.append({"role": "user", "content": user_parts})
        # Local models have no token cost — allow generous output.
        max_tokens = 4096 if self.provider == "local" else 1024
        resp = httpx.post(
            f"{self.endpoint}/chat/completions",
            headers={"Authorization": f"Bearer {self.api_key or 'local'}",
                     "Content-Type": "application/json"},
            json={"model": self.model, "messages": messages, "max_tokens": max_tokens},
            timeout=120,
        )
        resp.raise_for_status()
        body = resp.json()
        choice = body["choices"][0]
        text = choice["message"].get("content") or ""
        if choice.get("finish_reason") == "length":
            text += "\n\n*(Response was cut off — the model hit the token limit.)*"
        raw_usage = body.get("usage", {})
        usage = {
            "input_tokens":  raw_usage.get("prompt_tokens", 0),
            "output_tokens": raw_usage.get("completion_tokens", 0),
        }
        return text, usage


# ---------------------------------------------------------------------------
# Result display panel (non-modal)
# ---------------------------------------------------------------------------

class AiAdvisorResultPanel(QDialog):
    """Non-modal window that shows the AI advisor response."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("AI Advisor — Unified Fit")
        self.setModal(False)
        self.resize(600, 450)
        self._thread: Optional[AiAdvisorThread] = None

        layout = QVBoxLayout(self)

        # Status / spinner label
        self._status_label = QLabel("Asking AI advisor…")
        self._status_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        status_font = QFont()
        status_font.setItalic(True)
        self._status_label.setFont(status_font)
        layout.addWidget(self._status_label)

        # Response area — QTextBrowser renders HTML (bold, italic, lists)
        self._text_edit = QTextBrowser()
        self._text_edit.setOpenExternalLinks(False)
        # Match the application's default font (same as the rest of the GUI)
        app_font = QApplication.font()
        self._text_edit.setFont(app_font)
        self._text_edit.document().setDefaultStyleSheet(
            "p { margin: 4px 0; } "
            "ul { margin: 4px 0; padding-left: 20px; } "
            "li { margin-bottom: 5px; }"
        )
        layout.addWidget(self._text_edit)

        # Token usage footer (populated after response arrives)
        self._usage_label = QLabel("")
        usage_font = QFont()
        usage_font.setPointSize(max(8, QApplication.font().pointSize() - 1))
        usage_font.setItalic(True)
        self._usage_label.setFont(usage_font)
        self._usage_label.setStyleSheet("color: #555;")
        layout.addWidget(self._usage_label)

        # Close button
        btn_box = QDialogButtonBox(QDialogButtonBox.StandardButton.Close)
        btn_box.rejected.connect(self.close)
        layout.addWidget(btn_box)

        # Spinner animation
        self._dots = 0
        self._timer = QTimer(self)
        self._timer.timeout.connect(self._animate)
        self._timer.start(400)

    def _animate(self):
        if self._text_edit.toPlainText():
            self._timer.stop()
            return
        self._dots = (self._dots + 1) % 4
        self._status_label.setText("Asking AI advisor" + "." * self._dots)

    def set_thread(self, thread: AiAdvisorThread):
        self._thread = thread
        self._model = thread.model
        thread.result_ready.connect(self._on_result)
        thread.error_occurred.connect(self._on_error)
        thread.start()

    def _on_result(self, text: str, usage: dict):
        self._timer.stop()
        self._status_label.setText("AI advisor response:")
        if not text:
            self._text_edit.setPlainText("[No response content received from the model.]")
            return
        self._text_edit.setHtml(_format_response(text))

        # Usage footer
        inp = usage.get("input_tokens", 0)
        out = usage.get("output_tokens", 0)
        cost = _estimate_cost(getattr(self, "_model", ""), inp, out)
        cost_str = f" | est. cost ~${cost:.4f}" if cost is not None else ""
        self._usage_label.setText(
            f"Tokens: {inp:,} in + {out:,} out{cost_str}"
        )

    def _on_error(self, msg: str):
        self._timer.stop()
        self._status_label.setText("Error from AI advisor:")
        self._text_edit.setPlainText(msg)

    def closeEvent(self, event):
        if self._thread and self._thread.isRunning():
            self._thread.quit()
            self._thread.wait(3000)
        super().closeEvent(event)


# ---------------------------------------------------------------------------
# Configuration dialog
# ---------------------------------------------------------------------------

class AiAdvisorConfigDialog(QDialog):
    """Settings dialog: per-provider model/URL/key + shared instructions."""

    # Default settings for each provider
    _PROVIDER_DEFAULTS = {
        "anthropic": {"model": "claude-opus-4-7",     "url": ""},
        "openai":    {"model": "gpt-4o",               "url": "https://api.openai.com/v1"},
        "local":     {"model": "gemma-3-27b-it",       "url": "http://localhost:1234/v1"},
    }
    _URL_LABELS = {
        "anthropic": ("Custom API URL:",
                      "Leave blank for default (https://api.anthropic.com)"),
        "openai":    ("API base URL:",
                      "https://api.openai.com/v1"),
        "local":     ("Endpoint URL:",
                      "http://localhost:1234/v1"),
    }

    def __init__(self, state_manager, parent=None):
        super().__init__(parent)
        self.state_manager = state_manager
        self.setWindowTitle("AI Advisor — Configure")
        self.setModal(True)
        self.setMinimumWidth(500)
        # Per-provider in-memory cache: {provider: {model, url}}
        self._cache: dict[str, dict] = {p: dict(d) for p, d in self._PROVIDER_DEFAULTS.items()}
        self._active_provider: str = "anthropic"
        self._build_ui()
        self._load_state()

    def _build_ui(self):
        layout = QVBoxLayout(self)

        # --- Provider & model group ---
        prov_group = QGroupBox("Provider and model")
        prov_layout = QVBoxLayout(prov_group)

        radio_row = QHBoxLayout()
        self._radio_anthropic = QRadioButton("Anthropic (Claude)")
        self._radio_openai    = QRadioButton("OpenAI (ChatGPT)")
        self._radio_local     = QRadioButton("Local (LM Studio / Ollama)")
        radio_row.addWidget(self._radio_anthropic)
        radio_row.addWidget(self._radio_openai)
        radio_row.addWidget(self._radio_local)
        prov_layout.addLayout(radio_row)

        form = QFormLayout()

        self._model_edit = QLineEdit()
        form.addRow("Model name:", self._model_edit)

        self._url_label = QLabel("API URL:")
        self._url_edit  = QLineEdit()
        form.addRow(self._url_label, self._url_edit)

        self._key_edit = QLineEdit()
        self._key_edit.setEchoMode(QLineEdit.EchoMode.Password)
        self._key_edit.setPlaceholderText("Leave blank to use env var")
        form.addRow("API key:", self._key_edit)

        prov_layout.addLayout(form)

        self._test_btn   = QPushButton("Test connection")
        self._test_label = QLabel("")
        test_row = QHBoxLayout()
        test_row.addWidget(self._test_btn)
        test_row.addWidget(self._test_label)
        test_row.addStretch()
        prov_layout.addLayout(test_row)

        layout.addWidget(prov_group)

        # Wire radio buttons — save current fields before switching
        self._radio_anthropic.toggled.connect(
            lambda chk: self._on_provider_radio(chk, "anthropic"))
        self._radio_openai.toggled.connect(
            lambda chk: self._on_provider_radio(chk, "openai"))
        self._radio_local.toggled.connect(
            lambda chk: self._on_provider_radio(chk, "local"))
        self._test_btn.clicked.connect(self._test_connection)

        # --- Instructions group ---
        inst_group = QGroupBox("Instructions (sent to AI with every request)")
        inst_layout = QVBoxLayout(inst_group)

        inst_layout.addWidget(QLabel("General instructions (personal preferences):"))
        self._general_edit = QPlainTextEdit()
        self._general_edit.setMaximumHeight(80)
        self._general_edit.setPlaceholderText(
            "e.g. Always note if χ² < 2; flag parameters close to their bounds."
        )
        inst_layout.addWidget(self._general_edit)

        inst_layout.addWidget(QLabel("Sample / project context (material description):"))
        self._context_edit = QPlainTextEdit()
        self._context_edit.setMaximumHeight(80)
        self._context_edit.setPlaceholderText(
            "e.g. PEG-PCL diblock copolymer in water, expected feature sizes 5–50 nm."
        )
        inst_layout.addWidget(self._context_edit)

        layout.addWidget(inst_group)

        btn_box = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel
        )
        btn_box.accepted.connect(self._save_and_accept)
        btn_box.rejected.connect(self.reject)
        layout.addWidget(btn_box)

    # ------------------------------------------------------------------

    def _provider_from_radios(self) -> str:
        if self._radio_openai.isChecked():
            return "openai"
        if self._radio_local.isChecked():
            return "local"
        return "anthropic"

    def _save_active_to_cache(self):
        """Flush current field values back into _cache for the active provider."""
        self._cache[self._active_provider] = {
            "model": self._model_edit.text().strip(),
            "url":   self._url_edit.text().strip(),
        }

    def _populate_fields_from_cache(self, provider: str):
        """Fill model/URL fields from cache; load API key from keyring."""
        c = self._cache[provider]
        self._model_edit.setText(c.get("model", ""))
        self._url_edit.setText(c.get("url", ""))
        self._key_edit.setText(_get_api_key(provider))
        label, placeholder = self._URL_LABELS[provider]
        self._url_label.setText(label)
        self._url_edit.setPlaceholderText(placeholder)

    def _on_provider_radio(self, checked: bool, provider: str):
        if not checked:
            return
        if self._active_provider != provider:
            # True provider switch — save current fields before changing
            self._save_active_to_cache()
            self._active_provider = provider
        self._populate_fields_from_cache(provider)

    def _load_state(self):
        cfg = self.state_manager.get("ai_advisor") or {}
        # Seed cache from saved state
        for p in ("anthropic", "openai", "local"):
            saved = cfg.get(p, {})
            defaults = self._PROVIDER_DEFAULTS[p]
            # "url" key covers both base_url (anthropic/openai) and endpoint (local)
            self._cache[p] = {
                "model": saved.get("model",    saved.get("endpoint",
                         defaults["model"])),
                "url":   saved.get("base_url", saved.get("endpoint",
                         defaults["url"])),
            }
        # Set radio (triggers _on_provider_radio)
        provider = cfg.get("provider", "anthropic")
        self._active_provider = provider
        if provider == "openai":
            self._radio_openai.setChecked(True)
        elif provider == "local":
            self._radio_local.setChecked(True)
        else:
            self._radio_anthropic.setChecked(True)
        self._populate_fields_from_cache(provider)
        self._general_edit.setPlainText(cfg.get("user_instructions", ""))
        self._context_edit.setPlainText(cfg.get("project_context", ""))

    def _save_and_accept(self):
        self._save_active_to_cache()
        provider = self._provider_from_radios()
        api_key  = self._key_edit.text().strip()
        if api_key:
            _set_api_key(provider, api_key)

        self.state_manager.set("ai_advisor", "provider", provider)
        # Save all three provider caches
        for p, c in self._cache.items():
            block = {"model": c["model"]}
            if p == "local":
                block["endpoint"] = c["url"]
            else:
                block["base_url"] = c["url"]
            self.state_manager.set("ai_advisor", p, block)
        self.state_manager.set("ai_advisor", "user_instructions",
                                self._general_edit.toPlainText())
        self.state_manager.set("ai_advisor", "project_context",
                                self._context_edit.toPlainText())
        self.state_manager.save()
        self.accept()

    def _test_connection(self):
        self._test_label.setText("Testing…")
        provider = self._provider_from_radios()
        api_key  = self._key_edit.text().strip() or _get_api_key(provider)
        model    = self._model_edit.text().strip()
        url      = self._url_edit.text().strip()

        class _TestThread(QThread):
            done = Signal(str)

            def __init__(self, provider, api_key, model, url, parent=None):
                super().__init__(parent)
                self._p, self._k, self._m, self._u = provider, api_key, model, url

            def run(self):
                try:
                    if self._p == "anthropic":
                        import anthropic
                        kwargs: dict = {"api_key": self._k}
                        if self._u:
                            kwargs["base_url"] = self._u
                        c = anthropic.Anthropic(**kwargs)
                        r = c.messages.create(
                            model=self._m, max_tokens=10,
                            messages=[{"role": "user", "content": "Hi"}],
                        )
                        self.done.emit(f"OK — model: {r.model}")
                    else:
                        import httpx
                        ep = self._u.rstrip("/")
                        if not ep:
                            raise ValueError("No endpoint URL set.")
                        r = httpx.post(
                            f"{ep}/chat/completions",
                            headers={"Authorization": f"Bearer {self._k or 'local'}"},
                            json={"model": self._m,
                                  "messages": [{"role": "user", "content": "Hi"}],
                                  "max_tokens": 10},
                            timeout=15,
                        )
                        r.raise_for_status()
                        self.done.emit(f"OK — {r.json()['model']}")
                except Exception as exc:
                    self.done.emit(f"Error: {exc}")

        self._test_thread = _TestThread(provider, api_key, model, url, self)
        self._test_thread.done.connect(self._test_label.setText)
        self._test_thread.start()


# ---------------------------------------------------------------------------
# Entry point: collect state and launch advisor
# ---------------------------------------------------------------------------

def _capture_graph_image(panel) -> str:
    """Grab the graph window as a base64-encoded PNG string."""
    if panel.graph_window is None:
        return ""
    try:
        pixmap = panel.graph_window.grab()
        buf = QBuffer()
        buf.open(QIODevice.OpenModeFlag.WriteOnly)
        pixmap.save(buf, "PNG")
        img_bytes = bytes(buf.data())
        buf.close()
        return base64.b64encode(img_bytes).decode("ascii")
    except Exception:
        return ""


def launch_unified_fit_advisor(panel) -> None:
    """Collect current Unified Fit state and open the AI advisor dialog.

    Called by the "Ask AI advisor" button in UnifiedFitPanel.
    *panel* is the UnifiedFitPanel instance.
    """
    # --- collect state ---
    num_levels = panel.num_levels_spin.value()
    levels = [panel.level_widgets[i].get_parameters() for i in range(num_levels)]

    try:
        background = float(panel.background_value.text())
    except (ValueError, AttributeError):
        background = 0.0

    fit_result = getattr(panel, "fit_result", None)
    chi_sq = fit_result.get("reduced_chi_squared") if fit_result else None
    file_path = (panel.data.get("filepath") if panel.data else None)

    # --- image ---
    img_b64 = _capture_graph_image(panel)

    # --- LLM config from state_manager ---
    cfg      = panel.state_manager.get("ai_advisor") or {}
    provider = cfg.get("provider", "anthropic")
    prov_cfg = cfg.get(provider, {})
    model    = prov_cfg.get("model", "claude-opus-4-7")
    # "url" covers base_url (anthropic/openai) and endpoint (local)
    endpoint = prov_cfg.get("base_url", prov_cfg.get("endpoint", ""))
    user_inst = cfg.get("user_instructions", "")
    proj_ctx  = cfg.get("project_context", "")
    api_key   = _get_api_key(provider)

    if not api_key and provider in ("anthropic", "openai"):
        try:
            from PySide6.QtWidgets import QMessageBox
        except ImportError:
            try:
                from PyQt6.QtWidgets import QMessageBox
            except ImportError:
                from PyQt5.QtWidgets import QMessageBox
        QMessageBox.warning(
            panel, "AI Advisor",
            f"No API key found for {provider}.\n\n"
            "Click 'Configure AI' to enter your key, or set the "
            "ANTHROPIC_API_KEY / OPENAI_API_KEY environment variable.",
        )
        return

    # --- build prompts (includes skills file for this tool) ---
    system_prompt = _build_system_prompt("unified_fit", user_inst, proj_ctx)
    param_text    = _format_param_table(levels, background, num_levels,
                                        chi_sq, file_path)

    # --- launch ---
    thread = AiAdvisorThread(
        provider=provider,
        api_key=api_key,
        model=model,
        endpoint=endpoint,
        system_prompt=system_prompt,
        param_text=param_text,
        img_b64=img_b64,
        parent=panel,
    )

    result_panel = AiAdvisorResultPanel(parent=panel)
    result_panel.set_thread(thread)
    result_panel.show()
