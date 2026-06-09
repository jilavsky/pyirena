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
import os
import traceback
from pathlib import Path
from typing import Optional

try:
    from PySide6.QtCore import Qt, QThread, QTimer, Signal
    from PySide6.QtGui import QFont
    from PySide6.QtWidgets import (
        QDialog, QDialogButtonBox, QFormLayout, QGroupBox, QHBoxLayout,
        QLabel, QLineEdit, QPlainTextEdit, QPushButton, QRadioButton,
        QSizePolicy, QTextEdit, QVBoxLayout, QWidget,
    )
    from PySide6.QtCore import QBuffer, QByteArray, QIODevice
except ImportError:
    try:
        from PyQt6.QtCore import Qt, QThread, QTimer, pyqtSignal as Signal
        from PyQt6.QtGui import QFont
        from PyQt6.QtWidgets import (
            QDialog, QDialogButtonBox, QFormLayout, QGroupBox, QHBoxLayout,
            QLabel, QLineEdit, QPlainTextEdit, QPushButton, QRadioButton,
            QSizePolicy, QTextEdit, QVBoxLayout, QWidget,
        )
        from PyQt6.QtCore import QBuffer, QByteArray, QIODevice
    except ImportError:
        from PyQt5.QtCore import Qt, QThread, QTimer, pyqtSignal as Signal
        from PyQt5.QtGui import QFont
        from PyQt5.QtWidgets import (
            QDialog, QDialogButtonBox, QFormLayout, QGroupBox, QHBoxLayout,
            QLabel, QLineEdit, QPlainTextEdit, QPushButton, QRadioButton,
            QSizePolicy, QTextEdit, QVBoxLayout, QWidget,
        )
        from PyQt5.QtCore import QBuffer, QByteArray, QIODevice


# ---------------------------------------------------------------------------
# API key helpers (OS keyring ↔ env-var fallback)
# ---------------------------------------------------------------------------

_KEYRING_SERVICE = "pyirena-ai"
_KEYRING_KEYS = {
    "anthropic": "anthropic_api_key",
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
    return os.environ.get("OPENAI_API_KEY", "")


def _set_api_key(provider: str, key: str) -> None:
    """Store API key in OS keyring."""
    try:
        import keyring
        keyring.set_password(_KEYRING_SERVICE, _KEYRING_KEYS.get(provider, provider), key)
    except Exception:
        pass  # silently fall back to env var usage


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


def _build_system_prompt(user_instructions: str, project_context: str) -> str:
    parts = [_TOOL_SYSTEM_PROMPT]
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
        lines.append("| Parameter | Value | Fixed | Lo | Hi |")
        lines.append("|-----------|-------|-------|----|-----|")
        for param in ("Rg", "G", "B", "P", "ETA", "PACK", "RgCutoff"):
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
        if lv.get("correlated"):
            extras.append("correlations ON")
        if lv.get("estimate_B"):
            extras.append("estimate_B ON")
        if lv.get("link_rgco"):
            extras.append("link_RgCO ON")
        if extras:
            lines.append(f"Flags: {', '.join(extras)}")
        lines.append("")

    lines.append("### Background")
    lines.append(f"Value: {background:.4g}  |  Fixed: No")
    return "\n".join(lines)


class AiAdvisorThread(QThread):
    """Background thread that calls the LLM and emits the response."""

    result_ready   = Signal(str)
    error_occurred = Signal(str)

    def __init__(
        self,
        provider:          str,
        api_key:           str,
        model:             str,
        endpoint:          str,
        system_prompt:     str,
        param_text:        str,
        img_b64:           str,
        anthropic_base_url: str = "",
        parent=None,
    ):
        super().__init__(parent)
        self.provider           = provider
        self.api_key            = api_key
        self.model              = model
        self.endpoint           = endpoint.rstrip("/")
        self.system_prompt      = system_prompt
        self.param_text         = param_text
        self.img_b64            = img_b64
        self.anthropic_base_url = anthropic_base_url.strip()

    def run(self):
        try:
            if self.provider == "anthropic":
                response_text = self._call_anthropic()
            else:
                response_text = self._call_local()
            self.result_ready.emit(response_text)
        except Exception as exc:
            self.error_occurred.emit(f"{type(exc).__name__}: {exc}\n\n{traceback.format_exc()}")

    def _call_anthropic(self) -> str:
        import anthropic  # noqa: PLC0415 — optional dep, imported on use
        kwargs = {"api_key": self.api_key}
        if self.anthropic_base_url:
            kwargs["base_url"] = self.anthropic_base_url
        client = anthropic.Anthropic(**kwargs)
        content = []
        if self.img_b64:
            content.append({
                "type": "image",
                "source": {
                    "type":       "base64",
                    "media_type": "image/png",
                    "data":       self.img_b64,
                },
            })
        content.append({"type": "text", "text": self.param_text})
        msg = client.messages.create(
            model=self.model,
            max_tokens=1024,
            system=self.system_prompt,
            messages=[{"role": "user", "content": content}],
        )
        return msg.content[0].text

    def _call_local(self) -> str:
        """Call an OpenAI-compatible local endpoint (LM Studio / Ollama)."""
        import httpx  # noqa: PLC0415
        messages = [{"role": "system", "content": self.system_prompt}]
        user_parts = []
        if self.img_b64:
            user_parts.append({
                "type":      "image_url",
                "image_url": {"url": f"data:image/png;base64,{self.img_b64}"},
            })
        user_parts.append({"type": "text", "text": self.param_text})
        messages.append({"role": "user", "content": user_parts})

        resp = httpx.post(
            f"{self.endpoint}/chat/completions",
            headers={"Authorization": f"Bearer {self.api_key or 'local'}"},
            json={"model": self.model, "messages": messages, "max_tokens": 1024},
            timeout=120,
        )
        resp.raise_for_status()
        return resp.json()["choices"][0]["message"]["content"]


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

        # Response area
        self._text_edit = QTextEdit()
        self._text_edit.setReadOnly(True)
        mono = QFont("Courier New", 10)
        self._text_edit.setFont(mono)
        self._text_edit.setPlaceholderText("Response will appear here…")
        layout.addWidget(self._text_edit)

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
        thread.result_ready.connect(self._on_result)
        thread.error_occurred.connect(self._on_error)
        thread.start()

    def _on_result(self, text: str):
        self._timer.stop()
        self._status_label.setText("AI advisor response:")
        self._text_edit.setPlainText(text)

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
    """Settings dialog: provider, model, API key, instructions."""

    def __init__(self, state_manager, parent=None):
        super().__init__(parent)
        self.state_manager = state_manager
        self.setWindowTitle("AI Advisor — Configure")
        self.setModal(True)
        self.setMinimumWidth(480)
        self._build_ui()
        self._load_state()

    def _build_ui(self):
        layout = QVBoxLayout(self)

        # --- Provider & model group ---
        prov_group = QGroupBox("Provider and model")
        prov_layout = QVBoxLayout(prov_group)

        radio_row = QHBoxLayout()
        self._radio_anthropic = QRadioButton("Anthropic (Claude)")
        self._radio_local     = QRadioButton("Local (LM Studio / Ollama)")
        radio_row.addWidget(self._radio_anthropic)
        radio_row.addWidget(self._radio_local)
        prov_layout.addLayout(radio_row)

        form = QFormLayout()

        self._model_edit = QLineEdit()
        self._model_edit.setPlaceholderText("e.g. claude-opus-4-7")
        form.addRow("Model name:", self._model_edit)

        self._base_url_edit = QLineEdit()
        self._base_url_edit.setPlaceholderText(
            "Leave blank for default (https://api.anthropic.com)"
        )
        self._base_url_label = QLabel("Custom API URL:")
        form.addRow(self._base_url_label, self._base_url_edit)

        self._endpoint_edit = QLineEdit()
        self._endpoint_edit.setPlaceholderText("http://localhost:1234/v1")
        self._endpoint_label = QLabel("Endpoint URL:")
        form.addRow(self._endpoint_label, self._endpoint_edit)

        self._key_edit = QLineEdit()
        self._key_edit.setEchoMode(QLineEdit.EchoMode.Password)
        self._key_edit.setPlaceholderText("Leave blank to use env var")
        form.addRow("API key:", self._key_edit)

        prov_layout.addLayout(form)

        self._test_btn = QPushButton("Test connection")
        self._test_btn.clicked.connect(self._test_connection)
        self._test_label = QLabel("")
        test_row = QHBoxLayout()
        test_row.addWidget(self._test_btn)
        test_row.addWidget(self._test_label)
        test_row.addStretch()
        prov_layout.addLayout(test_row)

        layout.addWidget(prov_group)

        # Update endpoint visibility on radio change
        self._radio_anthropic.toggled.connect(self._update_endpoint_visibility)
        self._radio_local.toggled.connect(self._update_endpoint_visibility)

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

        # --- Dialog buttons ---
        btn_box = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel
        )
        btn_box.accepted.connect(self._save_and_accept)
        btn_box.rejected.connect(self.reject)
        layout.addWidget(btn_box)

    def _update_endpoint_visibility(self):
        local = self._radio_local.isChecked()
        # Local endpoint — shown for local provider only
        self._endpoint_label.setVisible(local)
        self._endpoint_edit.setVisible(local)
        # Custom Anthropic base URL — shown for Anthropic only
        self._base_url_label.setVisible(not local)
        self._base_url_edit.setVisible(not local)

    def _load_state(self):
        cfg = self.state_manager.get("ai_advisor") or {}
        provider = cfg.get("provider", "anthropic")
        if provider == "local":
            self._radio_local.setChecked(True)
        else:
            self._radio_anthropic.setChecked(True)
        self._update_endpoint_visibility()

        self._model_edit.setText(cfg.get("model", "claude-opus-4-7"))
        self._base_url_edit.setText(cfg.get("anthropic_base_url", ""))
        self._endpoint_edit.setText(cfg.get("local_endpoint", "http://localhost:1234/v1"))
        self._key_edit.setText(_get_api_key(provider))
        self._general_edit.setPlainText(cfg.get("user_instructions", ""))
        self._context_edit.setPlainText(cfg.get("project_context", ""))

    def _save_and_accept(self):
        provider = "local" if self._radio_local.isChecked() else "anthropic"
        api_key  = self._key_edit.text().strip()
        if api_key:
            _set_api_key(provider, api_key)

        self.state_manager.set("ai_advisor", "provider",             provider)
        self.state_manager.set("ai_advisor", "model",              self._model_edit.text().strip())
        self.state_manager.set("ai_advisor", "anthropic_base_url", self._base_url_edit.text().strip())
        self.state_manager.set("ai_advisor", "local_endpoint",     self._endpoint_edit.text().strip())
        self.state_manager.set("ai_advisor", "user_instructions", self._general_edit.toPlainText())
        self.state_manager.set("ai_advisor", "project_context",   self._context_edit.toPlainText())
        self.state_manager.save()
        self.accept()

    def _test_connection(self):
        self._test_label.setText("Testing…")
        provider = "local" if self._radio_local.isChecked() else "anthropic"
        api_key  = self._key_edit.text().strip() or _get_api_key(provider)
        model    = self._model_edit.text().strip()
        endpoint = self._endpoint_edit.text().strip()

        class _TestThread(QThread):
            done = Signal(str)

            def __init__(self, provider, api_key, model, endpoint, parent=None):
                super().__init__(parent)
                self._p, self._k, self._m, self._e = provider, api_key, model, endpoint

            def run(self):
                try:
                    if self._p == "anthropic":
                        import anthropic
                        c = anthropic.Anthropic(api_key=self._k)
                        r = c.messages.create(
                            model=self._m, max_tokens=10,
                            messages=[{"role": "user", "content": "Hi"}],
                        )
                        self.done.emit(f"OK — model: {r.model}")
                    else:
                        import httpx
                        r = httpx.post(
                            f"{self._e.rstrip('/')}/chat/completions",
                            headers={"Authorization": f"Bearer {self._k or 'local'}"},
                            json={"model": self._m,
                                  "messages": [{"role": "user", "content": "Hi"}],
                                  "max_tokens": 10},
                            timeout=15,
                        )
                        r.raise_for_status()
                        self.done.emit("OK — local endpoint responded")
                except Exception as exc:
                    self.done.emit(f"Error: {exc}")

        self._test_thread = _TestThread(provider, api_key, model, endpoint, self)
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
    cfg = panel.state_manager.get("ai_advisor") or {}
    provider         = cfg.get("provider", "anthropic")
    model            = cfg.get("model", "claude-opus-4-7")
    endpoint         = cfg.get("local_endpoint", "http://localhost:1234/v1")
    anthropic_base_url = cfg.get("anthropic_base_url", "")
    user_inst        = cfg.get("user_instructions", "")
    proj_ctx         = cfg.get("project_context", "")
    api_key          = _get_api_key(provider)

    if not api_key and provider == "anthropic":
        from PySide6.QtWidgets import QMessageBox
        QMessageBox.warning(
            panel, "AI Advisor",
            "No Anthropic API key found.\n\n"
            "Click 'Configure AI' to enter your key, or set the "
            "ANTHROPIC_API_KEY environment variable.",
        )
        return

    # --- build prompts ---
    system_prompt = _build_system_prompt(user_inst, proj_ctx)
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
        anthropic_base_url=anthropic_base_url,
        parent=panel,
    )

    result_panel = AiAdvisorResultPanel(parent=panel)
    result_panel.set_thread(thread)
    result_panel.show()
