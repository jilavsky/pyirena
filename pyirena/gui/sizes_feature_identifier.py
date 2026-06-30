"""Feature Identifier dialog for the Size Distribution panel.

Subclasses the Unified Fit :class:`FeatureIdentifierDialog` so it shares the
exact same log-log slope segmentation, marker overlays, and clear/remove
behaviour.  On top of that it appends the **size-distribution recommendation**
— the very output the AI control tool ``suggest_sizes_setup`` returns — so a
user can see (and sanity-check) what the agent would be told: suggested radius
grid, inversion Q-range, power-law / flat-background windows, and a suitability
verdict with warnings.

Both the GUI and the AI tool call the same core function
(:func:`pyirena.core.sizes.recommend_sizes_setup`), so the displayed
recommendation matches what the agent receives (with default segmentation
parameters).
"""
from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from pyirena.gui.feature_identifier import FeatureIdentifierDialog

if TYPE_CHECKING:
    from pyirena.gui.sizes_panel import SizesFitPanel


class SizesFeatureIdentifierDialog(FeatureIdentifierDialog):
    """Non-modal feature identifier tailored for the Size Distribution panel."""

    def __init__(self, parent: "SizesFitPanel"):
        super().__init__(parent)

    # ------------------------------------------------------------------
    # Presentation overrides
    # ------------------------------------------------------------------

    def _window_title(self) -> str:
        return "Feature Identifier — Size-distribution setup"

    def _help_text(self) -> str:
        return (
            "<b>Visualisation only</b> — never modifies your fit. The coloured "
            "regions are the same power-law slope segments as in Unified Fit. "
            "Below the segment list is the <b>size-distribution recommendation</b> "
            "(radius grid, inversion Q-range, background windows) and a "
            "suitability verdict — the same output the AI assistant's "
            "<i>suggest_sizes_setup</i> tool produces. Use it to judge whether a "
            "size distribution is appropriate for this data and where to place "
            "the cursors / background windows."
        )

    def _extra_summary_lines(self, result) -> list:
        """Append the size-distribution recommendation to the summary."""
        from pyirena.core.sizes import recommend_sizes_setup  # noqa: PLC0415

        data = getattr(self._panel, "data", None)
        if not data or data.get("Q") is None or data.get("Intensity") is None:
            return []

        q = np.asarray(data["Q"], dtype=float)
        I = np.asarray(data["Intensity"], dtype=float)
        err = data.get("Error")
        err = np.asarray(err, dtype=float) if err is not None else None

        try:
            # Use the dialog's current segmentation config so the recommendation
            # is consistent with the markers shown above.
            rec = recommend_sizes_setup(q, I, sigma=err, config=self._current_cfg())
        except Exception as exc:  # pragma: no cover - defensive
            return ["", f"Size recommendation failed: {exc}"]

        r = rec["recommended"]
        lines: list[str] = ["", "─" * 48, "SIZE-DISTRIBUTION RECOMMENDATION"]

        verdict = "✓ suitable" if rec["suitable"] else "✗ NOT a clean candidate"
        lines.append(f"  Suitability: {verdict}")
        lines.append("")
        lines.append(f"  Radius grid:     r ∈ [{_fmt(r['r_min'])}, "
                     f"{_fmt(r['r_max'])}] Å")
        lines.append(f"  Inversion Q:     [{_fmt(r['inversion_q_min'])}, "
                     f"{_fmt(r['inversion_q_max'])}] Å⁻¹")
        if r["power_law_q_min"] is not None:
            lines.append(f"  Power-law win.:  [{_fmt(r['power_law_q_min'])}, "
                         f"{_fmt(r['power_law_q_max'])}] Å⁻¹ (low-Q)")
        else:
            lines.append("  Power-law win.:  none detected")
        if r["background_q_min"] is not None:
            lines.append(f"  Flat-bg window:  [{_fmt(r['background_q_min'])}, "
                         f"{_fmt(r['background_q_max'])}] Å⁻¹ (high-Q)")
        else:
            lines.append("  Flat-bg window:  none detected")

        if rec["warnings"]:
            lines.append("")
            lines.append("  Warnings:")
            for w in rec["warnings"]:
                lines.append(f"    • {w}")
        return lines


def _fmt(v) -> str:
    """Compact number formatting that tolerates None/NaN."""
    if v is None:
        return "—"
    try:
        f = float(v)
    except (TypeError, ValueError):
        return str(v)
    if not np.isfinite(f):
        return "—"
    return f"{f:.4g}"


__all__ = ["SizesFeatureIdentifierDialog"]
