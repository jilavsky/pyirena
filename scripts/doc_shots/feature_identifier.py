"""Screenshots for docs/tools/feature-identifier.md.

Dataset: validationData/unified/unified_two_level.h5 loaded in Unified Fit
(two levels, Rg ≈ 120 Å and 1200 Å), then Identify Features… → Detect
segments with the default parameters.
"""

from __future__ import annotations

from ._common import Context, draw_callouts, shot, widget_rect
from .unified_fit import _open as _open_unified

DIALOG_SIZE = (440, 760)  # default 440×600; taller so the advanced parameters fit


def _detected(ctx: Context, advanced: bool):
    panel = _open_unified(ctx)
    panel._open_feature_identifier()
    dlg = panel._feature_dialog
    dlg.show_advanced.setChecked(advanced)
    ctx.show(dlg, DIALOG_SIZE)
    dlg._on_detect()
    ctx.settle()
    return panel, dlg


@shot("feature-identifier-01-window")
def window(ctx: Context):
    _, dlg = _detected(ctx, advanced=True)
    r = lambda w: widget_rect(w, dlg)  # noqa: E731
    marks = [
        (1, r(dlg.show_advanced), "r", (14, 0)),
        (2, r(dlg.detect_btn), "tl", (2, 2)),
        (3, r(dlg.clear_btn), "tr", (-2, 2)),
        (4, r(dlg.summary), "tr", (-16, 16)),
        (5, r(dlg._advanced_box), "tr", (-16, 10)),
    ]
    return draw_callouts(dlg.grab(), marks)


@shot("feature-identifier-02-markers")
def markers(ctx: Context):
    panel, _ = _detected(ctx, advanced=False)
    return panel.graph_window.grab()
