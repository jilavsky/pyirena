"""Tests for pyirena.api.plotting (headless matplotlib)."""
from __future__ import annotations

from pathlib import Path

import pytest


pytest.importorskip("matplotlib")


def test_plot_iq_returns_png(synth_nxcansas_file, tmp_path):
    from pyirena.api import plot_iq
    out = tmp_path / "iq.png"
    result = plot_iq([str(synth_nxcansas_file)], output_path=str(out),
                      return_base64=False)
    assert Path(result["path"]).exists()
    assert result["format"] == "png"
    assert result["n_files"] == 1
    assert Path(result["path"]).stat().st_size > 0


def test_plot_iq_base64(synth_nxcansas_file, tmp_path):
    from pyirena.api import plot_iq
    result = plot_iq([str(synth_nxcansas_file)],
                      output_path=str(tmp_path / "x.png"),
                      return_base64=True)
    assert result["base64_png"]
    # Quick sanity: base64 PNG starts with iVBOR (PNG header)
    assert result["base64_png"].startswith("iVBOR")


def test_plot_parameter_trend(synth_folder_multi, tmp_path):
    from pyirena.api import plot_parameter_trend
    out = tmp_path / "trend.png"
    result = plot_parameter_trend(
        folder=str(synth_folder_multi),
        tool="unified_fit",
        parameter="Rg",
        subgroup_index=1,
        output_path=str(out),
        return_base64=False,
    )
    assert Path(result["path"]).exists()
    assert result["n_files"] == 3
