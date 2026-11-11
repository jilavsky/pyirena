"""Euler-Poincaré characteristic in core/morphology.py.

The implementation replaced an optional ``skimage.measure.euler_number`` call
whose absence was reported as χ = 0 — a real topology (a torus). These cases
have analytically known answers under 6-connectivity (skimage
``connectivity=1``); the numpy version was also cross-checked against skimage
on 300 random volumes with zero mismatches when it was written (issue #38).
"""

import numpy as np
import pytest

from pyirena.core.morphology import _euler_number_6conn, compute_morphology_metrics


def _empty(n=7):
    return np.zeros((n, n, n), dtype=bool)


def test_empty_is_zero():
    assert _euler_number_6conn(_empty()) == 0


def test_single_voxel_and_solid_block_are_one():
    m = _empty()
    m[3, 3, 3] = True
    assert _euler_number_6conn(m) == 1
    m = _empty()
    m[1:6, 2:5, 1:4] = True
    assert _euler_number_6conn(m) == 1


def test_block_touching_array_edge_is_one():
    """Outside the array counts as background, so a full array is one blob."""
    assert _euler_number_6conn(np.ones((4, 5, 6), dtype=bool)) == 1


def test_shell_with_cavity_is_two():
    m = _empty()
    m[1:6, 1:6, 1:6] = True
    m[2:5, 2:5, 2:5] = False
    assert _euler_number_6conn(m) == 2


def test_ring_is_zero():
    m = _empty()
    m[1:6, 1:6, 3] = True
    m[2:5, 2:5, 3] = False
    assert _euler_number_6conn(m) == 0


def test_two_rings_is_minus_two_handles_plus_one_component():
    """A thick plate with two through-holes: χ = 1 − 2 = −1."""
    m = np.zeros((11, 7, 3), dtype=bool)
    m[1:10, 1:6, :] = True
    m[3, 3, :] = False
    m[7, 3, :] = False
    assert _euler_number_6conn(m) == -1


def test_diagonal_voxels_are_separate_components():
    """6-connectivity: voxels sharing only an edge or a corner do not connect."""
    m = _empty()
    m[2, 2, 2] = True
    m[3, 3, 2] = True          # edge-neighbour
    m[5, 5, 5] = True
    m[4, 4, 4] = True          # corner-neighbour
    assert _euler_number_6conn(m) == 4


def test_rejects_non_3d():
    with pytest.raises(ValueError):
        _euler_number_6conn(np.zeros((4, 4), dtype=bool))


def test_metrics_report_euler_for_ring_structure():
    """End to end: a ring minority phase must not come back as a default."""
    vox = np.zeros((9, 9, 9), dtype=np.uint8)
    vox[1:8, 1:8, 4] = 1
    vox[2:7, 2:7, 4] = 0
    vox[1:8, 1:8, 1] = 1        # plus one solid plate → χ = 0 + 1
    metrics = compute_morphology_metrics(vox, voxel_pitch_A=10.0)
    assert metrics.euler_number == 1
