#!/usr/bin/env python3
"""
Test script for the Unified Fit GUI.

This script creates test data and launches the Unified Fit GUI to verify
that the Graph Unified button works correctly.
"""

import numpy as np
import sys

# Test if imports work
try:
    from pyirena.gui.unified_fit import UnifiedFitPanel
    from pyirena.core.unified import UnifiedFitModel, UnifiedLevel
    print("✓ Imports successful")
except ImportError as e:
    print(f"✗ Import error: {e}")
    sys.exit(1)

try:
    from PySide6.QtWidgets import QApplication
except ImportError:
    try:
        from PyQt6.QtWidgets import QApplication
    except ImportError:
        print("✗ Neither PySide6 nor PyQt6 found")
        sys.exit(1)

print("✓ Qt imports successful")

# Create test data
print("\nCreating test data...")
q = np.logspace(-3, 0, 100)

# Create a simple 2-level unified fit
# Level 1: Large particles (Rg=100, G=1e10, P=4, B=1e6)
level1 = UnifiedLevel(Rg=100, G=1e10, P=4.0, B=1e6)

# Level 2: Small particles (Rg=10, G=1e8, P=3.5, B=1e4)
level2 = UnifiedLevel(Rg=10, G=1e8, P=3.5, B=1e4)

# Create model and calculate
model = UnifiedFitModel()
model.levels = [level1, level2]
model.background = 1e5

intensity = model.calculate_intensity(q)

# Add some noise
np.random.seed(42)
noise = 0.03 * intensity
intensity_noisy = intensity + np.random.normal(0, noise)
error = noise

print(f"✓ Created test data: {len(q)} points")
print(f"  Q range: {q.min():.2e} to {q.max():.2e}")
print(f"  Intensity range: {intensity_noisy.min():.2e} to {intensity_noisy.max():.2e}")

# Launch GUI
print("\nLaunching Unified Fit GUI...")
app = QApplication(sys.argv)
app.setStyle('Fusion')

window = UnifiedFitPanel()
window.set_data(q, intensity_noisy, error, "Test Data (2 levels)")

# Set some initial parameters
window.num_levels_spin.setValue(2)

# Level 1 initial guess (slightly off from true values)
window.level_widgets[0].g_value.setText("5e9")
window.level_widgets[0].rg_value.setText("120")
window.level_widgets[0].b_value.setText("5e5")
window.level_widgets[0].p_value.setText("4.0")

# Level 2 initial guess
window.level_widgets[1].g_value.setText("5e7")
window.level_widgets[1].rg_value.setText("12")
window.level_widgets[1].b_value.setText("5e3")
window.level_widgets[1].p_value.setText("3.5")

window.background_value.setText("5e4")

window.show()

print("\n" + "="*60)
print("GUI launched successfully!")
print("="*60)
print("\nTest Instructions:")
print("1. Click 'Graph Unified' button (green) to calculate model")
print("2. Verify that a fit curve appears in the top graph")
print("3. Verify that residuals appear in the bottom graph")
print("4. Try adjusting parameters and clicking 'Graph Unified' again")
print("5. Try clicking 'Fit' button to run optimization")
print("\nTrue parameters:")
print("  Level 1: G=1e10, Rg=100, B=1e6, P=4.0")
print("  Level 2: G=1e8, Rg=10, B=1e4, P=3.5")
print("  Background: 1e5")
print("="*60)

sys.exit(app.exec())
