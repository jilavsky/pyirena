#!/usr/bin/env python3
"""
Test script for PyQtGraph Unified Fit GUI.
"""

import numpy as np
import sys

# Test imports
try:
    from pyirena.gui.unified_fit import UnifiedFitPanel
    from pyirena.core.unified import UnifiedFitModel, UnifiedLevel
    import pyqtgraph as pg
    print("✓ All imports successful (pyqtgraph version)")
except ImportError as e:
    print(f"✗ Import error: {e}")
    sys.exit(1)

try:
    from PySide6.QtWidgets import QApplication
except ImportError:
    try:
        from PyQt6.QtWidgets import QApplication
    except ImportError:
        try:
            from PyQt5.QtWidgets import QApplication
        except ImportError:
            print("✗ No Qt library found (PySide6, PyQt6, or PyQt5)")
            sys.exit(1)

print("✓ Qt imports successful")

# Create test data
print("\nCreating test data...")
q = np.logspace(-3, 0, 100)

# Create a simple 2-level unified fit
level1 = UnifiedLevel(Rg=100, G=1e10, P=4.0, B=1e6)
level2 = UnifiedLevel(Rg=10, G=1e8, P=3.5, B=1e4)

# Create model and calculate
model = UnifiedFitModel()
model.levels = [level1, level2]
model.background = 1e5
model.num_levels = 2

intensity = model.calculate_intensity(q)

# Add noise
np.random.seed(42)
noise = 0.03 * intensity
intensity_noisy = intensity + np.random.normal(0, noise)
error = noise

print(f"✓ Created test data: {len(q)} points")
print(f"  Q range: {q.min():.2e} to {q.max():.2e}")
print(f"  Intensity range: {intensity_noisy.min():.2e} to {intensity_noisy.max():.2e}")

# Launch GUI
print("\nLaunching PyQtGraph Unified Fit GUI...")
app = QApplication(sys.argv)
app.setStyle('Fusion')

window = UnifiedFitPanel()
window.set_data(q, intensity_noisy, error, "Test Data (2 levels)")

# Set initial parameters
window.num_levels_spin.setValue(2)

# Level 1 initial guess
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
print("PyQtGraph GUI launched successfully!")
print("="*60)
print("\nTest Instructions:")
print("1. Check that the plot shows data with error bars")
print("2. Check that RED and BLUE cursor lines appear")
print("3. Try dragging the cursors - should be smooth!")
print("4. Click 'Graph Unified' - should see fit curve instantly")
print("5. Click 'Fit' - should optimize and update parameters")
print("6. Try 'Update automatically?' checkbox")
print("7. Change a parameter - should update in real-time!")
print("\nExpected performance:")
print("  - Cursor dragging: Smooth, no lag")
print("  - Graph updates: Near-instant (<50ms)")
print("  - Auto-update: Real-time response")
print("\nTrue parameters:")
print("  Level 1: G=1e10, Rg=100, B=1e6, P=4.0")
print("  Level 2: G=1e8, Rg=10, B=1e4, P=3.5")
print("  Background: 1e5")
print("="*60)

sys.exit(app.exec())
