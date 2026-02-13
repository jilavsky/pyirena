#!/usr/bin/env python
"""
Create synthetic test data for pyIrena GUI testing.

This script generates sample HDF5 files in NXcanSAS format.
"""

import numpy as np
import h5py
import os


def create_nxcansas_file(filename, q, intensity, error=None, dq=None):
    """
    Create a simplified NXcanSAS HDF5 file.

    Args:
        filename: Output file path
        q: Q values (1/Å)
        intensity: Intensity values (cm⁻¹)
        error: Error values (optional)
        dq: dQ values (optional)
    """
    with h5py.File(filename, 'w') as f:
        # Create NXcanSAS structure
        entry = f.create_group('entry1')
        entry.attrs['NX_class'] = 'NXentry'

        data_group = entry.create_group('data1')
        data_group.attrs['NX_class'] = 'NXdata'

        # Store datasets
        data_group.create_dataset('Q', data=q)
        data_group.create_dataset('I', data=intensity)

        if error is not None:
            data_group.create_dataset('Idev', data=error)

        if dq is not None:
            data_group.create_dataset('Qdev', data=dq)

        # Set attributes
        data_group['Q'].attrs['units'] = '1/angstrom'
        data_group['I'].attrs['units'] = '1/cm'

        print(f"Created: {filename}")


def main():
    """Generate test data files."""

    # Create testData directory
    os.makedirs('testData', exist_ok=True)

    # Test 1: Simple Porod scattering (sharp interface)
    print("\nGenerating test data files...")
    print("=" * 50)

    q = np.logspace(-3, 0, 150)  # 0.001 to 1.0 Å⁻¹

    # File 1: Porod scattering (spherical particles)
    I1 = 1000 * np.exp(-(q * 50)**2 / 3) + 1e-3 * q**(-4) + 0.01
    error1 = 0.05 * I1
    create_nxcansas_file('testData/spheres_Rg50.h5', q, I1, error1)

    # File 2: Larger particles
    I2 = 5000 * np.exp(-(q * 100)**2 / 3) + 1e-4 * q**(-4) + 0.01
    error2 = 0.05 * I2
    create_nxcansas_file('testData/spheres_Rg100.h5', q, I2, error2)

    # File 3: Mass fractal
    I3 = 500 * np.exp(-(q * 20)**2 / 3) + 0.1 * q**(-2.5) + 0.01
    error3 = 0.05 * I3
    create_nxcansas_file('testData/fractal_aggregate.h5', q, I3, error3)

    # File 4: Multi-level structure
    I4_level1 = 100 * np.exp(-(q * 15)**2 / 3) + 1e-2 * q**(-4)
    I4_level2 = 2000 * np.exp(-(q * 150)**2 / 3) + 1e-3 * q**(-3)
    I4 = I4_level1 + I4_level2 + 0.01
    error4 = 0.05 * I4
    create_nxcansas_file('testData/hierarchical_structure.h5', q, I4, error4)

    # File 5: Surface fractal
    I5 = 800 * np.exp(-(q * 60)**2 / 3) + 5e-3 * q**(-3.5) + 0.01
    error5 = 0.05 * I5
    create_nxcansas_file('testData/surface_fractal.h5', q, I5, error5)

    # Create a complex unified fit example (mentioned in requirements)
    I_complex = (
        # Level 1: small particles
        50 * np.exp(-(q * 10)**2 / 3) + 1e-2 * q**(-4) +
        # Level 2: aggregates
        500 * np.exp(-(q * 80)**2 / 3) + 1e-4 * q**(-2.8) +
        # Level 3: large clusters
        5000 * np.exp(-(q * 200)**2 / 3) +
        # Background
        0.01
    )
    error_complex = 0.03 * I_complex
    create_nxcansas_file('testData/complexUnified.h5', q, I_complex, error_complex)

    print("=" * 50)
    print(f"\nSuccess! Created 6 test files in testData/")
    print("\nYou can now:")
    print("1. Run the GUI: pyirena-gui")
    print("2. Select the testData folder")
    print("3. Select and plot the files")
    print("\nOr test directly:")
    print("  python -m pyirena.gui.launch")


if __name__ == "__main__":
    main()
