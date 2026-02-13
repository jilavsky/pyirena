#!/usr/bin/env python
"""
Convert HDF5 files to .dat text files.

This script reads all .h5 files in testData/ and creates corresponding .dat files
with Q, Intensity, and Error columns.
"""

import os
import h5py
import numpy as np
from pathlib import Path


def read_h5_file(filepath):
    """
    Read Q, I, Error from HDF5 file (supports both NXcanSAS and simple structure).

    Args:
        filepath: Path to HDF5 file

    Returns:
        tuple: (Q, Intensity, Error) arrays, or (None, None, None) if failed
    """
    try:
        with h5py.File(filepath, 'r') as f:
            # Try common simple structures
            possible_paths = [
                ('entry1/data1', 'Q', 'I'),
                ('entry/data', 'Q', 'I'),
                ('data', 'Q', 'I'),
                ('', 'Q', 'I'),  # Root level
                ('', 'Q', 'Intensity'),
            ]

            for base, q_name, i_name in possible_paths:
                q_path = f"{base}/{q_name}".strip('/') if base else q_name
                i_path = f"{base}/{i_name}".strip('/') if base else i_name

                if q_path in f and i_path in f:
                    Q = f[q_path][()]
                    I = f[i_path][()]

                    # Try to find error data
                    error = None
                    for err_name in ['Idev', 'Error', 'error', 'I_error']:
                        err_path = f"{base}/{err_name}".strip('/') if base else err_name
                        if err_path in f:
                            error = f[err_path][()]
                            break

                    # If no error found, create 5% error
                    if error is None:
                        error = 0.05 * np.abs(I)

                    print(f"  Read data from path: {base if base else 'root'}")
                    return Q, I, error

            print(f"  Warning: No recognizable data structure found")
            return None, None, None

    except Exception as e:
        print(f"  Error reading file: {e}")
        return None, None, None


def write_dat_file(filepath, Q, I, error):
    """
    Write Q, I, Error to .dat text file.

    Args:
        filepath: Output file path
        Q: Q values
        I: Intensity values
        error: Error values
    """
    with open(filepath, 'w') as f:
        # Write header
        f.write("# Small-Angle Scattering Data\n")
        f.write("# Converted from HDF5 format\n")
        f.write("#\n")
        f.write("# Column 1: Q (Å⁻¹)\n")
        f.write("# Column 2: Intensity (cm⁻¹)\n")
        f.write("# Column 3: Error (cm⁻¹)\n")
        f.write("#\n")
        f.write(f"# {'Q':>14s} {'Intensity':>15s} {'Error':>15s}\n")

        # Write data
        for q, intensity, err in zip(Q, I, error):
            f.write(f"{q:15.6e} {intensity:15.6e} {err:15.6e}\n")

    print(f"  Wrote {len(Q)} data points")


def main():
    """Convert all HDF5 files in testData to .dat format."""
    testdata_dir = Path('testData')

    if not testdata_dir.exists():
        print("Error: testData directory not found!")
        return

    print("\nConverting HDF5 files to .dat format...")
    print("=" * 60)

    # Find all .h5 files
    h5_files = list(testdata_dir.glob('*.h5'))

    if not h5_files:
        print("No .h5 files found in testData/")
        return

    converted_count = 0

    for h5_file in sorted(h5_files):
        print(f"\nProcessing: {h5_file.name}")

        # Read HDF5 file
        Q, I, error = read_h5_file(h5_file)

        if Q is None:
            print(f"  Skipping (could not read data)")
            continue

        # Create .dat filename
        dat_file = testdata_dir / (h5_file.stem + '.dat')

        # Write .dat file
        write_dat_file(dat_file, Q, I, error)
        print(f"  Created: {dat_file.name}")

        converted_count += 1

    print("=" * 60)
    print(f"\nSuccess! Converted {converted_count} files to .dat format")
    print("\nYou can now:")
    print("1. Run the GUI: pyirena-gui")
    print("2. Select the testData folder")
    print("3. Switch file type to 'Text Files (.txt, .dat)'")
    print("4. See both .h5 and .dat files")


if __name__ == "__main__":
    main()
