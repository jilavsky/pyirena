# HDF5 NXcanSAS Structure Reference

This document provides instructions for agents to locate reduced Intensity ($I$), $Q$, and uncertainties ($dI$) within a pyirena-generated HDF5 file following the NXcanSAS standard.

## 1. Locating Primary Data (Reduced Intensity)

The primary scattering data is stored within the `entry/` group. To find the correct sample data, you must first identify the default sample name.

### Step 1: Identify Sample Name
Read the `default` attribute of the `entry/` group.
- **Path:** `entry/`
- **Attribute:** `@default`
- **Value:** This will be the name of the primary (desmeared) sample.

### Step 2: Access the `sasdata` Group
Once you have the sample name (e.g., `<sample_name>`), navigate to its `sasdata` group.
- **Path:** `entry/<sample_name>/sasdata/`

### Step 3: Read Datasets
The following datasets contain the reduced scattering data:

| Dataset | Description | Units | Attributes to Note |
|---------|-------------|--------|---------------------|
| `Q`     | Scattering vector | $1/\text{\AA}$ | `@long_name`: "Q (A^-1)", `@resolutions`: "Qdev" |
| `I`     | Reduced Intensity | $1/\text{cm}$ | `@long_name`: "Intensity", `@uncertainties`: "Idev" |
| `Idev`  | Uncertainties ($dI$) | $1/\text{cm}$ or $\text{cm}^2/\text{cm}^3$ | `@long_name`: "Uncertainties" |

**Note on $dI$:** If the `Idev` dataset is missing from the `sasdata` group, it means uncertainties were not recorded or are unavailable for that specific dataset.

## 2. Handling USAXS (SMR vs. Desmeared)
For USAXS files, two variants of the sample group may exist:
1. **Desmeared (Primary):** `entry/<sample_name>/`
2. **Slit-Smeared (SMR):** `entry/<sample_name>_SMR/`

**Rule:** Always prefer the desmeared variant unless specifically instructed otherwise. The `entry/@default` attribute automatically points to the desmeared variant.

## 3. Summary of Paths for Agents
To extract the core data, use these relative paths from the root:

- **Sample Name:** `entry/@default`
- **Q Array:** `entry/<sample_name>/sasdata/Q`
- **I Array:** `entry/<sample_name>/sasdata/I`
- **dI Array:** `entry/<sample_name>/sasdata/Idev`

## 4. Example Python Snippet (h5py)
```python
import h5py

with h5py.File("file.h5", "r") as f:
    # 1. Get the default sample name
    sample_name = f["entry"].attrs["default"]
    
    # 2. Access the sasdata group
    sd = f[f"entry/{sample_name}/sasdata"]
    
    # 3. Extract data
    Q  = sd["Q"][:]
    I  = sd["I"][:]
    dI = sd["Idev"][:] if "Idev" in sd else None
```
