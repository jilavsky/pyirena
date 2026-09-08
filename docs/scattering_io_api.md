# Scattering discovery and loading API

PyIrena 1.1.0 provides a stable, Qt-independent boundary for applications
that need its tested data readers without depending on GUI internals.

```python
from pyirena.io import discover_scattering, load_scattering

locations = discover_scattering("sample.h5")
for location in locations:
    record = load_scattering(location, q_unit="1/A", error_fraction=0.05)
    print(location.internal_path, record.q, record.intensity)
```

`ScatteringLocation` records the source path, internal HDF5 path, display name,
variant, and summary metadata. `ScatteringRecord` returns Q, intensity,
optional intensity uncertainty and dQ, units, metadata, and provenance.

Supported inputs are NXcanSAS/simple HDF5 (`.h5`, `.hdf5`, `.hdf`, `.nxs`) and
text curves (`.dat`, `.txt`, `.csv`). Discovery returns every selectable curve,
including slit-smeared variants. Text loading is read-only: it never creates a
converted HDF5 sibling.

Companion Qt plotting applications can request the small shared dependency set
with:

```bash
python -m pip install "pyirena[qtplot]>=1.1.0"
```

Igor H5XP data writers are also public:

```python
from pyirena.io import create_h5xp, write_iq_data
```
