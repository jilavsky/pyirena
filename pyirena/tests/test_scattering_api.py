import h5py
import numpy as np

from pyirena.io import (
    ScatteringLocation,
    create_h5xp,
    discover_scattering,
    load_scattering,
    write_iq_data,
)


def test_text_api_is_read_only(tmp_path):
    source = tmp_path / "curve.dat"
    source.write_text("# q I\n1 10\n2 5\n", encoding="utf-8")
    locations = discover_scattering(source)
    assert locations == [ScatteringLocation(source, display_name=source.name, variant="text")]
    record = load_scattering(locations[0], q_unit="1/nm")
    np.testing.assert_allclose(record.q, [0.1, 0.2])
    assert not (tmp_path / "curve.h5").exists()


def test_hdf5_api_discovers_all_sasdata(tmp_path):
    source = tmp_path / "curves.h5"
    with h5py.File(source, "w") as handle:
        for path in ("entry/data", "entry/data_SMR"):
            group = handle.create_group(path)
            group.attrs["canSAS_class"] = "SASdata"
            group.attrs["signal"] = "I"
            group.attrs["I_axes"] = "Q"
            group["Q"] = [0.1, 0.2]
            group["I"] = [10.0, 5.0]
    locations = discover_scattering(source)
    assert len(locations) == 2
    assert {item.variant for item in locations} == {"default", "slit-smeared"}
    np.testing.assert_allclose(load_scattering(locations[0]).intensity, [10, 5])


def test_h5xp_writer_is_part_of_public_io_api(tmp_path):
    path = tmp_path / "public.h5xp"
    with create_h5xp(path, overwrite=True) as handle:
        write_iq_data(handle, "sample", [0.1, 0.2], [10.0, 5.0])
    assert path.is_file()
