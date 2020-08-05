"""Testing Snap."""

import pathlib

import numpy as np
import pytest

import plonk
from plonk.snap.utils import get_array_in_code_units

from .stubdata.phantom_snapshot import (
    array_name_map,
    mean_array_values,
    properties,
    std_array_values,
)

TEST_FILE = pathlib.Path(__file__).parent / 'stubdata/phantom_00000.h5'


def test_load_phantom_snap():
    """Testing reading Phantom HDF5 snapshots."""
    # Read from pathlib.Path
    snap = plonk.load_snap(TEST_FILE)
    snap.close_file()
    # Read from str
    snap = plonk.load_snap(str(TEST_FILE))
    snap.close_file()
    # Not exists
    with pytest.raises(FileNotFoundError):
        plonk.load_snap('does_not_exist.h5')


def test_read_particle_arrays_from_phantom():
    """Testing reading Phantom HDF5 snapshot particle arrays."""
    snap = plonk.load_snap(TEST_FILE)

    for array in mean_array_values.keys():
        np.testing.assert_allclose(
            get_array_in_code_units(snap, array_name_map[array]).mean(),
            mean_array_values[array],
        )

    for array in std_array_values.keys():
        np.testing.assert_allclose(
            get_array_in_code_units(snap, array_name_map[array]).std(),
            std_array_values[array],
        )

    snap.close_file()


def test_read_properties_from_phantom():
    """Testing reading Phantom HDF5 snapshot properties."""
    snap = plonk.load_snap(TEST_FILE)

    for key, value in properties.items():
        if isinstance(snap.properties[key], plonk.units.Quantity):
            snap_value = snap.properties[key].magnitude
        else:
            snap_value = snap.properties[key]
        np.testing.assert_allclose(snap_value, value)

    snap.close_file()
