"""Testing Snap."""

from pathlib import Path

import numpy as np
import pytest
from scipy.spatial.transform import Rotation

import plonk

from .stubdata.phantom_snapshot import (
    array_name_map,
    mean_array_values,
    properties,
    std_array_values,
)

TEST_FILE = Path(__file__).parent / 'stubdata/phantom_00000.h5'
AVAILABLE_ARRAYS = [
    'angular_momentum',
    'angular_velocity',
    'azimuthal_angle',
    'density',
    'dust_to_gas_ratio',
    'eccentricity',
    'id',
    'inclination',
    'keplerian_frequency',
    'kinetic_energy',
    'mass',
    'momentum',
    'polar_angle',
    'position',
    'pressure',
    'radius_cylindrical',
    'radius_spherical',
    'semi_major_axis',
    'smoothing_length',
    'sound_speed',
    'specific_angular_momentum',
    'stokes_number',
    'stopping_time',
    'sub_type',
    'temperature',
    'timestep',
    'type',
    'velocity',
    'velocity_divergence',
    'velocity_radial_cylindrical',
    'velocity_radial_spherical',
]


def test_load_phantom_snap():
    """Testing reading Phantom HDF5 snapshots."""
    # Read from Path
    snap = plonk.load_snap(TEST_FILE)
    snap.close_file()
    # Read from str
    snap = plonk.load_snap(str(TEST_FILE))
    snap.close_file()
    # Not exists
    with pytest.raises(FileNotFoundError):
        plonk.load_snap('does_not_exist.h5')


def test_get_item():
    """Testing getting items from Snap."""
    snap = plonk.load_snap(TEST_FILE)

    position = snap['position']
    assert position.shape == (2000, 3)

    subsnap = snap['gas']
    assert type(subsnap) == plonk.snap.snap.SubSnap

    subsnap = snap[:10]
    assert type(subsnap) == plonk.snap.snap.SubSnap
    assert len(subsnap) == 10

    subsnap = snap[[0, 1, 2, 3, 4, 5, 6, 7, 8, 9]]
    assert type(subsnap) == plonk.snap.snap.SubSnap
    assert len(subsnap) == 10

    subsnap = snap[(0, 1, 2, 3, 4, 5, 6, 7, 8, 9)]
    assert type(subsnap) == plonk.snap.snap.SubSnap
    assert len(subsnap) == 10

    snap.close_file()


def test_read_particle_arrays_from_phantom():
    """Testing reading Phantom HDF5 snapshot particle arrays."""
    snap = plonk.load_snap(TEST_FILE)
    _check_arrays(snap)
    snap.close_file()


def test_read_properties_from_phantom():
    """Testing reading Phantom HDF5 snapshot properties."""
    snap = plonk.load_snap(TEST_FILE)

    for key, value in properties.items():
        if isinstance(snap.properties[key], plonk.units.Quantity):
            snap_value = snap.properties[key].magnitude
            numpy_array = True
        elif isinstance(snap.properties[key], str):
            snap_value = snap.properties[key]
            numpy_array = False
        else:
            snap_value = snap.properties[key]
            numpy_array = True
        if numpy_array:
            np.testing.assert_allclose(snap_value, value)
        else:
            assert snap_value == value

    snap.close_file()


def test_available_loaded_arrays():
    """Testing seeing available/loaded arrays on Snap."""
    snap = plonk.load_snap(TEST_FILE)

    assert snap.available_arrays() == AVAILABLE_ARRAYS

    for arr in [
        'position_x',
        'position_y',
        'position_z',
        'position_mag',
        'h',
        'dust_to_gas_ratio_tot',
    ]:
        snap[arr]

    assert snap.loaded_arrays() == ['dust_to_gas_ratio', 'position', 'smoothing_length']

    snap.close_file()


def test_array_code_unit():
    """Testing getting array code unit."""
    snap = plonk.load_snap(TEST_FILE)

    position_unit = 149600000000.0 * plonk.units['meter']
    assert snap.array_code_unit('position') == position_unit

    for arr in ['position', 'position_x', 'x']:
        snap.array_code_unit(arr)

    with pytest.raises(ValueError):
        snap.array_code_unit('does_not_exist')

    snap.close_file()


def test_rotate_snap():
    """Testing rotating Snap."""
    snap = plonk.load_snap(TEST_FILE)

    snap['position']
    snap['radius_cylindrical']
    snap.sinks['position']
    snap.rotate(axis=(1, 2, 3), angle=np.pi)
    snap.rotate(axis=(1, 2, 3), angle=-np.pi)
    _check_arrays(snap)

    snap.rotate(axis=(1, 2, 3), angle=np.pi)
    snap.reset()
    _check_arrays(snap)

    rot = np.array([1, 2, 3])
    rot = rot / np.linalg.norm(rot)
    rot *= np.pi
    rotation = Rotation.from_rotvec(rot)
    snap.rotate(rotation=rotation)
    _check_arrays(snap)

    snap.close_file()


def test_translate_snap():
    """Testing translating Snap."""
    snap = plonk.load_snap(TEST_FILE)

    snap['position']
    snap.sinks['position']
    snap.translate(translation=(100, 200, 300), unit='au')
    snap.translate(translation=(-100, -200, -300), unit='au')
    _check_arrays(snap)

    snap.translate(translation=(100, 200, 300), unit='au')
    snap.reset()
    _check_arrays(snap)

    with pytest.raises(ValueError):
        snap.translate(translation=(100, 200, 300))

    with pytest.raises(ValueError):
        snap.translate(translation=(100, 200))

    snap.translate(translation=(100, 200, 300) * plonk.units['au'], unit='au')
    snap.translate(translation=(-100, -200, -300) * plonk.units['au'], unit='au')
    _check_arrays(snap)

    snap.close_file()


def test_write_to_dataframe():
    """Testing writing Snap to DataFrame."""
    snap = plonk.load_snap(TEST_FILE)

    columns = ['position', 'density', 'smoothing_length']
    snap.to_dataframe(columns=columns)

    columns = ['position', 'density', 'smoothing_length']
    snap.to_dataframe(columns=columns, units=['au', 'g/cm^3', 'au'])

    snap.close_file()


def test_subsnap():
    """Testing getting SubSnap."""
    snap = plonk.load_snap(TEST_FILE)

    gas = snap['gas']
    assert len(gas) == 1000

    subsnap = snap[0:100]
    assert len(subsnap) == 100

    subsnap = snap[[0, 1, 2]]
    assert len(subsnap) == 3

    subsnap = snap[(0, 1, 2)]
    assert len(subsnap) == 3

    subsnap = snap[0]
    assert len(subsnap) == 1

    snap.close_file()


def test_sinks():
    """Testing getting sink particles."""
    snap = plonk.load_snap(TEST_FILE)

    sinks = snap.sinks

    assert snap.num_sinks == 1
    assert len(sinks) == 1

    snap.close_file()


def test_set_array():
    """Testing setting array on particles."""
    snap = plonk.load_snap(TEST_FILE)

    particle_array = np.arange(len(snap)) * plonk.units['dimensionless']
    snap['array'] = particle_array
    np.testing.assert_allclose(snap['array'].m, particle_array.m)

    sink_array = np.arange(len(snap.sinks)) * plonk.units['dimensionless']
    snap.sinks['array'] = sink_array
    np.testing.assert_allclose(snap.sinks['array'].m, sink_array.m)

    snap.close_file()


def test_bulk_load():
    """Testing bulk loading arrays."""
    snap = plonk.load_snap(TEST_FILE)
    snap.bulk_load()

    snap.close_file()


def test_read_write_extra():
    """Testing read write extra arrays."""
    snap = plonk.load_snap(TEST_FILE)

    filename = Path('tmp.h5')

    arr = np.arange(len(snap)) * plonk.units['dimensionless']
    snap['my_array'] = arr
    snap.write_extra_arrays(arrays=['my_array'], filename=filename)
    snap = None

    snap = plonk.load_snap(TEST_FILE)
    snap.read_extra_arrays(filename=filename)
    np.allclose(snap['my_array'], arr)

    filename.unlink()

    snap.close_file()


def test_plot_as_methods():
    """Testing plot methods."""
    snap = plonk.load_snap(TEST_FILE)

    snap.image('density', number_of_pixels=(16, 16))
    snap.plot()

    sinks = snap.sinks
    sinks.plot()

    subsnap = snap['gas']
    subsnap.image('density', number_of_pixels=(16, 16))
    subsnap.plot()

    snap.close_file()


def test_context():
    """Testing cache context manager."""
    snap = plonk.load_snap(TEST_FILE)

    snap.cache_arrays = True
    with snap.context(cache=False):
        assert snap.loaded_arrays() == []
        with snap.context(cache=True):
            snap.image('density', number_of_pixels=(16, 16))
            assert snap.loaded_arrays() == [
                'density',
                'mass',
                'position',
                'smoothing_length',
            ]
        assert snap.loaded_arrays() == []
    assert snap.loaded_arrays() == []

    snap.close_file()


def _check_arrays(snap):
    for array in mean_array_values.keys():
        np.testing.assert_allclose(
            snap.array_in_code_units(array_name_map[array]).mean(),
            mean_array_values[array],
        )

    for array in std_array_values.keys():
        np.testing.assert_allclose(
            snap.array_in_code_units(array_name_map[array]).std(),
            std_array_values[array],
        )
