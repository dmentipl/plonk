"""Testing Snap."""

from pathlib import Path

import numpy as np
import pytest
from scipy.spatial.transform import Rotation

import plonk

from .data.phantom import adiabatic, dustmixture, dustseparate, mhd

SNAPTYPES = [adiabatic, dustmixture, dustseparate, mhd]
DIR = Path(__file__).parent / 'data/phantom'
RTOL = 1e-6


@pytest.mark.parametrize('snaptype', SNAPTYPES)
def test_load_phantom_snap(snaptype):
    """Testing reading Phantom HDF5 snapshots."""
    # Read from Path
    filename = DIR / snaptype.filename
    snap = plonk.load_snap(filename)
    snap.close_file()
    # Read from str
    snap = plonk.load_snap(str(filename))
    snap.close_file()
    # Not exists
    with pytest.raises(FileNotFoundError):
        plonk.load_snap('does_not_exist.h5')


@pytest.mark.parametrize('snaptype', SNAPTYPES)
def test_get_item(snaptype):
    """Testing getting items from Snap."""
    filename = DIR / snaptype.filename
    snap = plonk.load_snap(filename)

    position = snap['position']
    assert position.shape == snaptype.position_shape

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


@pytest.mark.parametrize('snaptype', SNAPTYPES)
def test_read_particle_arrays_from_phantom(snaptype):
    """Testing reading Phantom HDF5 snapshot particle arrays."""
    filename = DIR / snaptype.filename
    snap = plonk.load_snap(filename)
    _check_arrays(
        snap,
        snaptype.array_name_map,
        snaptype.mean_array_values,
        snaptype.std_array_values,
    )
    snap.close_file()


@pytest.mark.parametrize('snaptype', SNAPTYPES)
def test_read_properties_from_phantom(snaptype):
    """Testing reading Phantom HDF5 snapshot properties."""
    filename = DIR / snaptype.filename
    snap = plonk.load_snap(filename)

    for key, value in snaptype.properties.items():
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
            np.testing.assert_allclose(snap_value, value, rtol=RTOL)
        else:
            assert snap_value == value

    snap.close_file()


@pytest.mark.parametrize('snaptype', SNAPTYPES)
def test_available_loaded_arrays(snaptype):
    """Testing seeing available/loaded arrays on Snap."""
    filename = DIR / snaptype.filename
    snap = plonk.load_snap(filename)

    assert snap.available_arrays() == snaptype.available_arrays

    for arr in [
        'position_x',
        'position_y',
        'position_z',
        'position_mag',
        'h',
        'angular_momentum',
    ]:
        snap[arr]

    assert snap.loaded_arrays() == snaptype.loaded_arrays

    snap.close_file()


@pytest.mark.parametrize('snaptype', SNAPTYPES)
def test_array_code_unit(snaptype):
    """Testing getting array code unit."""
    filename = DIR / snaptype.filename
    snap = plonk.load_snap(filename)

    position_unit = snaptype.length_unit * plonk.units('meter')
    assert snap.array_code_unit('position') == position_unit

    for arr in ['position', 'position_x', 'x']:
        snap.array_code_unit(arr)

    with pytest.raises(ValueError):
        snap.array_code_unit('does_not_exist')

    snap.close_file()


@pytest.mark.parametrize('snaptype', SNAPTYPES)
def test_rotate_snap(snaptype):
    """Testing rotating Snap."""
    filename = DIR / snaptype.filename
    snap = plonk.load_snap(filename)

    snap['position']
    snap['radius_cylindrical']
    if snap.num_sinks > 0:
        snap.sinks['position']
    snap.rotate(axis=(1, 2, 3), angle=np.pi)
    snap.rotate(axis=(1, 2, 3), angle=-np.pi)
    _check_arrays(
        snap,
        snaptype.array_name_map,
        snaptype.mean_array_values,
        snaptype.std_array_values,
    )

    snap.rotate(axis=(1, 2, 3), angle=np.pi)
    snap.reset()
    _check_arrays(
        snap,
        snaptype.array_name_map,
        snaptype.mean_array_values,
        snaptype.std_array_values,
    )

    rot = np.array([1, 2, 3])
    rot = rot / np.linalg.norm(rot)
    rot *= 2 * np.pi
    rotation = Rotation.from_rotvec(rot)
    snap.rotate(rotation=rotation)
    _check_arrays(
        snap,
        snaptype.array_name_map,
        snaptype.mean_array_values,
        snaptype.std_array_values,
    )

    snap.close_file()


@pytest.mark.parametrize('snaptype', SNAPTYPES)
def test_translate_snap(snaptype):
    """Testing translating Snap."""
    filename = DIR / snaptype.filename
    snap = plonk.load_snap(filename)

    unit = f"{snap.code_units['length'].m} {snap.code_units['length'].u}"
    snap['position']
    if snap.num_sinks > 0:
        snap.sinks['position']
    snap.translate(translation=(100, 200, 300), unit=unit)
    snap.translate(translation=(-100, -200, -300), unit=unit)
    _check_arrays(
        snap,
        snaptype.array_name_map,
        snaptype.mean_array_values,
        snaptype.std_array_values,
    )

    snap.translate(translation=(100, 200, 300), unit=unit)
    snap.reset()
    _check_arrays(
        snap,
        snaptype.array_name_map,
        snaptype.mean_array_values,
        snaptype.std_array_values,
    )

    with pytest.raises(ValueError):
        snap.translate(translation=(100, 200, 300))

    with pytest.raises(ValueError):
        snap.translate(translation=(100, 200))

    snap.translate(translation=(100, 200, 300) * plonk.units(unit), unit=unit)
    snap.translate(translation=(-100, -200, -300) * plonk.units(unit), unit=unit)
    _check_arrays(
        snap,
        snaptype.array_name_map,
        snaptype.mean_array_values,
        snaptype.std_array_values,
    )

    snap.close_file()


@pytest.mark.parametrize('snaptype', SNAPTYPES)
def test_write_to_dataframe(snaptype):
    """Testing writing Snap to DataFrame."""
    filename = DIR / snaptype.filename
    snap = plonk.load_snap(filename)

    columns = ['position', 'density', 'smoothing_length']
    snap.to_dataframe(columns=columns)

    columns = ['position', 'density', 'smoothing_length']
    snap.to_dataframe(columns=columns, units=['au', 'g/cm^3', 'au'])

    snap.close_file()


@pytest.mark.parametrize('snaptype', SNAPTYPES)
def test_subsnap(snaptype):
    """Testing getting SubSnap."""
    filename = DIR / snaptype.filename
    snap = plonk.load_snap(filename)

    gas = snap['gas']
    assert len(gas) == snaptype.len_gas

    subsnap = snap[0:100]
    assert len(subsnap) == 100

    subsnap = snap[[0, 1, 2]]
    assert len(subsnap) == 3

    subsnap = snap[(0, 1, 2)]
    assert len(subsnap) == 3

    subsnap = snap[0]
    assert len(subsnap) == 1

    snap.close_file()


@pytest.mark.parametrize('snaptype', SNAPTYPES)
def test_sinks(snaptype):
    """Testing getting sink particles."""
    filename = DIR / snaptype.filename
    snap = plonk.load_snap(filename)

    if snap.num_sinks > 0:
        sinks = snap.sinks

        assert snap.num_sinks == snaptype.num_sinks
        assert len(sinks) == snaptype.num_sinks

    snap.close_file()


@pytest.mark.parametrize('snaptype', SNAPTYPES)
def test_set_array(snaptype):
    """Testing setting array on particles."""
    filename = DIR / snaptype.filename
    snap = plonk.load_snap(filename)

    particle_array = np.arange(len(snap)) * plonk.units('dimensionless')
    snap['array'] = particle_array
    np.testing.assert_allclose(snap['array'].m, particle_array.m, rtol=RTOL)

    sink_array = np.arange(len(snap.sinks)) * plonk.units('dimensionless')
    snap.sinks['array'] = sink_array
    np.testing.assert_allclose(snap.sinks['array'].m, sink_array.m, rtol=RTOL)

    snap.close_file()


@pytest.mark.parametrize('snaptype', SNAPTYPES)
def test_bulk_load(snaptype):
    """Testing bulk loading arrays."""
    filename = DIR / snaptype.filename
    snap = plonk.load_snap(filename)
    snap.bulk_load()

    snap.close_file()


@pytest.mark.parametrize('snaptype', SNAPTYPES)
def test_read_write_extra(snaptype):
    """Testing read write extra arrays."""
    filename = DIR / snaptype.filename
    snap = plonk.load_snap(filename)

    _filename = Path('tmp.h5')

    arr = np.arange(len(snap)) * plonk.units('dimensionless')
    snap['my_array'] = arr
    snap.write_extra_arrays(arrays=['my_array'], filename=_filename)
    snap = None

    snap = plonk.load_snap(filename)
    snap.read_extra_arrays(filename=_filename)
    np.testing.assert_allclose(snap['my_array'], arr, rtol=RTOL)

    _filename.unlink()

    snap.close_file()


@pytest.mark.parametrize('snaptype', SNAPTYPES)
def test_plot_as_methods(snaptype):
    """Testing plot methods."""
    filename = DIR / snaptype.filename
    snap = plonk.load_snap(filename)

    snap.image('density', num_pixels=(16, 16))
    snap.plot()

    if snap.num_sinks > 0:
        sinks = snap.sinks
        sinks.plot()

    subsnap = snap['gas']
    subsnap.image('density', num_pixels=(16, 16))
    subsnap.plot()

    snap.close_file()


@pytest.mark.parametrize('snaptype', SNAPTYPES)
def test_context(snaptype):
    """Testing cache context manager."""
    filename = DIR / snaptype.filename
    snap = plonk.load_snap(filename)

    snap.cache_arrays = True
    with snap.context(cache=False):
        assert snap.loaded_arrays() == []
        with snap.context(cache=True):
            snap.image('density', num_pixels=(16, 16))
            assert snap.loaded_arrays() == [
                'density',
                'mass',
                'position',
                'smoothing_length',
            ]
        assert snap.loaded_arrays() == []
    assert snap.loaded_arrays() == []

    snap.close_file()


def _check_arrays(snap, array_name_map, mean_array_values, std_array_values):
    for array in mean_array_values.keys():
        np.testing.assert_allclose(
            snap.array_in_code_units(array_name_map[array]).mean(),
            mean_array_values[array],
            rtol=RTOL,
        )

    for array in std_array_values.keys():
        np.testing.assert_allclose(
            snap.array_in_code_units(array_name_map[array]).std(),
            std_array_values[array],
            rtol=RTOL,
        )
