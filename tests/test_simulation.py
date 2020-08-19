"""Testing Simulation."""

from pathlib import Path

import numpy as np
import pytest

import plonk


def test_init_simulation():
    """Testing initialising simulation."""
    dir_path = Path(__file__).parent / 'stubdata'
    plonk.load_sim(prefix='phantom', directory=dir_path)
    with pytest.raises(ValueError):
        plonk.load_sim(
            prefix='phantom', directory=dir_path, data_source='not_available'
        )
    with pytest.raises(FileNotFoundError):
        plonk.load_sim(prefix='does_not_exist', directory=dir_path)


def test_sim_data():
    """Testing data in simulation."""
    dir_path = Path(__file__).parent / 'stubdata'
    sim = plonk.load_sim(prefix='phantom', directory=dir_path)

    snaps = sim.snaps
    assert len(snaps) == 1
    assert len(snaps[0]) == 2000

    properties = {
        'adiabatic_index': 1.0,
        'dust_method': 'dust as separate sets of particles',
        'equation_of_state': 'locally isothermal disc',
        'grain_density': [3000.0] * plonk.units('kg/m^3'),
        'grain_size': [0.01] * plonk.units('m'),
        'smoothing_length_factor': 1.0,
        'time': [0.0] * plonk.units('s'),
    }

    for key, val in sim.properties.items():
        if isinstance(val, plonk._units.Quantity):
            assert np.allclose(val.m, properties[key].m)
        else:
            assert sim.properties[key] == properties[key]

    assert sim.paths['time_series_global'][0].name == 'phantom01.ev'


def test_simulation_visualization():
    """Test simulation visualization."""
    dir_path = Path(__file__).parent / 'stubdata'
    sim = plonk.load_sim(prefix='phantom', directory=dir_path)

    viz = sim.visualize(kind='particle', x='x', y='y')
    viz.next()
    viz.prev()


def test_to_array():
    """Testing to_array method."""
    dir_path = Path(__file__).parent / 'stubdata'
    sim = plonk.load_sim(prefix='phantom', directory=dir_path)
    sim.to_array(quantity='density', indices=[0, 1, 2])


def test_set_units_time_series():
    """Test set/unset units time series."""
    dir_path = Path(__file__).parent / 'stubdata'
    sim = plonk.load_sim(prefix='phantom', directory=dir_path)
    sim.unset_units_on_time_series()
    sim.set_units_on_time_series()
