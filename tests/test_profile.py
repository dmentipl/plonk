"""Test Profile."""

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import plonk

TEST_FILE = Path(__file__).parent / 'stubdata/phantom_00000.h5'
CSV_FILE = Path(__file__).parent / 'stubdata/phantom_00000_profile.csv'


def test_load_profile():
    """Test loading Profile."""
    snap = plonk.load_snap(TEST_FILE)

    prof = plonk.load_profile(snap=snap)
    for p in ['aspect_ratio', 'angular_momentum_phi', 'angular_momentum_theta']:
        prof[p]
    with pytest.raises(ValueError):
        prof['does_not_exist']
    plonk.load_profile(snap=snap, ndim=3, cmin='10 au', cmax='100 au', n_bins=30)
    plonk.load_profile(snap=snap, spacing='log', ignore_accreted=False)
    plonk.load_profile(snap=snap, ndim=1, coordinate='x')
    plonk.load_profile(snap=snap, ndim=1, coordinate='y')
    plonk.load_profile(snap=snap, ndim=1, coordinate='z')
    with pytest.raises(ValueError):
        plonk.load_profile(snap=snap, ndim=1, coordinate='does_not_exist')
    with pytest.raises(ValueError):
        plonk.load_profile(snap=snap, cmin=10, cmax=100)

    snap.close_file()


def test_check_data():
    """Test Profile data accuracy."""
    snap = plonk.load_snap(TEST_FILE)
    prof = plonk.load_profile(snap)

    columns = [
        'density',
        'mass',
        'number',
        'radius',
        'scale_height',
        'size',
        'smoothing_length',
        'sound_speed',
    ]

    df = prof.to_dataframe(columns=columns)

    units = ['g/cm^3', 'g', '', 'au', 'au', 'au^2', 'au', 'km/s']
    df = prof.to_dataframe(columns=columns, units=units)
    pd.testing.assert_frame_equal(df, pd.read_csv(CSV_FILE, index_col=0))

    snap.close_file()


def test_set_data():
    """Test setting array on Profile."""
    snap = plonk.load_snap(TEST_FILE)

    prof = plonk.load_profile(snap=snap)
    prof['array'] = np.arange(len(prof)) * plonk.units.au
    with pytest.raises(ValueError):
        prof['array'] = np.arange(len(prof) - 1)
    with pytest.raises(ValueError):
        prof['array'] = 1.0

    snap.close_file()


def test_profile_plot():
    """Test loading Profile."""
    snap = plonk.load_snap(TEST_FILE)

    prof = plonk.load_profile(snap=snap)
    prof.plot(x='radius', y='surface_density')
    prof.plot(
        x='radius',
        y='density',
        units={'position': 'au', 'density': 'g/cm^3'},
        std_dev_shading=True,
    )

    snap.close_file()


def test_to_function():
    """Test to_function."""
    snap = plonk.load_snap(TEST_FILE)
    prof = plonk.load_profile(snap=snap)

    fn = prof.to_function(profile='scale_height')
    assert np.allclose(fn(prof['radius']), prof['scale_height'])
    snap.close_file()
