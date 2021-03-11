"""Test Profile."""

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import plonk

from .data.phantom import adiabatic, dustmixture, dustseparate, mhd

SNAPTYPES = [adiabatic, dustmixture, dustseparate, mhd]
DIR = Path(__file__).parent / 'data/phantom'


@pytest.mark.parametrize('snaptype', SNAPTYPES)
def test_load_profile(snaptype):
    """Test loading Profile."""
    filename = DIR / snaptype.filename
    snap = plonk.load_snap(filename)

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


@pytest.mark.parametrize('snaptype', SNAPTYPES)
def test_check_data(snaptype):
    """Test Profile data accuracy."""
    filename = DIR / snaptype.filename
    snap = plonk.load_snap(filename)
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
    profile_file = DIR / snaptype.profile_file
    pd.testing.assert_frame_equal(df, pd.read_csv(profile_file, index_col=0))

    snap.close_file()


@pytest.mark.parametrize('snaptype', SNAPTYPES)
def test_set_data(snaptype):
    """Test setting array on Profile."""
    filename = DIR / snaptype.filename
    snap = plonk.load_snap(filename)

    prof = plonk.load_profile(snap=snap)
    prof['array'] = np.arange(len(prof)) * plonk.units.au
    with pytest.raises(ValueError):
        prof['array'] = np.arange(len(prof) - 1)
    with pytest.raises(ValueError):
        prof['array'] = 1.0

    snap.close_file()


@pytest.mark.parametrize('snaptype', SNAPTYPES)
def test_profile_plot(snaptype):
    """Test loading Profile."""
    filename = DIR / snaptype.filename
    snap = plonk.load_snap(filename)

    prof = plonk.load_profile(snap=snap)
    prof.plot(x='radius', y='surface_density')
    prof.plot(
        x='radius',
        y='density',
        units={'position': 'au', 'density': 'g/cm^3'},
        std='shading',
    )

    snap.close_file()


@pytest.mark.parametrize('snaptype', SNAPTYPES)
def test_to_function(snaptype):
    """Test to_function."""
    filename = DIR / snaptype.filename
    snap = plonk.load_snap(filename)
    prof = plonk.load_profile(snap=snap)

    fn = prof.to_function(profile='scale_height')
    assert np.allclose(fn(prof['radius']), prof['scale_height'])
    snap.close_file()


@pytest.mark.parametrize('snaptype', SNAPTYPES)
def test_alias(snaptype):
    """Test using profile aliases."""
    filename = DIR / snaptype.filename
    snap = plonk.load_snap(filename)
    prof = plonk.load_profile(snap=snap)

    prof.add_alias('scale_height', 'H')
    prof['H']

    snap.close_file()
