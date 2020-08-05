"""Testing Analysis."""

import pathlib

import pandas as pd

import plonk

TEST_FILE = pathlib.Path(__file__).parent / 'stubdata/phantom_00000.h5'
CSV_FILE = pathlib.Path(__file__).parent / 'stubdata/phantom_00000_profile.csv'


def test_profile():
    """Test Profile."""
    snap = plonk.load_snap(TEST_FILE)
    prof = plonk.analysis.Profile(snap)
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
    units = ['g/cm^3', 'g', '', 'au', 'au', 'au^2', 'au', 'km/s']
    df = prof.to_dataframe(columns=columns, units=units)

    pd.testing.assert_frame_equal(df, pd.read_csv(CSV_FILE, index_col=0))
