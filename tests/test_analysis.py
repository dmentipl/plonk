"""Testing Analysis."""

import pathlib

import pandas as pd

import plonk

TEST_FILE = pathlib.Path(__file__).parent / 'stubdata/phantom_00000.h5'
CSV_FILE = pathlib.Path(__file__).parent / 'stubdata/phantom_00000_profile.csv'


def test_profile():
    """Test Profile."""
    snap = plonk.load_snap(TEST_FILE)
    p = plonk.analysis.Profile(snap)
    for key in (
        'angular_momentum_magnitude',
        'angular_momentum_x',
        'density',
        'mass',
        'number',
        'radius',
        'scale_height',
        'sound_speed',
        'size',
        'smoothing_length',
    ):
        p[key]
    df = p.to_dataframe()

    pd.testing.assert_frame_equal(
        df, pd.read_csv(CSV_FILE, index_col=0), check_less_precise=2
    )
