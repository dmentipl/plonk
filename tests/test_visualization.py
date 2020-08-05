"""Testing visualization."""

import pathlib

import plonk

TEST_FILE = pathlib.Path(__file__).parent / 'stubdata/phantom_00000.h5'


def test_initialization():
    """Test plot function."""
    snap = plonk.load_snap(TEST_FILE)
    plonk.plot(
        snap=snap,
        quantity='density',
        x='x',
        y='y',
        extent=(-150, 150, -150, 150),
        norm='linear',
        cmap='gist_heat',
        number_of_pixels=(512, 512),
    )
