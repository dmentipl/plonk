"""
Testing Visualization.
"""

import pathlib

import plonk

TEST_FILE = pathlib.Path(__file__).parent / 'stubdata/phantom_00000.h5'


def test_initialization():
    """Test plot function."""
    snap = plonk.load_snap(TEST_FILE)

    extent = (-150, 150, -150, 150)
    scalar_options = {'norm': 'lin', 'cmap': 'gist_heat'}
    interpolation_options = {
        'number_of_pixels': (512, 512),
        'cross_section': None,
        'density_weighted': False,
    }

    plonk.visualize.plot(
        scalar_data=snap['rho'],
        x_coordinate=snap['x'],
        y_coordinate=snap['y'],
        extent=extent,
        particle_mass=snap['m'],
        smoothing_length=snap['h'],
        scalar_options=scalar_options,
        interpolation_options=interpolation_options,
    )
