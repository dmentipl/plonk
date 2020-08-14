"""Testing visualization."""

from pathlib import Path

import plonk
from plonk import visualize

AU = plonk.units['au']
TEST_FILE = Path(__file__).parent / 'stubdata/phantom_00000.h5'


def test_plot():
    """Test particle plot."""
    snap = plonk.load_snap(TEST_FILE)
    plonk.plot(snap=snap)

    snap.close_file()


def test_plot_with_kwargs():
    """Test particle plot with kwargs."""
    snap = plonk.load_snap(TEST_FILE)
    plonk.plot(
        snap=snap,
        x='x',
        y='y',
        c='density',
        units={'position': 'au', 'density': 'g/cm^3', 'projection': 'cm'},
    )

    snap.close_file()


def test_image_projection():
    """Test image projection."""
    snap = plonk.load_snap(TEST_FILE)
    plonk.image(snap=snap, quantity='density', number_of_pixels=(32, 32))

    snap.close_file()


def test_image_projection_with_kwargs():
    """Test image projection with kwargs."""
    snap = plonk.load_snap(TEST_FILE)
    plonk.image(
        snap=snap,
        quantity='density',
        x='x',
        y='y',
        units={'position': 'au', 'density': 'g/cm^3', 'projection': 'cm'},
        extent=(-150, 150, -150, 150) * AU,
        norm='linear',
        cmap='gist_heat',
        number_of_pixels=(32, 32),
    )

    snap.close_file()


def test_image_on_snap():
    """Test image projection as method on Snap."""
    snap = plonk.load_snap(TEST_FILE)
    snap.image(quantity='density', number_of_pixels=(32, 32))

    snap.close_file()


def test_image_slice():
    """Test image slice."""
    snap = plonk.load_snap(TEST_FILE)
    plonk.image(
        snap=snap, quantity='density', interp='slice', number_of_pixels=(32, 32)
    )

    snap.close_file()


def test_image_slice_with_kwargs():
    """Test image slice with kwargs."""
    snap = plonk.load_snap(TEST_FILE)
    plonk.image(
        snap=snap,
        quantity='density',
        interp='slice',
        x='x',
        y='y',
        units={'position': 'au', 'density': 'g/cm^3', 'projection': 'cm'},
        extent=(-150, 150, -150, 150) * AU,
        norm='linear',
        cmap='gist_heat',
        number_of_pixels=(32, 32),
    )

    snap.close_file()


def test_plot_smoothing_length():
    """Test plot smoothing length as circle."""
    snap = plonk.load_snap(TEST_FILE)

    visualize.plot_smoothing_length(snap=snap, indices=[0, 1])

    snap.close_file()


def test_get_extent():
    """Test getting extent from percentile."""
    snap = plonk.load_snap(TEST_FILE)

    visualize.get_extent_from_percentile(snap=snap, x='x', y='y')

    snap.close_file()
