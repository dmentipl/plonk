"""Testing visualization."""

import pathlib

import plonk

AU = plonk.units['au']
TEST_FILE = pathlib.Path(__file__).parent / 'stubdata/phantom_00000.h5'


def test_plot():
    """Test particle plot."""
    snap = plonk.load_snap(TEST_FILE)
    plonk.plot(snap=snap)


def test_plot_with_kwargs():
    """Test particle plot with kwargs."""
    snap = plonk.load_snap(TEST_FILE)
    plonk.plot(
        snap=snap,
        x='x',
        y='y',
        c='density',
        units={'x': 'au', 'y': 'au', 'c': 'g/cm^3'},
    )


def test_image_projection():
    """Test image projection."""
    snap = plonk.load_snap(TEST_FILE)
    plonk.image(snap=snap, quantity='density', number_of_pixels=(32, 32))


def test_image_projection_with_kwargs():
    """Test image projection with kwargs."""
    snap = plonk.load_snap(TEST_FILE)
    plonk.image(
        snap=snap,
        quantity='density',
        x='x',
        y='y',
        units={'extent': 'au', 'quantity': 'g/cm^3'},
        extent=(-150, 150, -150, 150) * AU,
        norm='linear',
        cmap='gist_heat',
        number_of_pixels=(32, 32),
    )


def test_image_on_snap():
    """Test image projection as method on Snap."""
    snap = plonk.load_snap(TEST_FILE)
    snap.image(quantity='density', number_of_pixels=(32, 32))


def test_image_slice():
    """Test image slice."""
    snap = plonk.load_snap(TEST_FILE)
    plonk.image(
        snap=snap, quantity='density', interp='slice', number_of_pixels=(32, 32)
    )


def test_image_slice_with_kwargs():
    """Test image slice with kwargs."""
    snap = plonk.load_snap(TEST_FILE)
    plonk.image(
        snap=snap,
        quantity='density',
        interp='slice',
        x='x',
        y='y',
        units={'extent': 'au', 'quantity': 'g/cm^3'},
        extent=(-150, 150, -150, 150) * AU,
        norm='linear',
        cmap='gist_heat',
        number_of_pixels=(32, 32),
    )
