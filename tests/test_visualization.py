"""Testing visualization."""

from pathlib import Path

import pytest

import plonk
from plonk.utils import visualize

from .data.phantom import adiabatic, dustmixture, dustseparate, mhd

SNAPTYPES = [adiabatic, dustmixture, dustseparate, mhd]
DIR = Path(__file__).parent / 'data/phantom'
AU = plonk.units('au')


@pytest.mark.parametrize('snaptype', SNAPTYPES)
def test_plot(snaptype):
    """Test particle plot."""
    filename = DIR / snaptype.filename
    snap = plonk.load_snap(filename)
    plonk.plot(snap=snap)

    snap.close_file()


@pytest.mark.parametrize('snaptype', SNAPTYPES)
def test_plot_with_kwargs(snaptype):
    """Test particle plot with kwargs."""
    filename = DIR / snaptype.filename
    snap = plonk.load_snap(filename)
    plonk.plot(
        snap=snap,
        x='x',
        y='y',
        c='density',
        units={'position': 'au', 'density': 'g/cm^3', 'projection': 'cm'},
    )

    snap.close_file()


@pytest.mark.parametrize('snaptype', SNAPTYPES)
def test_image_projection(snaptype):
    """Test image projection."""
    filename = DIR / snaptype.filename
    snap = plonk.load_snap(filename)
    plonk.image(snap=snap, quantity='density', num_pixels=(32, 32))

    snap.close_file()


@pytest.mark.parametrize('snaptype', SNAPTYPES)
def test_image_projection_with_kwargs(snaptype):
    """Test image projection with kwargs."""
    filename = DIR / snaptype.filename
    snap = plonk.load_snap(filename)
    plonk.image(
        snap=snap,
        quantity='density',
        x='x',
        y='y',
        units={'position': 'au', 'density': 'g/cm^3', 'projection': 'cm'},
        extent=(-150, 150, -150, 150) * AU,
        norm='linear',
        cmap='gist_heat',
        num_pixels=(32, 32),
    )

    snap.close_file()


@pytest.mark.parametrize('snaptype', SNAPTYPES)
def test_image_on_snap(snaptype):
    """Test image projection as method on Snap."""
    filename = DIR / snaptype.filename
    snap = plonk.load_snap(filename)
    snap.image(quantity='density', num_pixels=(32, 32))

    snap.close_file()


@pytest.mark.parametrize('snaptype', SNAPTYPES)
def test_image_slice(snaptype):
    """Test image slice."""
    filename = DIR / snaptype.filename
    snap = plonk.load_snap(filename)
    plonk.image(snap=snap, quantity='density', interp='slice', num_pixels=(32, 32))

    snap.close_file()


@pytest.mark.parametrize('snaptype', SNAPTYPES)
def test_image_slice_with_kwargs(snaptype):
    """Test image slice with kwargs."""
    filename = DIR / snaptype.filename
    snap = plonk.load_snap(filename)
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
        num_pixels=(32, 32),
    )

    snap.close_file()


@pytest.mark.parametrize('snaptype', SNAPTYPES)
def test_plot_smoothing_length(snaptype):
    """Test plot smoothing length as circle."""
    filename = DIR / snaptype.filename
    snap = plonk.load_snap(filename)

    visualize.plot_smoothing_length(snap=snap, indices=[0, 1])

    snap.close_file()


@pytest.mark.parametrize('snaptype', SNAPTYPES)
def test_get_extent(snaptype):
    """Test getting extent from percentile."""
    filename = DIR / snaptype.filename
    snap = plonk.load_snap(filename)

    visualize.get_extent_from_percentile(snap=snap, x='x', y='y')

    snap.close_file()
