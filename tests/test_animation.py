"""Test animation functions."""

import pathlib

import pytest

import plonk


# This test is hanging with NUMBA_DISABLE_JIT=1
@pytest.mark.skip
def test_animation():
    """Test animation of images."""
    dir_path = pathlib.Path(__file__).parent / 'stubdata'
    sim = plonk.load_sim(prefix='phantom', directory=dir_path)

    snaps = [sim.snaps[0], sim.snaps[0], sim.snaps[0]]
    filename = pathlib.Path('animation.mp4')
    plonk.animation(
        filename=filename,
        snaps=snaps,
        quantity='density',
        units={'extent': 'au', 'quantity': 'g/cm^3'},
        adaptive_colorbar=False,
    )
    filename.unlink()


def test_animation_profiles():
    """Test animation of profiles."""
    dir_path = pathlib.Path(__file__).parent / 'stubdata'
    sim = plonk.load_sim(prefix='phantom', directory=dir_path)

    snaps = [sim.snaps[0], sim.snaps[0], sim.snaps[0]]
    profiles = [plonk.load_profile(snap) for snap in snaps]

    filename = pathlib.Path('animation.mp4')
    plonk.animation_profiles(
        filename=filename,
        profiles=profiles,
        x='radius',
        y='surface_density',
        units={'x': 'au', 'y': 'g/cm^2'},
    )
    filename.unlink()


def test_animation_particles():
    """Test animation of particle plots."""
    dir_path = pathlib.Path(__file__).parent / 'stubdata'
    sim = plonk.load_sim(prefix='phantom', directory=dir_path)

    snaps = [sim.snaps[0], sim.snaps[0], sim.snaps[0]]
    filename = pathlib.Path('animation.mp4')
    plonk.animation_particles(
        filename=filename,
        snaps=snaps,
        x='x',
        y='density',
        units={'x': 'au', 'y': 'g/cm^3'},
        adaptive_limits=False,
    )
    filename.unlink()
