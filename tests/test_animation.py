"""Test animation functions."""

from pathlib import Path

import plonk


def test_animation():
    """Test animation of images."""
    dir_path = Path(__file__).parent / 'stubdata'
    sim = plonk.load_sim(prefix='phantom', directory=dir_path)

    snaps = [sim.snaps[0], sim.snaps[0], sim.snaps[0]]
    filename = Path('animation.mp4')
    plonk.animation(
        filename=filename,
        snaps=snaps,
        quantity='density',
        units={'position': 'au', 'density': 'g/cm^3'},
        adaptive_colorbar=False,
        number_of_pixels=(32, 32),
    )
    filename.unlink()


def test_animation_profiles():
    """Test animation of profiles."""
    dir_path = Path(__file__).parent / 'stubdata'
    sim = plonk.load_sim(prefix='phantom', directory=dir_path)

    snaps = [sim.snaps[0], sim.snaps[0], sim.snaps[0]]
    profiles = [plonk.load_profile(snap) for snap in snaps]

    filename = Path('animation.mp4')
    plonk.animation_profiles(
        filename=filename,
        profiles=profiles,
        x='radius',
        y='surface_density',
        units={'position': 'au', 'surface_density': 'g/cm^2'},
    )
    filename.unlink()


def test_animation_particles():
    """Test animation of particle plots."""
    dir_path = Path(__file__).parent / 'stubdata'
    sim = plonk.load_sim(prefix='phantom', directory=dir_path)

    snaps = [sim.snaps[0], sim.snaps[0], sim.snaps[0]]
    filename = Path('animation.mp4')
    plonk.animation_particles(
        filename=filename,
        snaps=snaps,
        x='x',
        y='density',
        units={'position': 'au', 'density': 'g/cm^3'},
        adaptive_limits=False,
    )
    filename.unlink()
