"""Test animation functions."""

from pathlib import Path

import plonk
from plonk import visualize

DIR_PATH = Path(__file__).parent / 'data/phantom'
PREFIX = 'dustseparate'


def test_animate():
    """Test animate."""
    sim = plonk.load_simulation(prefix=PREFIX, directory=DIR_PATH)

    snaps = [sim.snaps[0], sim.snaps[0], sim.snaps[0]]
    filename = Path('animation.mp4')
    plonk.animate(
        filename=filename,
        snaps=snaps,
        quantity='density',
        units={'position': 'au', 'density': 'g/cm^3'},
        adaptive_colorbar=False,
        num_pixels=(32, 32),
    )
    filename.unlink()


def test_animation_images():
    """Test animation of images."""
    sim = plonk.load_simulation(prefix=PREFIX, directory=DIR_PATH)

    snaps = [sim.snaps[0], sim.snaps[0], sim.snaps[0]]
    filename = Path('animation.mp4')
    visualize.animation_images(
        filename=filename,
        snaps=snaps,
        quantity='density',
        units={'position': 'au', 'density': 'g/cm^3'},
        adaptive_colorbar=False,
        num_pixels=(32, 32),
    )
    filename.unlink()


def test_animation_profiles():
    """Test animation of profiles."""
    sim = plonk.load_simulation(prefix=PREFIX, directory=DIR_PATH)

    snaps = [sim.snaps[0], sim.snaps[0], sim.snaps[0]]
    profiles = [plonk.load_profile(snap) for snap in snaps]

    filename = Path('animation.mp4')
    visualize.animation_profiles(
        filename=filename,
        profiles=profiles,
        x='radius',
        y='surface_density',
        units={'position': 'au', 'surface_density': 'g/cm^2'},
    )
    filename.unlink()


def test_animation_particles():
    """Test animation of particle plots."""
    sim = plonk.load_simulation(prefix=PREFIX, directory=DIR_PATH)

    snaps = [sim.snaps[0], sim.snaps[0], sim.snaps[0]]
    filename = Path('animation.mp4')
    visualize.animation_particles(
        filename=filename,
        snaps=snaps,
        x='x',
        y='density',
        units={'position': 'au', 'density': 'g/cm^3'},
        adaptive_limits=False,
    )
    filename.unlink()
