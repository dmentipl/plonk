"""Calculate extra quantities on the sinks."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numba
import numpy as np

from .._units import Quantity
from .._units import units as plonk_units

if TYPE_CHECKING:
    from ..snap.snap import Sinks

G = (1 * plonk_units.newtonian_constant_of_gravitation).to_base_units()


def kinetic_energy(sinks: Sinks) -> Quantity:
    """Calculate the total kinetic energy.

    Parameters
    ----------
    sinks
        The Sinks object.

    Returns
    -------
    Quantity
    """
    m, v = sinks['mass'], sinks['velocity']
    return np.sum(1 / 2 * m * _norm(v) ** 2)


def gravitational_potential_energy(sinks: Sinks) -> Quantity:
    """Calculate the total gravitational potential energy.

    Parameters
    ----------
    sinks
        The Sinks object.

    Returns
    -------
    Quantity
    """
    position = sinks['position']
    mass = sinks['mass']
    potential = _potential(position=position.magnitude, mass=mass.magnitude)
    return (potential * G * mass.units ** 2 / position.units).to_base_units()


def specific_orbital_energy(sinks: Sinks) -> Quantity:
    """Calculate the specific orbital energy for two bodies.

    Parameters
    ----------
    sinks
        The Sinks object. Must have length 2.

    Returns
    -------
    Quantity
    """
    if len(sinks) != 2:
        raise ValueError('sinks must have length 2')
    m1, m2 = sinks['mass']
    mu = m1 * m2 / (m1 + m2)
    ke = kinetic_energy(sinks=sinks)
    pe = gravitational_potential_energy(sinks=sinks)
    return (ke + pe) / mu


def specific_angular_momentum(sinks: Sinks) -> Quantity:
    """Calculate the specific orbital energy for two bodies.

    Parameters
    ----------
    sinks
        The Sinks object. Must have length 2.

    Returns
    -------
    Quantity
        The specific angular momentum.
    """
    if len(sinks) != 2:
        raise ValueError('sinks must have length 2')
    r = sinks['position'][1] - sinks['position'][0]
    v = sinks['velocity'][1] - sinks['velocity'][0]
    return np.cross(r, v)


def eccentricity(sinks: Sinks) -> Quantity:
    """Calculate the eccentricity.

    Parameters
    ----------
    sinks
        The Sinks object. Must have length 2.

    Returns
    -------
    Quantity
    """
    if len(sinks) != 2:
        raise ValueError('sinks must have length 2')
    mu = G * np.sum(sinks['mass'])
    eps = specific_orbital_energy(sinks)
    h = specific_angular_momentum(sinks)
    h_mag = _norm(h)
    return np.sqrt(1 + 2 * eps * h_mag ** 2 / mu ** 2)


def semi_major_axis(sinks: Sinks) -> Quantity:
    """Calculate the semi-major axis for two bodies.

    Parameters
    ----------
    sinks
        The Sinks object. Must have length 2.

    Returns
    -------
    Quantity
    """
    if len(sinks) != 2:
        raise ValueError('sinks must have length 2')
    mu = G * np.sum(sinks['mass'])
    h = specific_angular_momentum(sinks)
    h_mag = _norm(h)
    e = eccentricity(sinks)
    return h_mag ** 2 / (mu * (1 - e ** 2))


def inclination(sinks: Sinks) -> Quantity:
    """Calculate the inclination for two bodies.

    Parameters
    ----------
    sinks
        The Sinks object. Must have length 2.

    Returns
    -------
    Quantity
    """
    if len(sinks) != 2:
        raise ValueError('sinks must have length 2')
    h = specific_angular_momentum(sinks)
    h_z = h[..., 2]
    h_mag = _norm(h)
    return np.arccos(h_z / h_mag)


def orbital_period(sinks: Sinks) -> Quantity:
    """Calculate the orbital period for two bodies.

    Parameters
    ----------
    sinks
        The Sinks object. Must have length 2.

    Returns
    -------
    Quantity
    """
    if len(sinks) != 2:
        raise ValueError('sinks must have length 2')
    a = semi_major_axis(sinks)
    mu = G * np.sum(sinks['mass'])
    return 2 * np.pi * np.sqrt(a ** 3 / mu)


def mean_motion(sinks: Sinks) -> Quantity:
    """Calculate the mean motion for two bodies.

    Parameters
    ----------
    sinks
        The Sinks object. Must have length 2.

    Returns
    -------
    Quantity
    """
    if len(sinks) != 2:
        raise ValueError('sinks must have length 2')
    P = orbital_period(sinks)
    return 2 * np.pi / P


def Roche_sphere(sinks: Sinks) -> Quantity:
    """Calculate an estimate of the Roche sphere for two bodies.

    The Roche sphere radius is calculated around the first of the two
    sinks. Uses the formula from Eggleton (1983) ApJ 268, 368-369.

    Parameters
    ----------
    sinks
        The Sinks object. Must have length 2.

    Returns
    -------
    Quantity
    """
    if len(sinks) != 2:
        raise ValueError('sinks must have length 2')
    separation = _norm(sinks['position'][1] - sinks['position'][0])
    q = sinks['mass'][0] / sinks['mass'][1]
    return (
        separation
        * 0.49
        * q ** (2 / 3)
        / (0.6 * q ** (2 / 3) + np.log(1.0 + q ** (1 / 3)))
    )


def Hill_radius(primary: Sinks, secondary: Sinks) -> Quantity:
    """Calculate the Hill radius.

    Parameters
    ----------
    primary
        The primary, i.e. heavy object, as a Sinks object. It must have
        one sink.
    secondary
        The secondary, i.e. light objects, as a Sinks object. It can
        have more than one sink, e.g. when calculating the Hill radiius
        for multiple planets orbiting a star.

    Returns
    -------
    Quantity
    """
    if len(primary) != 1:
        raise ValueError('primary must have length 1')
    Hill = [_Hill_radius(primary + sink) for sink in secondary]
    return [H.magnitude for H in Hill] * Hill[0].units


def _Hill_radius(sinks):
    if len(sinks) != 2:
        raise ValueError('sinks must have length 2')
    M, m = sinks['mass']
    q = m / M if m / M < 1 else M / m
    a = semi_major_axis(sinks)
    e = eccentricity(sinks)
    return a * (1 - e) * (q / 3) ** (1 / 3)


def _norm(x):
    return np.sqrt(x[..., 0] ** 2 + x[..., 1] ** 2 + x[..., 2] ** 2)


@numba.njit
def _potential(position, mass):
    """Get gravitational potential on particles.

    Parameters
    ----------
    position
        Particle positions.
    mass
        Particle masses.

    Returns
    -------
    potential
        The total gravitational potential.
    """
    number_of_particles = len(mass)
    potential = 0.0

    # Loop over all particles...
    for i in range(number_of_particles):
        # ...and all neighbours
        phi = 0.0
        for j in range(number_of_particles):
            # Ignore self
            if j == i:
                continue
            dx = position[i, :] - position[j, :]
            r = np.sqrt(dx[0] ** 2 + dx[1] ** 2 + dx[2] ** 2)
            phi += -mass[j] / r
        potential += mass[i] * phi

    # Account for double counting
    potential /= 2

    return potential
