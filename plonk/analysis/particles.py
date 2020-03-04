"""Particle quantities.

Calculate various quantities on the particles.
"""

from typing import Tuple, Union

import numpy as np
from numpy import ndarray

from ..snap.snap import Snap, SubSnap

SnapLike = Union[Snap, SubSnap]


def momentum(snap: SnapLike, ignore_accreted: bool = False) -> ndarray:
    """Calculate the momentum.

    Parameters
    ----------
    snap
        The Snap object.
    ignore_accreted : optional
        Ignore accreted particles. Default is False.

    Returns
    -------
    ndarray
        The linear momentum on the particles.
    """
    if ignore_accreted:
        h: ndarray = snap['smooth']
        mass: ndarray = snap['mass'][h > 0]
        vel: ndarray = snap['velocity'][h > 0]
    else:
        mass = snap['mass']
        vel = snap['velocity']

    return mass[:, np.newaxis] * vel


def angular_momentum(
    snap: SnapLike,
    origin: Union[ndarray, Tuple[float, float, float]] = (0.0, 0.0, 0.0),
    ignore_accreted: bool = False,
) -> ndarray:
    """Calculate the angular momentum.

    Parameters
    ----------
    snap
        The Snap object.
    origin : optional
        The origin around which to compute the angular momentum as a
        ndarray or tuple (x, y, z). Default is (0, 0, 0).
    ignore_accreted : optional
        Ignore accreted particles. Default is False.

    Returns
    -------
    ndarray
        The angular momentum on the particles.
    """
    if ignore_accreted:
        h: ndarray = snap['smooth']
        mass: ndarray = snap['mass'][h > 0]
        pos: ndarray = snap['position'][h > 0]
        vel: ndarray = snap['velocity'][h > 0]
    else:
        mass = snap['mass']
        pos = snap['position']
        vel = snap['velocity']

    origin = np.array(origin)
    pos = pos - origin

    return mass[:, np.newaxis] * np.cross(pos, vel)


def specific_angular_momentum(
    snap: SnapLike,
    origin: Union[ndarray, Tuple[float, float, float]] = (0.0, 0.0, 0.0),
    ignore_accreted: bool = False,
) -> ndarray:
    """Calculate the specific angular momentum.

    Parameters
    ----------
    snap
        The Snap object.
    origin : optional
        The origin around which to compute the angular momentum as a
        ndarray or tuple (x, y, z). Default is (0, 0, 0).
    ignore_accreted : optional
        Ignore accreted particles. Default is False.

    Returns
    -------
    ndarray
        The specific angular momentum on the particles.
    """
    if ignore_accreted:
        h: ndarray = snap['smooth']
        pos: ndarray = snap['position'][h > 0]
        vel: ndarray = snap['velocity'][h > 0]
    else:
        pos = snap['position']
        vel = snap['velocity']

    origin = np.array(origin)
    pos = pos - origin

    return np.cross(pos, vel)


def kinetic_energy(snap: SnapLike, ignore_accreted: bool = False) -> ndarray:
    """Calculate the kinetic energy.

    Parameters
    ----------
    snap
        The Snap object.
    ignore_accreted : optional
        Ignore accreted particles. Default is False.

    Returns
    -------
    ndarray
        The kinetic energy on the particles.
    """
    if ignore_accreted:
        h: ndarray = snap['smooth']
        mass: ndarray = snap['mass'][h > 0]
        vel: ndarray = snap['velocity'][h > 0]
    else:
        mass = snap['mass']
        vel = snap['velocity']

    return 1 / 2 * mass * np.linalg.norm(vel, axis=1) ** 2


def specific_kinetic_energy(snap: SnapLike, ignore_accreted: bool = False) -> ndarray:
    """Calculate the specific kinetic energy.

    Parameters
    ----------
    snap
        The Snap object.
    ignore_accreted : optional
        Ignore accreted particles. Default is False.

    Returns
    -------
    ndarray
        The specific kinetic energy on the particles.
    """
    if ignore_accreted:
        h: ndarray = snap['smooth']
        vel: ndarray = snap['velocity'][h > 0]
    else:
        vel = snap['velocity']

    return 1 / 2 * np.linalg.norm(vel, axis=1) ** 2


def semi_major_axis(
    snap: SnapLike,
    gravitational_parameter: float,
    origin: Union[ndarray, Tuple[float, float, float]] = (0.0, 0.0, 0.0),
    ignore_accreted: bool = False,
) -> ndarray:
    """Calculate the semi-major axis.

    The semi-major axis of particles around a mass specified by
    gravitational parameter with an optional to specify the position of
    the mass.

    Parameters
    ----------
    snap
        The Snap object.
    gravitational_parameter
        The gravitational parameter (G*M).
    origin : optional
        The origin around which to compute the angular momentum as a
        ndarray or tuple (x, y, z). Default is (0, 0, 0).
    ignore_accreted : optional
        Ignore accreted particles. Default is False.

    Returns
    -------
    ndarray
        The semi-major axis on the particles.
    """
    if ignore_accreted:
        h: ndarray = snap['smooth']
        pos: ndarray = snap['position'][h > 0]
        vel: ndarray = snap['velocity'][h > 0]
    else:
        pos = snap['position']
        vel = snap['velocity']

    origin = np.array(origin)
    pos = pos - origin

    mu = gravitational_parameter

    radius = np.linalg.norm(pos, axis=1)

    specific_angular_momentum = np.cross(pos, vel)
    specific_angular_momentum_magnitude = np.linalg.norm(
        specific_angular_momentum, axis=1
    )

    specific_kinetic_energy = 1 / 2 * np.linalg.norm(vel, axis=1) ** 2
    specific_potential_energy = -mu / radius
    specific_energy = specific_kinetic_energy + specific_potential_energy

    eccentricity = np.sqrt(
        1 + 2 * specific_energy * (specific_angular_momentum_magnitude / mu) ** 2
    )

    return specific_angular_momentum_magnitude ** 2 / (mu * (1 - eccentricity ** 2))


def eccentricity(
    snap: SnapLike,
    gravitational_parameter: float,
    origin: Union[ndarray, Tuple[float, float, float]] = (0.0, 0.0, 0.0),
    ignore_accreted: bool = False,
) -> ndarray:
    """Calculate the eccentricity.

    The eccentricity of particles around a mass specified by
    gravitational parameter with an optional to specify the position of
    the mass.

    Parameters
    ----------
    snap
        The Snap object.
    gravitational_parameter
        The gravitational parameter (G*M).
    origin : optional
        The origin around which to compute the angular momentum as a
        ndarray or tuple (x, y, z). Default is (0, 0, 0).
    ignore_accreted : optional
        Ignore accreted particles. Default is False.

    Returns
    -------
    ndarray
        The eccentricity on the particles.
    """
    if ignore_accreted:
        h: ndarray = snap['smooth']
        pos: ndarray = snap['position'][h > 0]
        vel: ndarray = snap['velocity'][h > 0]
    else:
        pos = snap['position']
        vel = snap['velocity']

    origin = np.array(origin)
    pos = pos - origin

    mu = gravitational_parameter

    radius = np.linalg.norm(pos, axis=1)

    specific_angular_momentum = np.cross(pos, vel)
    specific_angular_momentum_magnitude = np.linalg.norm(
        specific_angular_momentum, axis=1
    )

    specific_kinetic_energy = 1 / 2 * np.linalg.norm(vel, axis=1) ** 2
    specific_potential_energy = -mu / radius
    specific_energy = specific_kinetic_energy + specific_potential_energy

    return np.sqrt(
        1 + 2 * specific_energy * (specific_angular_momentum_magnitude / mu) ** 2
    )


def inclination(snap: SnapLike, ignore_accreted: bool = False) -> ndarray:
    """Calculate the inclination with respect to the xy-plane.

    The inclination is calculated by taking the angle between the
    angular momentum vector and the z-axis, with the angular momentum
    calculated with respect to the center of mass.

    Parameters
    ----------
    snap
        The Snap object.
    ignore_accreted : optional
        Ignore accreted particles. Default is False.

    Returns
    -------
    ndarray
        The inclination on the particles.
    """
    if ignore_accreted:
        h: ndarray = snap['smooth']
        mass: ndarray = snap['mass'][h > 0]
        pos: ndarray = snap['position'][h > 0]
        vel: ndarray = snap['velocity'][h > 0]
    else:
        mass = snap['mass']
        pos = snap['position']
        vel = snap['velocity']

    origin = (mass[:, np.newaxis] * pos).sum(axis=0)
    pos = pos - origin

    specific_angular_momentum = np.cross(pos, vel)

    inclination = np.arccos(
        specific_angular_momentum[:, 2]
        / np.linalg.norm(specific_angular_momentum, axis=1)
    )

    return inclination


def gas_mass(snap: SnapLike, ignore_accreted: bool = False) -> ndarray:
    """Calculate the gas mass from the dust fraction.

    Parameters
    ----------
    snap
        The Snap object.
    ignore_accreted : optional
        Ignore accreted particles. Default is False.

    Returns
    -------
    ndarray
        The gas mass on the particles.
    """
    if ignore_accreted:
        h: ndarray = snap['smooth']
        dustfrac: ndarray = snap['dustfrac'][h > 0]
        mass: ndarray = snap['mass'][h > 0]
    else:
        dustfrac = snap['dustfrac']
        mass = snap['mass']

    gasfrac = 1 - dustfrac.sum(axis=1)
    return gasfrac * mass


def dust_mass(snap: SnapLike, ignore_accreted: bool = False) -> ndarray:
    """Calculate the dust mass per species from the dust fraction.

    Parameters
    ----------
    snap
        The Snap object.
    ignore_accreted : optional
        Ignore accreted particles. Default is False.

    Returns
    -------
    ndarray
        The dust mass per species on the particles.
    """
    if ignore_accreted:
        h: ndarray = snap['smooth']
        dustfrac: ndarray = snap['dustfrac'][h > 0]
        mass: ndarray = snap['mass'][h > 0]
    else:
        dustfrac = snap['dustfrac']
        mass = snap['mass']

    return dustfrac * mass[:, np.newaxis]


def gas_density(snap: SnapLike, ignore_accreted: bool = False) -> ndarray:
    """Calculate the gas density from the dust fraction.

    Parameters
    ----------
    snap
        The Snap object.
    ignore_accreted : optional
        Ignore accreted particles. Default is False.

    Returns
    -------
    ndarray
        The gas density on the particles.
    """
    if ignore_accreted:
        h: ndarray = snap['smooth']
        dustfrac: ndarray = snap['dustfrac'][h > 0]
        density: ndarray = snap['density'][h > 0]
    else:
        dustfrac = snap['dustfrac']
        density = snap['density']

    gasfrac = 1 - dustfrac.sum(axis=1)
    return gasfrac * density


def dust_density(snap: SnapLike, ignore_accreted: bool = False) -> ndarray:
    """Calculate the dust density per species from the dust fraction.

    Parameters
    ----------
    snap
        The Snap object.
    ignore_accreted : optional
        Ignore accreted particles. Default is False.

    Returns
    -------
    ndarray
        The dust density per species on the particles.
    """
    if ignore_accreted:
        h: ndarray = snap['smooth']
        dustfrac: ndarray = snap['dustfrac'][h > 0]
        density: ndarray = snap['density'][h > 0]
    else:
        dustfrac = snap['dustfrac']
        density = snap['density']

    return dustfrac * density[:, np.newaxis]
