"""Analysis for accretion discs."""

from __future__ import annotations

from typing import TYPE_CHECKING, List

import numpy as np
from numpy import ndarray

from .._units import Quantity
from .._units import units as plonk_units
from ..utils.math import cross, norm
from .total import angular_momentum, center_of_mass

if TYPE_CHECKING:
    from ..snap.snap import SnapLike

ORIGIN = (0, 0, 0) * plonk_units.au

# Derived quantities require some arrays already present
# Ignoring type, sub_type, (optional) smoothing_length
array_requires = {
    'eccentricity': ['position', 'velocity'],
    'inclination': ['position', 'velocity'],
    'keplerian_frequency': ['position'],
    'semi_major_axis': ['position', 'velocity'],
    'stokes_number': ['position', 'stopping_time'],
}

# Arrays which represent quantities with x-, y-, z-components in space
vector_arrays: List[str] = []

# Arrays which represent dust quantities with columns per dust species
dust_arrays = ['stokes_number']


def eccentricity(
    snap: SnapLike,
    gravitational_parameter: Quantity = None,
    origin: Quantity = None,
    ignore_accreted: bool = False,
) -> Quantity:
    """Calculate the eccentricity.

    The eccentricity of particles around a mass specified by
    gravitational parameter with an optional to specify the position of
    the mass.

    Parameters
    ----------
    snap
        The Snap object.
    gravitational_parameter
        The gravitational parameter (mu = G M) as a Pint quantity.
    origin : optional
        The origin around which to compute the eccentricity as a
        Quantity like (x, y, z) * au. Default is (0, 0, 0).
    ignore_accreted : optional
        Ignore accreted particles. Default is False.

    Returns
    -------
    Quantity
        The eccentricity on the particles.
    """
    if ignore_accreted:
        h: Quantity = snap['smoothing_length']
        pos: Quantity = snap['position'][h > 0]
        vel: Quantity = snap['velocity'][h > 0]
    else:
        pos = snap['position']
        vel = snap['velocity']

    origin = snap.translation if snap.translation is not None else ORIGIN
    pos = pos - origin

    mu = snap.properties.get('gravitational_parameter')
    if mu is None:
        raise ValueError(
            'must pass in gravitational_parameter or '
            'set on Snap with set_gravitational_parameter'
        )

    radius = norm(pos, axis=1)

    _specific_angular_momentum = cross(pos, vel)
    specific_angular_momentum_magnitude = norm(_specific_angular_momentum, axis=1)

    _specific_kinetic_energy = 1 / 2 * norm(vel, axis=1) ** 2
    specific_potential_energy = -mu / radius
    specific_energy = _specific_kinetic_energy + specific_potential_energy

    term = specific_energy * (specific_angular_momentum_magnitude / mu) ** 2
    return np.sqrt(1 + 2 * term)


def inclination(
    snap: SnapLike, origin: Quantity = None, ignore_accreted: bool = False,
) -> Quantity:
    """Calculate the inclination.

    The inclination is calculated by taking the angle between the
    angular momentum vector and the z-axis, with the angular momentum
    calculated with respect to the center of mass.

    Parameters
    ----------
    snap
        The Snap object.
    origin : optional
        The origin around which to compute the semi-major axis as a
        Quantity like (x, y, z) * au. Default is (0, 0, 0).
    ignore_accreted : optional
        Ignore accreted particles. Default is False.

    Returns
    -------
    Quantity
        The inclination on the particles.
    """
    if ignore_accreted:
        h: Quantity = snap['smoothing_length']
        pos: Quantity = snap['position'][h > 0]
        vel: Quantity = snap['velocity'][h > 0]
    else:
        pos = snap['position']
        vel = snap['velocity']

    origin = snap.translation if snap.translation is not None else ORIGIN
    pos = pos - origin

    _specific_angular_momentum = cross(pos, vel)

    return np.arccos(
        _specific_angular_momentum[:, 2] / norm(_specific_angular_momentum, axis=1)
    )


def keplerian_frequency(
    snap: SnapLike,
    gravitational_parameter: Quantity = None,
    origin: Quantity = None,
    ignore_accreted: bool = False,
) -> Quantity:
    """Calculate the Keplerian orbital frequency.

    The Keplerian orbital frequency of particles around a mass specified
    by gravitational parameter with an optional to specify the position
    of the mass.

    Parameters
    ----------
    snap
        The Snap object.
    gravitational_parameter
        The gravitational parameter (mu = G M) as a Pint quantity.
    origin : optional
        The origin around which to compute the Keplerian frequency as a
        Quantity like (x, y, z) * au. Default is (0, 0, 0).
    ignore_accreted : optional
        Ignore accreted particles. Default is False.

    Returns
    -------
    Quantity
        The Keplerian frequency on the particles.
    """
    if ignore_accreted:
        h: Quantity = snap['smoothing_length']
        pos: Quantity = snap['position'][h > 0]
    else:
        pos = snap['position']

    origin = snap.translation if snap.translation is not None else ORIGIN
    pos = pos - origin

    mu = snap.properties.get('gravitational_parameter')
    if mu is None:
        raise ValueError(
            'must pass in gravitational_parameter or '
            'set on Snap with set_gravitational_parameter'
        )

    radius = norm(pos, axis=1)
    return np.sqrt(mu / radius ** 3)


def semi_major_axis(
    snap: SnapLike,
    gravitational_parameter: Quantity = None,
    origin: Quantity = None,
    ignore_accreted: bool = False,
) -> Quantity:
    """Calculate the semi-major axis.

    The semi-major axis of particles around a mass specified by
    gravitational parameter with an optional to specify the position of
    the mass.

    Parameters
    ----------
    snap
        The Snap object.
    gravitational_parameter
        The gravitational parameter (mu = G M) as a Pint quantity.
    origin : optional
        The origin around which to compute the semi-major axis as a
        Quantity like (x, y, z) * au. Default is (0, 0, 0).
    ignore_accreted : optional
        Ignore accreted particles. Default is False.

    Returns
    -------
    Quantity
        The semi-major axis on the particles.
    """
    if ignore_accreted:
        h: Quantity = snap['smoothing_length']
        pos: Quantity = snap['position'][h > 0]
        vel: Quantity = snap['velocity'][h > 0]
    else:
        pos = snap['position']
        vel = snap['velocity']

    origin = snap.translation if snap.translation is not None else ORIGIN
    pos = pos - origin

    mu = snap.properties.get('gravitational_parameter')
    if mu is None:
        raise ValueError(
            'must pass in gravitational_parameter or '
            'set on Snap with set_gravitational_parameter'
        )

    radius = norm(pos, axis=1)

    _specific_angular_momentum = cross(pos, vel)
    specific_angular_momentum_magnitude = norm(_specific_angular_momentum, axis=1)

    _specific_kinetic_energy = 1 / 2 * norm(vel, axis=1) ** 2
    specific_potential_energy = -mu / radius
    specific_energy = _specific_kinetic_energy + specific_potential_energy

    term = specific_energy * (specific_angular_momentum_magnitude / mu) ** 2

    _eccentricity = np.sqrt(1 + 2 * term)

    return specific_angular_momentum_magnitude ** 2 / (mu * (1 - _eccentricity ** 2))


def stokes_number(
    snap: SnapLike,
    gravitational_parameter: Quantity = None,
    origin: Quantity = None,
    ignore_accreted: bool = False,
) -> Quantity:
    """Calculate the Stokes number.

    Parameters
    ----------
    snap
        The Snap object.
    gravitational_parameter
        The gravitational parameter (mu = G M) as a Pint quantity.
    origin : optional
        The origin around which to compute the Stokes number as a
        Quantity like (x, y, z) * au. Default is (0, 0, 0).
    ignore_accreted : optional
        Ignore accreted particles. Default is False.

    Returns
    -------
    Quantity
        The Stokes number on the particles.
    """
    if ignore_accreted:
        h: Quantity = snap['smoothing_length']
        pos: Quantity = snap['position'][h > 0]
        t_s: Quantity = snap['stopping_time'][h > 0]
    else:
        pos = snap['position']
        t_s = snap['stopping_time']

    origin = snap.translation if snap.translation is not None else ORIGIN
    pos = pos - origin

    mu = snap.properties.get('gravitational_parameter')
    if mu is None:
        raise ValueError(
            'must pass in gravitational_parameter or '
            'set on Snap with set_gravitational_parameter'
        )

    radius = norm(pos, axis=1)
    Omega_k = np.sqrt(mu / radius ** 3)

    Stokes = t_s * Omega_k[:, np.newaxis]
    return Stokes


def normal(snap: SnapLike, ignore_accreted: bool = True) -> ndarray:
    """Calculate unit normal to disc.

    Parameters
    ----------
    snap
        The Snap object.
    ignore_accreted : optional
        Ignore accreted particles. Default is True.

    Returns
    -------
    ndarray
        A unit normal to the plane of the disc.
    """
    origin = center_of_mass(snap=snap, ignore_accreted=ignore_accreted)
    L = angular_momentum(
        snap=snap, origin=origin, ignore_accreted=ignore_accreted
    ).magnitude
    return L / np.linalg.norm(L)


def rotate_face_on(snap: SnapLike, ignore_accreted: bool = True) -> SnapLike:
    """Rotate disc to face-on.

    Parameters
    ----------
    snap
        The Snap object.
    ignore_accreted : optional
        Ignore accreted particles. Default is True.

    Returns
    -------
    Snap
        The rotated Snap.
    """
    x, y, z = normal(snap=snap, ignore_accreted=ignore_accreted)
    axis = (-x, y, 0)
    angle = np.arctan(np.sqrt(x ** 2 + y ** 2) / z)
    return snap.rotate(axis=axis, angle=angle)


def rotate_edge_on(snap: SnapLike, ignore_accreted: bool = True) -> SnapLike:
    """Rotate disc to edge-on.

    Parameters
    ----------
    snap
        The Snap object.
    ignore_accreted : optional
        Ignore accreted particles. Default is True.

    Returns
    -------
    Snap
        The rotated Snap.
    """
    snap = rotate_face_on(snap=snap, ignore_accreted=ignore_accreted)
    return snap.rotate(axis=(1, 0, 0), angle=np.pi / 2)


def position_angle(snap: SnapLike, ignore_accreted: bool = True) -> Quantity:
    """Calculate the disc position angle.

    The position angle is taken from the x-axis in the xy-plane. It
    defines a unit vector around which the snap is inclined.

    Parameters
    ----------
    snap
        The Snap object.
    ignore_accreted : optional
        Ignore accreted particles. Default is True.

    Returns
    -------
    Quantity
        The disc position angle.
    """
    angmom = angular_momentum(snap=snap, ignore_accreted=ignore_accreted)
    if isinstance(angmom, Quantity):
        pi_2 = np.pi / 2 * Quantity('radian')
    else:
        pi_2 = np.pi / 2
    return np.arctan2(angmom[1], angmom[0]) + pi_2


def inclination_disc(snap: SnapLike, ignore_accreted: bool = True) -> Quantity:
    """Calculate the disc inclination.

    The inclination is calculated by taking the angle between the
    angular momentum vector and the z-axis, with the angular momentum
    calculated with respect to the center of mass.

    Parameters
    ----------
    snap
        The Snap object.
    ignore_accreted : optional
        Ignore accreted particles. Default is True.

    Returns
    -------
    Quantity
        The disc inclination.
    """
    angmom = angular_momentum(snap=snap, ignore_accreted=ignore_accreted)
    return np.arccos(angmom[2] / norm(angmom))
