"""Analysis for accretion discs."""

from __future__ import annotations

from typing import TYPE_CHECKING, Dict, List, Union

import numpy as np
from numpy import ndarray

from .._units import Quantity
from .._units import units as plonk_units
from ..utils.math import cross, norm
from .total import angular_momentum, center_of_mass

if TYPE_CHECKING:
    from ..snap.snap import SnapLike

G = (1 * plonk_units.newtonian_constant_of_gravitation).to_base_units()

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
    central_body: Dict[str, Quantity] = None,
    ignore_accreted: bool = False,
) -> Quantity:
    """Calculate the eccentricity.

    Parameters
    ----------
    snap
        The Snap object.
    central_body : optional
        A dictionary with the mass, position, and velocity (as Pint
        quantities) of the central body around which the particles are
        orbiting. If None, attempt to read from
        snap.properties['central_body'].
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

    if central_body is None:
        try:
            central_body = snap.properties['central_body']
        except KeyError:
            raise ValueError(
                'must pass in central_body or '
                'set both on snap with snap.set_central_body'
            )

    mu = (G * central_body['mass']).to_reduced_units()
    r = pos - central_body['position']
    v = vel - central_body['velocity']

    h = norm(cross(r, v), axis=1)
    eps_k = 1 / 2 * norm(v, axis=1) ** 2
    eps_p = -mu / norm(r, axis=1)
    eps = eps_k + eps_p

    return np.sqrt(1 + 2 * eps * (h / mu) ** 2)


def inclination(
    snap: SnapLike,
    central_body: Dict[str, Quantity] = None,
    ignore_accreted: bool = False,
) -> Quantity:
    """Calculate the inclination.

    Parameters
    ----------
    snap
        The Snap object.
    central_body : optional
        A dictionary with the mass, position, and velocity (as Pint
        quantities) of the central body around which the particles are
        orbiting. If None, attempt to read from
        snap.properties['central_body'].
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

    if central_body is None:
        try:
            central_body = snap.properties['central_body']
        except KeyError:
            raise ValueError(
                'must pass in central_body or '
                'set both on snap with snap.set_central_body'
            )

    r = pos - central_body['position']
    v = vel - central_body['velocity']
    h = cross(r, v)

    return np.arccos(h[:, 2] / norm(h, axis=1))


def keplerian_frequency(
    snap: SnapLike,
    central_body: Dict[str, Quantity] = None,
    ignore_accreted: bool = False,
) -> Quantity:
    """Calculate the Keplerian orbital frequency.

    Parameters
    ----------
    snap
        The Snap object.
    central_body : optional
        A dictionary with the mass, position, and velocity (as Pint
        quantities) of the central body around which the particles are
        orbiting. If None, attempt to read from
        snap.properties['central_body'].
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

    if central_body is None:
        try:
            central_body = snap.properties['central_body']
        except KeyError:
            raise ValueError(
                'must pass in central_body or '
                'set both on snap with snap.set_central_body'
            )

    mu = (G * central_body['mass']).to_reduced_units()
    r = norm(pos - central_body['position'], axis=1)

    return np.sqrt(mu / r ** 3)


def semi_major_axis(
    snap: SnapLike,
    central_body: Dict[str, Quantity] = None,
    ignore_accreted: bool = False,
) -> Quantity:
    """Calculate the semi-major axis.

    Parameters
    ----------
    snap
        The Snap object.
    central_body : optional
        A dictionary with the mass, position, and velocity (as Pint
        quantities) of the central body around which the particles are
        orbiting. If None, attempt to read from
        snap.properties['central_body'].
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

    if central_body is None:
        try:
            central_body = snap.properties['central_body']
        except KeyError:
            raise ValueError(
                'must pass in central_body or '
                'set both on snap with snap.set_central_body'
            )

    mu = (G * central_body['mass']).to_reduced_units()
    r = norm(pos - central_body['position'], axis=1)
    v = vel - central_body['velocity']

    eps_k = 1 / 2 * norm(v, axis=1) ** 2
    eps_p = -mu / r
    eps = eps_k + eps_p

    return -mu / (2 * eps)


def stokes_number(
    snap: SnapLike,
    central_body: Dict[str, Quantity] = None,
    ignore_accreted: bool = False,
) -> Quantity:
    """Calculate the Stokes number.

    Parameters
    ----------
    snap
        The Snap object.
    central_body : optional
        A dictionary with the mass, position, and velocity (as Pint
        quantities) of the central body around which the particles are
        orbiting. If None, attempt to read from
        snap.properties['central_body'].
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

    if central_body is None:
        try:
            central_body = snap.properties['central_body']
        except KeyError:
            raise ValueError(
                'must pass in central_body or '
                'set both on snap with snap.set_central_body'
            )

    mu = (G * central_body['mass']).to_reduced_units()
    r = norm(pos - central_body['position'], axis=1)

    Omega_k = np.sqrt(mu / r ** 3)

    return t_s * Omega_k[:, np.newaxis]


def unit_normal(snap: SnapLike, sinks: Union[bool, List[int]] = False) -> ndarray:
    """Calculate unit normal to plane of rotation.

    I.e. calculate a unit angular momentum vector.

    Parameters
    ----------
    snap
        The Snap object.
    sinks : optional
        Include sink particles specified by a list of indices, or a
        bool indicating all sinks or no sinks. Default is True (all
        sinks).

    Returns
    -------
    ndarray
        A unit angular momentum vector.
    """
    origin = center_of_mass(snap=snap, sinks=sinks)
    L = angular_momentum(snap=snap, sinks=sinks, origin=origin).magnitude
    return L / np.linalg.norm(L)


def rotate_face_on(snap: SnapLike, sinks: Union[bool, List[int]] = False) -> SnapLike:
    """Rotate to face-on with the angular momentum vector.

    I.e. rotate such that the angular momentum points in the
    z-direction.

    Parameters
    ----------
    snap
        The Snap object.
    sinks : optional
        Include sink particles specified by a list of indices, or a
        bool indicating all sinks or no sinks. Default is True (all
        sinks).

    Returns
    -------
    Snap
        The rotated Snap.
    """
    vec = unit_normal(snap=snap, sinks=sinks)
    z_axis = np.array([0, 0, 1])
    axis = cross(vec, z_axis)
    angle = np.arccos(np.dot(vec, z_axis))
    return snap.rotate(axis=axis, angle=angle)


def rotate_edge_on(snap: SnapLike, sinks: Union[bool, List[int]] = False) -> SnapLike:
    """Rotate to edge-on with the angular momentum vector.

    I.e. rotate such that the angular momentum points in the
    y-direction.

    Parameters
    ----------
    snap
        The Snap object.
    sinks : optional
        Include sink particles specified by a list of indices, or a
        bool indicating all sinks or no sinks. Default is True (all
        sinks).

    Returns
    -------
    Snap
        The rotated Snap.
    """
    vec = unit_normal(snap=snap, sinks=sinks)
    y_axis = np.array([0, 1, 0])
    axis = cross(vec, y_axis)
    angle = np.arccos(np.dot(vec, y_axis))
    return snap.rotate(axis=axis, angle=angle)


def position_angle(snap: SnapLike) -> Quantity:
    """Calculate the disc position angle.

    The position angle is taken from the x-axis in the xy-plane. It
    defines a unit vector around which the snap is inclined.

    Parameters
    ----------
    snap
        The Snap object.

    Returns
    -------
    Quantity
        The disc position angle.
    """
    angmom = angular_momentum(snap=snap, sinks=False)
    if isinstance(angmom, Quantity):
        pi_2 = np.pi / 2 * Quantity('radian')
    else:
        pi_2 = np.pi / 2
    return np.arctan2(angmom[1], angmom[0]) + pi_2


def inclination_angle(snap: SnapLike) -> Quantity:
    """Calculate the disc inclination.

    The inclination is calculated by taking the angle between the
    angular momentum vector and the z-axis, with the angular momentum
    calculated with respect to the center of mass.

    Parameters
    ----------
    snap
        The Snap object.

    Returns
    -------
    Quantity
        The disc inclination.
    """
    angmom = angular_momentum(snap=snap, sinks=False)
    return np.arccos(angmom[2] / norm(angmom))
