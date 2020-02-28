"""Analysis quantities.

Calculate various quantities on the particles.
"""

from typing import Union

import numpy as np
from numpy import ndarray

from ..snap.snap import Snap, SubSnap

SnapLike = Union[Snap, SubSnap]


def center_of_mass(snap: SnapLike) -> ndarray:
    """Calculate the center of mass on a snapshot.

    Parameters
    ----------
    snap
        The Snap object.

    Returns
    -------
    ndarray
        The center of mass as a vector (cx, cy, cz).
    """
    mass: ndarray = snap['mass']
    pos: ndarray = snap['position']
    return (mass[:, np.newaxis] * pos).sum(axis=0)


def momentum(snap: SnapLike) -> ndarray:
    """Calculate the momentum.

    Parameters
    ----------
    snap
        The Snap object.

    Returns
    -------
    ndarray
        The linear momentum on the particles.
    """
    mass: ndarray = snap['mass']
    vel: ndarray = snap['velocity']
    return mass[:, np.newaxis] * vel


def angular_momentum(snap: SnapLike) -> ndarray:
    """Calculate the angular momentum.

    Parameters
    ----------
    snap
        The Snap object.

    Returns
    -------
    ndarray
        The angular momentum on the particles.
    """
    mass: ndarray = snap['mass']
    pos: ndarray = snap['position']
    vel: ndarray = snap['velocity']
    return mass[:, np.newaxis] * np.cross(pos, vel)


def specific_angular_momentum(snap: SnapLike) -> ndarray:
    """Calculate the specific angular momentum.

    Parameters
    ----------
    snap
        The Snap object.

    Returns
    -------
    ndarray
        The angular momentum on the particles.
    """
    pos: ndarray = snap['position']
    vel: ndarray = snap['velocity']
    return np.cross(pos, vel)


def kinetic_energy(snap: SnapLike) -> ndarray:
    """Calculate the kinetic energy.

    Parameters
    ----------
    snap
        The Snap object.

    Returns
    -------
    ndarray
        The kinetic energy on the particles.
    """
    mass: ndarray = snap['mass']
    vel: ndarray = snap['velocity']
    return 1 / 2 * mass * np.linalg.norm(vel, axis=1) ** 2


def eccentricity(snap: SnapLike, gravitational_parameter: float) -> ndarray:
    """Calculate the eccentricity.

    The eccentricity of particles around some object specified at the
    origin with some gravitational parameter.


    Parameters
    ----------
    snap
        The Snap object.
    gravitational_parameter
        The gravitational parameter (G*M).

    Returns
    -------
    ndarray
        The eccentricity on the particles.
    """
    mu = gravitational_parameter

    pos: ndarray = snap['position']
    vel: ndarray = snap['velocity']

    r = np.sqrt(pos[:, 0] ** 2 + pos[:, 1] ** 2 + pos[:, 2] ** 2)
    v = np.sqrt(vel[:, 0] ** 2 + vel[:, 1] ** 2 + vel[:, 2] ** 2)

    h = np.cross(pos, vel)
    h_mag = np.sqrt(h[:, 0] ** 2 + h[:, 1] ** 2 + h[:, 2] ** 2)

    ke = 0.5 * v ** 2
    pe = -mu / r
    e = ke + pe
    term = 2 * e * h_mag ** 2 / mu ** 2
    ecc = np.sqrt(1 + term)

    return ecc
