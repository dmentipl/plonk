"""Analysis quantities.

Quantities include:
- center of mass
- momentum
- angular momentum
- kinetic energy
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
    """Calculate the total momentum vector on a snapshot.

    Parameters
    ----------
    snap
        The Snap object.

    Returns
    -------
    ndarray
        The total momentum as a vector (px, py, pz).
    """
    mass: ndarray = snap['mass']
    vel: ndarray = snap['velocity']
    return (mass[:, np.newaxis] * vel).sum(axis=0)


def angular_momentum(snap: SnapLike) -> ndarray:
    """Calculate the total angular momentum vector on a snapshot.

    Parameters
    ----------
    snap
        The Snap object.
    mask
        Mask the particle arrays. Default is None.

    Returns
    -------
    ndarray
        The total angular momentum as a vector (lx, ly, lz).
    """
    mass: ndarray = snap['mass']
    pos: ndarray = snap['position']
    vel: ndarray = snap['velocity']
    return (mass[:, np.newaxis] * np.cross(pos, vel)).sum(axis=0)


def kinetic_energy(snap: SnapLike) -> ndarray:
    """Calculate the kinetic energy on a snapshot.

    Parameters
    ----------
    snap
        The Snap object.
    mask
        Mask the particle arrays. Default is None.

    Returns
    -------
    ndarray
        The total kinetic energy.
    """
    mass: ndarray = snap['mass']
    vel: ndarray = snap['velocity']
    return (1 / 2 * mass * np.linalg.norm(vel, axis=1) ** 2).sum(axis=0)
