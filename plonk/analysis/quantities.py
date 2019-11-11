"""Analysis quantities.

Quantities include:
- center of mass
- momentum
- angular momentum
"""

from typing import Optional

import numpy as np
from numpy import ndarray

from ..snap.snap import Snap


def center_of_mass(snap: Snap, mask: Optional[ndarray] = None) -> ndarray:
    """Calculate the center of mass on a snapshot.

    Parameters
    ----------
    snap
        The Snap object.
    mask
        Mask the particle arrays. Default is None.

    Returns
    -------
    ndarray
        The center of mass as a vector (cx, cy, cz).
    """
    mass: ndarray = snap['mass']
    pos: ndarray = snap['position']
    if mask is None:
        return (mass * pos).sum(axis=0)
    return (mass * pos)[mask].sum(axis=0)


def momentum(snap: Snap, mask: Optional[ndarray] = None) -> ndarray:
    """Calculate the total momentum vector on a snapshot.

    Parameters
    ----------
    snap
        The Snap object.
    mask
        Mask the particle arrays. Default is None.

    Returns
    -------
    ndarray
        The total momentum as a vector (px, py, pz).
    """
    mass: ndarray = snap['mass']
    vel: ndarray = snap['velocity']
    if mask is None:
        return (mass * vel).sum(axis=0)
    return (mass * vel)[mask].sum(axis=0)


def angular_momentum(snap: Snap, mask: Optional[ndarray] = None) -> ndarray:
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
    if mask is None:
        return (mass * np.cross(pos, vel)).sum(axis=0)
    return (mass * np.cross(pos, vel))[mask].sum(axis=0)
