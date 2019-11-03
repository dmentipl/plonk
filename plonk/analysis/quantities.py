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
    mask, optional
        Mask the particle arrays. Default is None.

    Returns
    -------
    ndarray
        The center of mass as a vector (cx, cy, cz).
    """
    mass = snap.mass[:, np.newaxis]
    pos = snap.particles.arrays['xyz'][:]
    if mask is None:
        return (mass * pos).sum(axis=0)
    return (mass * pos)[mask].sum(axis=0)


def momentum(snap: Snap, mask: Optional[ndarray] = None) -> ndarray:
    """Calculate the total momentum vector on a snapshot.

    Parameters
    ----------
    snap
        The Snap object.
    mask, optional
        Mask the particle arrays. Default is None.

    Returns
    -------
    ndarray
        The total momentum as a vector (px, py, pz).
    """
    mass = snap.mass[:, np.newaxis]
    vel = snap.particles.arrays['vxyz'][:]
    if mask is None:
        return (mass * vel).sum(axis=0)
    return (mass * vel)[mask].sum(axis=0)


def angular_momentum(snap: Snap, mask: Optional[ndarray] = None) -> ndarray:
    """Calculate the total angular momentum vector on a snapshot.

    Parameters
    ----------
    snap
        The Snap object.
    mask, optional
        Mask the particle arrays. Default is None.

    Returns
    -------
    ndarray
        The total angular momentum as a vector (lx, ly, lz).
    """
    mass = snap.mass[:, np.newaxis]
    pos = snap.particles.arrays['xyz'][:]
    vel = snap.particles.arrays['vxyz'][:]
    if mask is None:
        return (mass * np.cross(pos, vel)).sum(axis=0)
    return (mass * np.cross(pos, vel))[mask].sum(axis=0)
