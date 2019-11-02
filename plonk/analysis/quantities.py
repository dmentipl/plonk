"""Analysis quantities.

Quantities include:
- center of mass
- momentum
- angular momentum
"""

from typing import Optional

import numpy as np
from numpy import ndarray

from ..dump.dump import Dump


def center_of_mass(dump: Dump, mask: Optional[ndarray] = None) -> ndarray:
    """Calculate the center of mass on a dump.

    Parameters
    ----------
    dump
        The Dump object.
    mask, optional
        Mask the particle arrays. Default is None.

    Returns
    -------
    ndarray
        The center of mass as a vector (cx, cy, cz).
    """
    mass = dump.mass[:, np.newaxis]
    pos = dump.particles.arrays['xyz'][:]
    if mask is None:
        return (mass * pos).sum(axis=0)
    return (mass * pos)[mask].sum(axis=0)


def momentum(dump: Dump, mask: Optional[ndarray] = None) -> ndarray:
    """Calculate the total momentum vector on a dump.

    Parameters
    ----------
    dump
        The Dump object.
    mask, optional
        Mask the particle arrays. Default is None.

    Returns
    -------
    ndarray
        The total momentum as a vector (px, py, pz).
    """
    mass = dump.mass[:, np.newaxis]
    vel = dump.particles.arrays['vxyz'][:]
    if mask is None:
        return (mass * vel).sum(axis=0)
    return (mass * vel)[mask].sum(axis=0)


def angular_momentum(dump: Dump, mask: Optional[ndarray] = None) -> ndarray:
    """Calculate the total angular momentum vector on a dump.

    Parameters
    ----------
    dump
        The Dump object.
    mask, optional
        Mask the particle arrays. Default is None.

    Returns
    -------
    ndarray
        The total angular momentum as a vector (lx, ly, lz).
    """
    mass = dump.mass[:, np.newaxis]
    pos = dump.particles.arrays['xyz'][:]
    vel = dump.particles.arrays['vxyz'][:]
    if mask is None:
        return (mass * np.cross(pos, vel)).sum(axis=0)
    return (mass * np.cross(pos, vel))[mask].sum(axis=0)
