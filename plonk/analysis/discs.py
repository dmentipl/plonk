"""Analysis for accretion discs."""

import numpy as np
from numpy import ndarray

from .._units import units as plonk_units
from ..snap.snap import Snap
from .total import angular_momentum, center_of_mass

ORIGIN = (0, 0, 0) * plonk_units.au


def normal(snap: Snap, ignore_accreted: bool = True) -> ndarray:
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


def rotate_face_on(snap: Snap, ignore_accreted: bool = True) -> Snap:
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


def rotate_edge_on(snap: Snap, ignore_accreted: bool = True) -> Snap:
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
