"""Filter particles to produce SubSnaps."""

from .._units import Quantity
from .._units import units as plonk_units
from ..snap.snap import SnapLike, SubSnap
from .particles import radial_distance

CENTER = (0, 0, 0) * plonk_units.au


def box(
    snap: SnapLike,
    xwidth: Quantity,
    ywidth: Quantity,
    zwidth: Quantity,
    center: Quantity = CENTER,
) -> SubSnap:
    """Particles within a box.

    Parameters
    ----------
    snap
        The Snap object.
    xwidth
        The x-width of the box.
    ywidth
        The y-width of the box.
    zwidth
        The z-width of the box.
    center : optional
        The center of the box as a Quantity like (x, y, z) * au.
        Default is (0, 0, 0).

    Returns
    -------
    SubSnap
        The SubSnap with particles in the box.
    """
    dx, dy, dz = xwidth / 2, ywidth / 2, zwidth / 2
    mask = (
        (snap['position_x'] > center[0] - dx)
        & (snap['position_x'] < center[0] + dx)
        & (snap['position_y'] > center[1] - dy)
        & (snap['position_y'] < center[1] + dy)
        & (snap['position_z'] > center[2] - dz)
        & (snap['position_z'] < center[2] + dz)
    )
    return snap[mask]


def cylinder(
    snap: SnapLike, radius: Quantity, height: Quantity, center: Quantity = CENTER
) -> SubSnap:
    """Particles within a cylinder.

    Parameters
    ----------
    snap
        The Snap object.
    radius
        The radius of the cylinder.
    height
        The height of the cylinder.
    center : optional
        The center of the cylinder as a Quantity like (x, y, z) * au.
        Default is (0, 0, 0).

    Returns
    -------
    SubSnap
        The SubSnap with particles in the cylinder.
    """
    dh = height / 2
    R = radial_distance(snap=snap, origin=center, coordinates='cylindrical')
    mask = (
        (R < radius)
        & (snap['position_z'] < center[2] + dh)
        & (snap['position_z'] > center[2] - dh)
    )
    return snap[mask]


def annulus(
    snap: SnapLike,
    radius_min: Quantity,
    radius_max: Quantity,
    height: Quantity,
    center: Quantity = CENTER,
) -> SubSnap:
    """Particles within an annulus.

    Parameters
    ----------
    snap
        The Snap object.
    radius_min
        The inner radius of the annulus.
    radius_max
        The outer radius of the annulus.
    height
        The height of the annulus.
    center : optional
        The center of the annulus as a Quantity like (x, y, z) * au.
        Default is (0, 0, 0).

    Returns
    -------
    SubSnap
        The SubSnap with particles in the annulus.
    """
    dh = height / 2
    R = radial_distance(snap=snap, origin=center, coordinates='cylindrical')
    mask = (
        (R > radius_min)
        & (R < radius_max)
        & (snap['position_z'] < center[2] + dh)
        & (snap['position_z'] > center[2] - dh)
    )
    return snap[mask]


def sphere(snap: SnapLike, radius: Quantity, center: Quantity = CENTER) -> SubSnap:
    """Particles within a sphere.

    Parameters
    ----------
    snap
        The Snap object.
    radius
        The radius of the sphere.
    center : optional
        The center of the sphere as a Quantity like (x, y, z) * au.
        Default is (0, 0, 0).

    Returns
    -------
    SubSnap
        The SubSnap with particles in the sphere.
    """
    R = radial_distance(snap=snap, origin=center, coordinates='spherical')
    mask = R < radius
    return snap[mask]


def shell(
    snap: SnapLike,
    radius_min: Quantity,
    radius_max: Quantity,
    center: Quantity = CENTER,
) -> SubSnap:
    """Particles within a spherical shell.

    Parameters
    ----------
    snap
        The Snap object.
    radius_min
        The inner radius of the shell.
    radius_max
        The outer radius of the shell.
    center : optional
        The center of the shell as a Quantity like (x, y, z) * au.
        Default is (0, 0, 0).

    Returns
    -------
    SubSnap
        The SubSnap with particles in the shell.
    """
    R = radial_distance(snap=snap, origin=center, coordinates='spherical')
    mask = (R > radius_min) & (R < radius_max)
    return snap[mask]
