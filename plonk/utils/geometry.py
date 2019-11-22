"""Geometry, and coordinate transformations."""

from typing import Optional, Tuple

import numpy as np
from numpy import ndarray


def coordinate_transform(
    *,
    position: ndarray,
    velocity: ndarray = None,
    geometry_from: str,
    geometry_to: str,
    in_place: bool = False,
) -> Optional[Tuple[ndarray, Optional[ndarray]]]:
    """Coordinate transformation.

    Transform 3d coordinates from one system to another. Coordinate
    systems supported: 'cartesian', 'cylindrical', 'spherical'.
    Performs transformation of positions and, optionally, velocities.

    Parameters
    ----------
    position
        The 3d coordinates to transform, as a (N, 3) ndarray.
    velocity
        The 3d velocity components to transform, as a (N, 3) ndarray.
    geometry_from
        The geometry that the coordinates are in.
    geometry_to
        The geometry to convert to.
    in_place
        If True, the coordinate transformation operates on the the
        array in place, and the function returns None. Default: False.

    Returns
    -------
    position : ndarray
        The coordinates transformed to the new coordinate system, if
        in_place is False. Otherwise, return None.
    velocity : ndarray
        The velocity components transformed to the new coordinate
        system, if in_place is False, and if velocity is not None.
        Otherwise, return None.
    """
    _GEOMETRIES = ('cartesian', 'cylindrical', 'spherical')

    geometry_from = geometry_from.lower()
    geometry_to = geometry_to.lower()
    if geometry_from not in _GEOMETRIES:
        raise ValueError('"geometry_from" not available')
    if geometry_to not in _GEOMETRIES:
        raise ValueError('"geometry_to" not available')

    if geometry_from == 'cartesian':
        if geometry_to == 'cylindrical':
            return _cartesian_to_cylindrical(position, velocity, in_place)
        elif geometry_to == 'spherical':
            return _cartesian_to_spherical(position, velocity, in_place)
    if geometry_from == 'spherical':
        if geometry_to == 'cartesian':
            return _spherical_to_cartesian(position, velocity, in_place)
        else:
            raise ValueError('Can only convert spherical to cartesian')
    if geometry_from == 'cylindrical':
        if geometry_to == 'cartesian':
            return _cylindrical_to_cartesian(position, velocity, in_place)
        else:
            raise ValueError('Can only convert cylindrical to cartesian')
    raise ValueError('Failed to perform coordinate transform')


def _cartesian_to_cylindrical(
    position: ndarray, velocity: ndarray = None, in_place: bool = False
) -> Optional[Tuple[ndarray, Optional[ndarray]]]:
    x, y, z = position[:, 0], position[:, 1], position[:, 2]
    r, phi = np.hypot(x, y), np.arctan2(y, x)
    phi[phi < 0] += 2 * np.pi
    if velocity is not None:
        vx, vy, vz = velocity[:, 0], velocity[:, 1], velocity[:, 2]
        vr = (x * vx + y * vy) / r
        vphi = (x * vy - y * vx) / r ** 2
    if in_place:
        position[:, 0] = r
        position[:, 1] = phi
        if velocity is not None:
            velocity[:, 0] = vr
            velocity[:, 1] = vphi
        return None
    _position = np.zeros(position.shape)
    _position[:, 0] = r
    _position[:, 1] = phi
    _position[:, 2] = z
    if velocity is not None:
        _velocity = np.zeros(velocity.shape)
        _velocity[:, 0] = vr
        _velocity[:, 1] = vphi
        _velocity[:, 2] = vz
        return _position, _velocity
    return _position, None


def _cylindrical_to_cartesian(
    position: ndarray, velocity: ndarray = None, in_place: bool = False
) -> Optional[Tuple[ndarray, Optional[ndarray]]]:
    r, phi, z = position[:, 0], position[:, 1], position[:, 2]
    x = r * np.cos(phi)
    y = r * np.sin(phi)
    if velocity is not None:
        vr, vphi, vz = velocity[:, 0], velocity[:, 1], velocity[:, 2]
        vx = vr * np.cos(phi) - r * vphi * np.sin(phi)
        vy = vr * np.sin(phi) + r * vphi * np.cos(phi)
    if in_place:
        position[:, 0] = x
        position[:, 1] = y
        if velocity is not None:
            velocity[:, 0] = vx
            velocity[:, 1] = vy
        return None
    _position = np.zeros(position.shape)
    _position[:, 0] = x
    _position[:, 1] = y
    _position[:, 2] = z
    if velocity is not None:
        _velocity = np.zeros(velocity.shape)
        _velocity[:, 0] = vx
        _velocity[:, 1] = vy
        _velocity[:, 2] = vz
        return _position, _velocity
    return _position, None


def _cartesian_to_spherical(
    position: ndarray, velocity: ndarray = None, in_place: bool = False
) -> Optional[Tuple[ndarray, Optional[ndarray]]]:
    x, y, z = position[:, 0], position[:, 1], position[:, 2]
    xy = np.hypot(x, y)
    r = np.hypot(xy, z)
    theta = np.arccos(z / r)
    phi = np.arctan2(y, x)
    if velocity is not None:
        vx, vy, vz = velocity[:, 0], velocity[:, 1], velocity[:, 2]
        vr = (x * vx + y * vy + z * vz) / r
        vtheta = (r * vz - vr * z) / r ** 2
        vphi = (x * vy - vx * y) / xy ** 2
    if in_place:
        position[:, 0] = r
        position[:, 1] = theta
        position[:, 2] = phi
        if velocity is not None:
            velocity[:, 0] = vr
            velocity[:, 1] = vtheta
            velocity[:, 2] = vphi
        return None
    _position = np.zeros(position.shape)
    _position[:, 0] = r
    _position[:, 1] = theta
    _position[:, 2] = phi
    if velocity is not None:
        _velocity = np.zeros(position.shape)
        _velocity[:, 0] = vr
        _velocity[:, 1] = vtheta
        _velocity[:, 2] = vphi
        return _position, _velocity
    return _position, None


def _spherical_to_cartesian(
    position: ndarray, velocity: ndarray = None, in_place: bool = False
) -> Optional[Tuple[ndarray, Optional[ndarray]]]:
    r, theta, phi = position[:, 0], position[:, 1], position[:, 2]
    x = r * np.sin(theta) * np.cos(phi)
    y = r * np.sin(theta) * np.sin(phi)
    z = r * np.cos(theta)
    if velocity is not None:
        vr, vtheta, vphi = velocity[:, 0], velocity[:, 1], velocity[:, 2]
        vx = (
            vr * np.sin(theta) * np.cos(phi)
            + r * vtheta * np.cos(theta) * np.cos(phi)
            - r * vphi * np.sin(theta) * np.sin(phi)
        )
        vy = (
            vr * np.sin(theta) * np.sin(phi)
            + r * vtheta * np.cos(theta) * np.sin(phi)
            + r * vphi * np.sin(theta) * np.cos(phi)
        )
        vz = vr * np.cos(theta) - r * vtheta * np.sin(theta)
    if in_place:
        position[:, 0] = x
        position[:, 1] = y
        position[:, 2] = z
        if velocity is not None:
            velocity[:, 0] = vx
            velocity[:, 1] = vy
            velocity[:, 2] = vz
        return None
    _position = np.zeros(position.shape)
    _position[:, 0] = x
    _position[:, 1] = y
    _position[:, 2] = z
    if velocity is not None:
        _velocity = np.zeros(velocity.shape)
        _velocity[:, 0] = vx
        _velocity[:, 1] = vy
        _velocity[:, 2] = vz
        return _position, _velocity
    return _position, None
