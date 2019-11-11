"""Interpolation to a pixel grid.

There are two functions: one for interpolation of scalar fields, and one
for interpolation of vector fields. They both make use of KDEpy.
"""

from typing import Optional, Tuple

import numpy as np
from KDEpy import FFTKDE
from numpy import ndarray
from scipy.interpolate import RectBivariateSpline

_H_FACT = 1.2
_C_NORM_3D = 1 / np.pi


def scalar_interpolation(
    *,
    data: ndarray,
    x_position: ndarray,
    y_position: ndarray,
    z_position: Optional[ndarray] = None,
    extent: Tuple[float, float, float, float],
    smoothing_length: ndarray,
    particle_mass: ndarray,
    number_of_pixels: Tuple[float, float],
    cross_section: Optional[float] = None,
    density_weighted: Optional[bool] = None,
) -> ndarray:
    """Interpolate scalar quantity to a pixel grid.

    Parameters
    ----------
    data
        A scalar quantity on the particles to interpolate.
    x_coordinate
        Particle coordinate for x-axis in interpolation.
    y_coordinate
        Particle coordinate for y-axis in interpolation.
    z_coordinate
        Particle coordinate for z-axis. Only required for cross section
        interpolation.
    extent
        The range in the x- and y-direction as (xmin, xmax, ymin, ymax).
    smoothing_length
        The smoothing length on each particle.
    particle_mass
        The particle mass on each particle.
    number_pixels
        The pixel grid to interpolate the scalar quantity to, as
        (npixx, npixy).
    cross_section
        Cross section slice position as a z-value. If None, cross
        section interpolation is turned off. Default is off.
    density_weighted
        Use density weighted interpolation. Default is off.

    Returns
    -------
    ndarray
        An array of scalar quantities interpolated to a pixel grid with
        shape (npixx, npixy).
    """
    return _interpolate(
        data=data,
        x_position=x_position,
        y_position=y_position,
        z_position=z_position,
        extent=extent,
        smoothing_length=smoothing_length,
        particle_mass=particle_mass,
        number_of_pixels=number_of_pixels,
        cross_section=cross_section,
        density_weighted=density_weighted,
    )


def vector_interpolation(
    *,
    x_data: ndarray,
    y_data: ndarray,
    x_position: ndarray,
    y_position: ndarray,
    z_position: Optional[ndarray] = None,
    extent: Tuple[float, float, float, float],
    smoothing_length: ndarray,
    particle_mass: ndarray,
    number_of_pixels: Tuple[float, float],
    cross_section: Optional[float] = None,
    density_weighted: Optional[bool] = None,
) -> ndarray:
    """Interpolate scalar quantity to a pixel grid.

    Parameters
    ----------
    x_data
        The x-component of a vector quantity to interpolate.
    y_data
        The y-component of a vector quantity to interpolate.
    x_coordinate
        Particle coordinate for x-axis in interpolation.
    y_coordinate
        Particle coordinate for y-axis in interpolation.
    z_coordinate
        Particle coordinate for z-axis. Only required for cross section
        interpolation.
    extent
        The range in the x- and y-direction as (xmin, xmax, ymin, ymax).
    smoothing_length
        The smoothing length on each particle.
    particle_mass
        The particle mass on each particle.
    number_pixels
        The pixel grid to interpolate the scalar quantity to, as
        (npixx, npixy).
    cross_section
        Cross section slice position as a z-value. If None, cross
        section interpolation is turned off. Default is off.
    density_weighted
        Use density weighted interpolation. Default is off.

    Returns
    -------
    ndarray
        An array of vector quantities interpolated to a pixel grid with
        shape (2, npixx, npixy).
    """
    vecsmoothx = _interpolate(
        data=x_data,
        x_position=x_position,
        y_position=y_position,
        z_position=z_position,
        extent=extent,
        smoothing_length=smoothing_length,
        particle_mass=particle_mass,
        number_of_pixels=number_of_pixels,
        cross_section=cross_section,
        density_weighted=density_weighted,
    )
    vecsmoothy = _interpolate(
        data=y_data,
        x_position=x_position,
        y_position=y_position,
        z_position=z_position,
        extent=extent,
        smoothing_length=smoothing_length,
        particle_mass=particle_mass,
        number_of_pixels=number_of_pixels,
        cross_section=cross_section,
        density_weighted=density_weighted,
    )
    return np.stack((np.array(vecsmoothx), np.array(vecsmoothy)))


def _interpolate(
    *,
    data: ndarray,
    x_position: ndarray,
    y_position: ndarray,
    z_position: Optional[ndarray] = None,
    extent: Tuple[float, float, float, float],
    smoothing_length: ndarray,
    particle_mass: ndarray,
    number_of_pixels: Tuple[float, float],
    cross_section: Optional[float] = None,
    density_weighted: Optional[bool] = None,
) -> ndarray:
    normalized = False
    if density_weighted is None:
        density_weighted = False
    if density_weighted:
        normalized = True

    mask = smoothing_length > 0.0
    mask = mask & (
        (x_position >= extent[0])
        & (x_position <= extent[1])
        & (y_position >= extent[2])
        & (y_position <= extent[3])
    )
    if cross_section is not None:
        if z_position is None:
            raise ValueError('Must specify z position for cross section')
        mask = mask & (np.abs(z_position - cross_section) < 2 * smoothing_length)

    xy = np.vstack((x_position[mask], y_position[mask])).T
    scalar = data[mask]
    h = smoothing_length[mask]
    m = particle_mass[mask]

    if density_weighted:
        if cross_section is not None:
            weights = scalar * m / h * _C_NORM_3D
        else:
            weights = scalar * m
    else:
        if cross_section is not None:
            weights = scalar * h ** 2 * _C_NORM_3D / _H_FACT ** 3
        else:
            weights = scalar * h ** 3 / _H_FACT ** 3
    if normalized:
        weights_norm = weights / scalar

    kde = FFTKDE(kernel='gaussian')
    grid, points = kde.fit(xy, weights=weights).evaluate(number_of_pixels)
    z = points.reshape(number_of_pixels)

    if normalized:
        _, points_norm = kde.fit(xy, weights=weights_norm).evaluate(number_of_pixels)
        z_norm = points_norm.reshape(number_of_pixels)
        z /= z_norm

    normalization = np.sum(weights)
    if normalized:
        normalization /= np.sum(m)
    z *= normalization

    x_grid = np.linspace(grid[0, 0], grid[-1, 0], number_of_pixels[0])
    y_grid = np.linspace(grid[0, 1], grid[-1, 1], number_of_pixels[1])
    spl = RectBivariateSpline(x_grid, y_grid, z)
    x_regrid = np.linspace(*extent[:2], number_of_pixels[0])
    y_regrid = np.linspace(*extent[2:], number_of_pixels[1])
    z_regrid = spl(x_regrid, y_regrid)

    return z_regrid.T
