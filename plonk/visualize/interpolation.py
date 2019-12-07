"""Interpolation to a pixel grid.

There are two functions: one for interpolation of scalar fields, and one
for interpolation of vector fields.
"""

from typing import Optional, Tuple

import numpy as np
from numpy import ndarray

from .splash import interpolate_cross_section, interpolate_projection


def scalar_interpolation(
    *,
    data: ndarray,
    x_coordinate: ndarray,
    y_coordinate: ndarray,
    z_coordinate: Optional[ndarray] = None,
    extent: Tuple[float, float, float, float],
    smoothing_length: ndarray,
    particle_mass: ndarray,
    hfact: float,
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
    hfact
        The smoothing length factor.
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
        x_coordinate=x_coordinate,
        y_coordinate=y_coordinate,
        z_coordinate=z_coordinate,
        extent=extent,
        smoothing_length=smoothing_length,
        particle_mass=particle_mass,
        hfact=hfact,
        number_of_pixels=number_of_pixels,
        cross_section=cross_section,
        density_weighted=density_weighted,
    )


def vector_interpolation(
    *,
    x_data: ndarray,
    y_data: ndarray,
    x_coordinate: ndarray,
    y_coordinate: ndarray,
    z_coordinate: Optional[ndarray] = None,
    extent: Tuple[float, float, float, float],
    smoothing_length: ndarray,
    particle_mass: ndarray,
    hfact: float,
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
    hfact
        The smoothing length factor.
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
        x_coordinate=x_coordinate,
        y_coordinate=y_coordinate,
        z_coordinate=z_coordinate,
        extent=extent,
        smoothing_length=smoothing_length,
        particle_mass=particle_mass,
        hfact=hfact,
        number_of_pixels=number_of_pixels,
        cross_section=cross_section,
        density_weighted=density_weighted,
    )
    vecsmoothy = _interpolate(
        data=y_data,
        x_coordinate=x_coordinate,
        y_coordinate=y_coordinate,
        z_coordinate=z_coordinate,
        extent=extent,
        smoothing_length=smoothing_length,
        particle_mass=particle_mass,
        hfact=hfact,
        number_of_pixels=number_of_pixels,
        cross_section=cross_section,
        density_weighted=density_weighted,
    )
    return np.stack((np.array(vecsmoothx), np.array(vecsmoothy)))


def _interpolate(
    *,
    data: ndarray,
    x_coordinate: ndarray,
    y_coordinate: ndarray,
    z_coordinate: Optional[ndarray] = None,
    extent: Tuple[float, float, float, float],
    smoothing_length: ndarray,
    particle_mass: ndarray,
    hfact: float,
    number_of_pixels: Tuple[float, float],
    cross_section: Optional[float] = None,
    density_weighted: Optional[bool] = None,
) -> ndarray:
    if cross_section is None:
        do_cross_section = False
    else:
        do_cross_section = True
        zslice = cross_section
    normalise = False
    if density_weighted is None:
        density_weighted = False
    if density_weighted:
        normalise = True

    npixx, npixy = number_of_pixels
    xmin, ymin = extent[0], extent[2]
    pixwidthx = (extent[1] - extent[0]) / npixx
    pixwidthy = (extent[3] - extent[2]) / npixy
    npart = len(smoothing_length)

    itype = np.ones(smoothing_length.shape)
    if density_weighted:
        weight = particle_mass / smoothing_length ** 3
    else:
        weight = hfact ** -3 * np.ones(smoothing_length.shape)

    if do_cross_section:
        data = interpolate_cross_section(
            x=x_coordinate,
            y=y_coordinate,
            z=z_coordinate,
            hh=smoothing_length,
            weight=weight,
            dat=data,
            itype=itype,
            npart=npart,
            xmin=xmin,
            ymin=ymin,
            zslice=zslice,
            npixx=npixx,
            npixy=npixy,
            pixwidthx=pixwidthx,
            pixwidthy=pixwidthy,
            normalise=normalise,
        )
    else:
        data = interpolate_projection(
            x=x_coordinate,
            y=y_coordinate,
            z=z_coordinate,
            hh=smoothing_length,
            weight=weight,
            dat=data,
            itype=itype,
            npart=npart,
            xmin=xmin,
            ymin=ymin,
            npixx=npixx,
            npixy=npixy,
            pixwidthx=pixwidthx,
            pixwidthy=pixwidthy,
            normalise=normalise,
        )

    return data
