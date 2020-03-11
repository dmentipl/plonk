"""Interpolation to a pixel grid.

There are two functions: one for interpolation of scalar fields, and one
for interpolation of vector fields.
"""

from typing import Optional, Tuple, Union

import numpy as np
from numpy import ndarray

from ..snap import SnapLike
from ..snap.snap import get_array_from_input, get_array_in_code_units
from .splash import interpolate_cross_section, interpolate_projection

Extent = Tuple[float, float, float, float]


def interpolate(
    *,
    snap: SnapLike,
    quantity: Union[str, ndarray],
    x: Union[str, ndarray] = 'x',
    y: Union[str, ndarray] = 'y',
    z: Optional[Union[str, ndarray]] = None,
    interp: 'str',
    z_slice: Optional[float] = None,
    extent: Extent,
    **kwargs,
) -> ndarray:
    """Interpolate a quantity on the snapshot to a pixel grid.

    Parameters
    ----------
    snap
        The Snap (or SubSnap) object.
    quantity
        The quantity to visualize. Can be a string to pass to Snap, or
        a 1d array (N,) of scalar data, or a 2d array (N, 3) of
        vector data. If quantity is 2d, only the first two components
        are interpolated, i.e. quantity[:, 0] and quantity[:, 1].
        Default is None.
    x
        The x-coordinate for the visualization. Can be a string to
        pass to Snap, or a 1d array (N,). Default is 'x'.
    y
        The y-coordinate for the visualization. Can be a string to
        pass to Snap, or a 1d array (N,). Default is 'y'.
    z
        The z-coordinate for the visualization. Can be a string to
        pass to Snap, or a 1d array (N,). This is only required for
        cross-section plots. Default is 'z'.
    interp
        The interpolation type. Default is 'projection'.

        - 'projection' : 2d interpolation via projection to xy-plane
        - 'cross_section' : 3d interpolation via cross-section in
          z-direction
    z_slice
        The z-coordinate value of the cross-section slice. Default
        is 0.0.
    extent
        The xy extent of the image as (xmin, xmax, ymin, ymax).
    **kwargs
        Additional keyword arguments to pass to scalar_interpolation
        and vector_interpolation.

    Returns
    -------
    ndarray
        The interpolated quantity on a pixel grid as an ndarray.
    """
    quantity = get_array_from_input(snap, quantity)
    x = get_array_from_input(snap, x, 'x')
    y = get_array_from_input(snap, y, 'y')
    z = get_array_from_input(snap, z, 'z')
    h = get_array_in_code_units(snap, 'smooth')
    m = get_array_in_code_units(snap, 'mass')

    if interp == 'projection':
        cross_section = None
    elif interp == 'cross_section':
        if z_slice is None:
            z_slice = 0.0
        cross_section = z_slice

    if quantity.ndim == 1:
        interpolated_data = scalar_interpolation(
            quantity=quantity,
            x_coordinate=x,
            y_coordinate=y,
            z_coordinate=z,
            extent=extent,
            smoothing_length=h,
            particle_mass=m,
            hfact=snap.properties['hfact'],
            cross_section=cross_section,
            **kwargs,
        )

    elif quantity.ndim == 2:
        interpolated_data = vector_interpolation(
            quantity_x=quantity[:, 0],
            quantity_y=quantity[:, 1],
            x_coordinate=x,
            y_coordinate=y,
            z_coordinate=z,
            extent=extent,
            smoothing_length=h,
            particle_mass=m,
            hfact=snap.properties['hfact'],
            cross_section=cross_section,
            **kwargs,
        )

    else:
        raise ValueError('quantity.ndim > 2: cannot determine quantity')

    return interpolated_data


def scalar_interpolation(
    *,
    quantity: ndarray,
    x_coordinate: ndarray,
    y_coordinate: ndarray,
    z_coordinate: Optional[ndarray] = None,
    extent: Extent,
    smoothing_length: ndarray,
    particle_mass: ndarray,
    hfact: float,
    number_of_pixels: Tuple[float, float] = (512, 512),
    cross_section: Optional[float] = None,
    density_weighted: Optional[bool] = None,
) -> ndarray:
    """Interpolate scalar quantity to a pixel grid.

    Parameters
    ----------
    quantity
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
    number_of_pixels
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
        quantity=quantity,
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
    quantity_x: ndarray,
    quantity_y: ndarray,
    x_coordinate: ndarray,
    y_coordinate: ndarray,
    z_coordinate: Optional[ndarray] = None,
    extent: Extent,
    smoothing_length: ndarray,
    particle_mass: ndarray,
    hfact: float,
    number_of_pixels: Tuple[float, float] = (512, 512),
    cross_section: Optional[float] = None,
    density_weighted: Optional[bool] = None,
) -> ndarray:
    """Interpolate scalar quantity to a pixel grid.

    Parameters
    ----------
    quantity_x
        The x-component of a vector quantity to interpolate.
    quantity_y
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
    number_of_pixels
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
        quantity=quantity_x,
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
        quantity=quantity_y,
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
    quantity: ndarray,
    x_coordinate: ndarray,
    y_coordinate: ndarray,
    z_coordinate: Optional[ndarray] = None,
    extent: Extent,
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
        if z_coordinate is None:
            raise ValueError('Cross section interpolation requires z_coordinate')
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
        interpolated_data = interpolate_cross_section(
            x=x_coordinate,
            y=y_coordinate,
            z=z_coordinate,
            hh=smoothing_length,
            weight=weight,
            dat=quantity,
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
        interpolated_data = interpolate_projection(
            x=x_coordinate,
            y=y_coordinate,
            hh=smoothing_length,
            weight=weight,
            dat=quantity,
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

    return interpolated_data
