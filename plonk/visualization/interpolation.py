"""
Splash interpolation wrapper routines.

There are two main routines. One for interpolation of scalar fields and
one for interpolation of vector fields. They both make calls to the
Splash Fortran libraries via Cython.
"""

from typing import Tuple

import numpy as np
from numpy import ndarray

try:
    from splash import splash
except ImportError:
    raise Exception('Cannot import Splash interpolation routines. See documentation.')

_IVERBOSE = -2
_PIXEL = (512, 512)


def scalar_interpolation(
    x_coordinate: ndarray,
    y_coordinate: ndarray,
    z_coordinate: ndarray,
    smoothing_length: ndarray,
    smoothing_length_factor: float,
    scalar_data: ndarray,
    particle_mass: ndarray,
    x_range: Tuple[float, float],
    y_range: Tuple[float, float],
    *,
    number_pixels: Tuple[float, float] = _PIXEL,
    cross_section: bool = False,
    slice_position: float = 0.0,
    density_weighted: bool = False,
    normalize: bool = False,
) -> ndarray:
    """
    Interpolate a scalar quantity to a pixel grid.

    Parameters
    ----------
    x_coordinate
        Particle coordinate for x-axis in interpolation.
    y_coordinate
        Particle coordinate for y-axis in interpolation.
    z_coordinate
        Particle coordinate for depth in interpolation.
    smoothing_length
        Particle smoothing length.
    smoothing_length_factor
        The smoothing length factor.
    scalar_data
        A scalar quantity on the particles to interpolate.
    particle_mass
        The particle mass on each particle.
    x_range
        The range in the x-direction as (xmin, xmax).
    y_range
        The vertical range as (ymin, ymax).

    Optional parameters
    -------------------
    number_pixels
        The pixel grid to interpolate the scalar quantity to, as
        (npixx, npixy).
    cross_section
        Turn on cross section rendering.
    slice_position
        Slice position as a z-value. Default is 0.0.
    density_weighted
        Use density weighted interpolation.
    normalize
        Use normalized interpolation.

    Returns
    -------
    ndarray
        An array of scalar quantities interpolated to a pixel grid with
        shape (npixx, npixy).
    """
    npixx = number_pixels[0]
    npixy = number_pixels[1]
    projection = not cross_section
    zslice = slice_position
    normalise = normalize
    npart = len(smoothing_length)
    xmin = x_range[0]
    ymin = y_range[0]
    xmax = x_range[1]
    ymax = y_range[1]
    pixwidthx = (xmax - xmin) / npixx
    pixwidthy = (ymax - ymin) / npixy
    if density_weighted:
        weights = particle_mass / smoothing_length ** 2
    else:
        weights = np.ones(npart) / smoothing_length_factor

    # Splash routines expect single precision
    x = np.array(x_coordinate, dtype=np.single)
    y = np.array(y_coordinate, dtype=np.single)
    z = np.array(z_coordinate, dtype=np.single)
    hh = np.array(smoothing_length, dtype=np.single)
    weight = np.array(weights, dtype=np.single)
    dat = np.array(scalar_data, dtype=np.single)
    itype = np.ones(npart, dtype=np.int32)

    if projection:
        datsmooth = splash.interpolate3d_projection(
            x=x,
            y=y,
            z=z,
            hh=hh,
            weight=weight,
            dat=dat,
            itype=itype,
            npart=npart,
            xmin=xmin,
            ymin=ymin,
            npixx=npixx,
            npixy=npixy,
            pixwidthx=pixwidthx,
            pixwidthy=pixwidthy,
            normalise=normalise,
            zobserver=0.0,
            dscreen=0.0,
            useaccelerate=False,
            iverbose=_IVERBOSE,
        )
    elif cross_section:
        datsmooth = splash.interpolate3d_fastxsec(
            x=x,
            y=y,
            z=z,
            hh=hh,
            weight=weight,
            dat=dat,
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
            iverbose=_IVERBOSE,
        )

    return np.array(datsmooth)


def vector_interpolation(
    x_coordinate: ndarray,
    y_coordinate: ndarray,
    z_coordinate: ndarray,
    smoothing_length: ndarray,
    smoothing_length_factor: float,
    x_vector_data: ndarray,
    y_vector_data: ndarray,
    particle_mass: ndarray,
    x_range: Tuple[float, float],
    y_range: Tuple[float, float],
    *,
    number_pixels: Tuple[float, float] = _PIXEL,
    cross_section: bool = False,
    slice_position: float = 0.0,
    density_weighted: bool = False,
    normalize: bool = False,
) -> ndarray:
    """
    Interpolate a vector quantity to a pixel grid.

    Parameters
    ----------
    x_coordinate
        Particle coordinate for x-axis in interpolation.
    y_coordinate
        Particle coordinate for y-axis in interpolation.
    z_coordinate
        Particle coordinate for depth in interpolation.
    smoothing_length
        Particle smoothing length.
    smoothing_length_factor
        The smoothing length factor.
    x_vector_data
        The x-component of a vector quantity on the particles to
        interpolate. This is the x-component with reference to the
        coordinate system.
    y_vector_data
        The y-component of a vector quantity on the particles to
        interpolate. This is the y-component with reference to the
        coordinate system.
    particle_mass
        The particle mass on each particle.
    x_range
        The range in the x-direction as (xmin, xmax).
    y_range
        The vertical range as (ymin, ymax).

    Optional parameters
    -------------------
    number_pixels
        The pixel grid to interpolate the scalar quantity to, as
        (npixx, npixy).
    cross_section
        Turn on cross section rendering.
    slice_position
        Slice position as a z-value.
    density_weighted
        Use density weighted interpolation.
    normalize
        Use normalized interpolation.

    Returns
    -------
    ndarray
        An array of vector quantities interpolated to a pixel grid with
        shape (2, npixx, npixy).
    """
    npixx = number_pixels[0]
    npixy = number_pixels[1]
    projection = not cross_section
    zslice = slice_position
    normalise = normalize
    npart = len(smoothing_length)
    xmin = x_range[0]
    ymin = y_range[0]
    xmax = x_range[1]
    ymax = y_range[1]
    pixwidthx = (xmax - xmin) / npixx
    pixwidthy = (ymax - ymin) / npixy
    if density_weighted:
        weights = particle_mass / smoothing_length ** 2
    else:
        weights = np.ones(npart) / smoothing_length_factor

    # Splash routines expect single precision
    x = np.array(x_coordinate, dtype=np.single)
    y = np.array(y_coordinate, dtype=np.single)
    z = np.array(z_coordinate, dtype=np.single)
    hh = np.array(smoothing_length, dtype=np.single)
    vecx = np.array(x_vector_data, dtype=np.single)
    vecy = np.array(y_vector_data, dtype=np.single)
    weight = np.array(weights, dtype=np.single)
    itype = np.ones(npart, dtype=np.int32)

    if projection:
        vecsmoothx, vecsmoothy = splash.interpolate3d_proj_vec(
            x=x,
            y=y,
            z=z,
            hh=hh,
            weight=weight,
            vecx=vecx,
            vecy=vecy,
            itype=itype,
            npart=npart,
            xmin=xmin,
            ymin=ymin,
            npixx=npixx,
            npixy=npixy,
            pixwidthx=pixwidthx,
            pixwidthy=pixwidthy,
            normalise=normalise,
            zobserver=0.0,
            dscreen=0.0,
            iverbose=_IVERBOSE,
        )
    elif cross_section:
        vecsmoothx, vecsmoothy = splash.interpolate3d_xsec_vec(
            x=x,
            y=y,
            z=z,
            hh=hh,
            weight=weight,
            vecx=vecx,
            vecy=vecy,
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
            iverbose=_IVERBOSE,
        )

    return np.stack((np.array(vecsmoothx), np.array(vecsmoothy)))
