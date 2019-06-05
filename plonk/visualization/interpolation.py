"""
Splash interpolation wrapper routines.

There are two main routines. One for interpolation of scalar fields and
one for interpolation of vector fields. They both make calls to the
Splash Fortran libraries via cython.
"""

import warnings

import numpy as np

try:
    from splash import splash
except ImportError:
    raise Exception(
        'Cannot import Splash interpolation routines. See documentation.'
    )


def scalar_interpolation(
    positions,
    smoothing_length,
    weights,
    scalar_data,
    particle_mass,
    horizontal_range,
    vertical_range,
    *,
    number_pixels=None,
    cross_section=None,
    slice_position=None,
    perspective=None,
    observer_distance=None,
    opacity=None,
    normalize=None,
    accelerate=None,
    integrated_z=None,
    **kwargs,
):
    """
    Interpolate a scalar quantity to a pixel grid.

    Parameters
    ----------
    positions : numpy.ndarray
        Particle positions where columns are Cartesian 'x', 'y', 'z'.
    smoothing_length : numpy.ndarray
        Particle smoothing length.
    weights : numpy.ndarray
        Interpolation weights.
    scalar_data : numpy.ndarray
        A scalar quantity on the particles to interpolate.
    particle_mass : numpy.ndarray
        The particle mass on each particle.
    horizontal_range : list of float
        The horizontal range as [xmin, xmax].
    vertical_range : list of float
        The vertical range as [ymin, ymax].

    Optional parameters
    -------------------
    number_pixels : list of float (default [512, 512])
        The pixel grid to interpolate the scalar quantity to, as
        [npixx, npixy].
    cross_section : bool (default False)
        Turn on cross section rendering.
    slice_position : float (default 0.0)
        Slice position as a z-value.
    perspective : bool (default False)
        Turn on perspective rendering.
    observer_distance : float (default 0.0)
        Distance from the screen to the observer. This turns on 3D
        perspective rendering.
    opacity : bool (default False)
        Turn on opacity rendering.
    normalize : bool (default False)
        Use normalized interpolation.
    accelerate : bool (default False)
        Use accelerated interpolation.
    integrated_z : float (default None)
        Z unit for projection plots.

    Returns
    -------
    numpy.ndarray
        An array of scalar quantities interpolated to a pixel grid with
        shape (npixx, npixy).
    """

    if number_pixels is None:
        npixx = 512
        npixy = 512
    else:
        npixx = number_pixels[0]
        npixy = number_pixels[1]

    projection = False
    if cross_section is None:
        cross_section = False
    if perspective is None:
        perspective = False
    if not cross_section and not perspective:
        projection = True

    if cross_section and slice_position is None:
        slice_position = 0.0
    zslice = slice_position

    if perspective and observer_distance is None:
        observer_distance = 0.0

    dscreen = observer_distance
    zobserver = observer_distance

    if opacity is None:
        opacity = False

    if normalize is None:
        normalise = False
    else:
        normalise = normalize

    if accelerate is None:
        useaccelerate = False
    else:
        useaccelerate = accelerate

    npart = len(smoothing_length)
    itype = np.ones(npart, dtype=np.int32)

    xmin = horizontal_range[0]
    ymin = vertical_range[0]
    xmax = horizontal_range[1]
    ymax = vertical_range[1]
    pixwidthx = (xmax - xmin) / npixx
    pixwidthy = (ymax - ymin) / npixy

    # Splash routines expect single precision
    x = np.array(positions[:, 0], dtype=np.single)
    y = np.array(positions[:, 1], dtype=np.single)
    z = np.array(positions[:, 2], dtype=np.single)
    pmass = np.array(particle_mass, dtype=np.single)
    hh = np.array(smoothing_length, dtype=np.single)
    dat = np.array(scalar_data, dtype=np.single)
    weight = np.array(weights, dtype=np.single)

    if projection or perspective:
        if not opacity:
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
                zobserver=zobserver,
                dscreen=dscreen,
                useaccelerate=useaccelerate,
            )
        else:
            ####################################################################
            # TODO: check this
            warnings.warn('Opacity rendering in Plonk is experimental.')
            npmass = npart
            zorig = z
            pixwidth = pixwidthx
            dscreenfromobserver = dscreen
            rkappa = np.pi * hh.mean() ** 2 / pmass[0]
            zcut = observer_distance
            zobserver = observer_distance
            ####################################################################
            datsmooth = splash.interp3d_proj_opacity(
                x=x,
                y=y,
                z=z,
                hh=hh,
                pmass=pmass,
                npmass=npmass,
                weight=weight,
                dat=dat,
                zorig=zorig,
                itype=itype,
                npart=npart,
                xmin=xmin,
                ymin=ymin,
                npixx=npixx,
                npixy=npixy,
                pixwidth=pixwidth,
                zobserver=zobserver,
                dscreenfromobserver=dscreenfromobserver,
                rkappa=rkappa,
                zcut=zcut,
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
        )

    if integrated_z is not None and projection:
        return integrated_z * np.array(datsmooth)
    return np.array(datsmooth)


def vector_interpolation(
    positions,
    smoothing_length,
    weights,
    vector_data,
    horizontal_range,
    vertical_range,
    *,
    number_pixels=None,
    cross_section=None,
    slice_position=None,
    perspective=None,
    observer_distance=None,
    normalize=None,
    integrated_z=None,
    **kwargs,
):
    """
    Interpolate a vector quantity to a pixel grid.

    Parameters
    ----------
    positions : numpy.ndarray
        Particle positions where columns are Cartesian 'x', 'y', 'z'.
    smoothing_length : numpy.ndarray
        Particle smoothing length.
    weights : numpy.ndarray
        Interpolation weights.
    vector_data : numpy.ndarray
        A vector quantity on the particles to interpolate.
    horizontal_range : list of float
        The horizontal range as [xmin, xmax].
    vertical_range : list of float
        The vertical range as [ymin, ymax].

    Optional parameters
    -------------------
    number_pixels : list of float (default [512, 512])
        The pixel grid to interpolate the scalar quantity to, as
        [npixx, npixy].
    cross_section : bool (default False)
        Turn on cross section rendering.
    slice_position : float (default 0.0)
        Slice position as a z-value.
    perspective : bool (default False)
        Turn on perspective rendering.
    observer_distance : float (default 0.0)
        Distance from the screen to the observer. This turns on 3D
        perspective rendering.
    normalize : bool
        Use normalized interpolation.
    integrated_z : float (default None)
        Z unit for projection plots.

    Returns
    -------
    numpy.ndarray
        An array of vector quantities interpolated to a pixel grid with
        shape (2, npixx, npixy).
    """

    if number_pixels is None:
        npixx = 512
        npixy = 512
    else:
        npixx = number_pixels[0]
        npixy = number_pixels[1]

    projection = False
    if cross_section is None:
        cross_section = False
    if perspective is None:
        perspective = False
    if not cross_section and not perspective:
        projection = True

    if cross_section and slice_position is None:
        slice_position = 0.0
    zslice = slice_position

    if perspective and observer_distance is None:
        observer_distance = 0.0

    dscreen = observer_distance
    zobserver = observer_distance

    if normalize is None:
        normalise = False
    else:
        normalise = normalize

    npart = len(smoothing_length)
    itype = np.ones(npart, dtype=np.int32)

    xmin = horizontal_range[0]
    ymin = vertical_range[0]
    xmax = horizontal_range[1]
    ymax = vertical_range[1]
    pixwidthx = (xmax - xmin) / npixx
    pixwidthy = (ymax - ymin) / npixy

    # Splash routines expect single precision
    x = np.array(positions[:, 0], dtype=np.single)
    y = np.array(positions[:, 1], dtype=np.single)
    z = np.array(positions[:, 2], dtype=np.single)
    hh = np.array(smoothing_length, dtype=np.single)
    vecx = np.array(vector_data[:, 0], dtype=np.single)
    vecy = np.array(vector_data[:, 1], dtype=np.single)
    weight = np.array(weights, dtype=np.single)

    if projection or perspective:
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
            zobserver=zobserver,
            dscreen=dscreen,
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
        )

    if integrated_z is not None and projection:
        return integrated_z * np.stack(
            (np.array(vecsmoothx), np.array(vecsmoothy))
        )
    return np.stack((np.array(vecsmoothx), np.array(vecsmoothy)))
