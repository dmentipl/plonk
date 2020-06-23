"""Splash interpolation functions.

Derived from Splash: https://github.com/danieljprice/splash.
"""

import numba
import numpy as np
from numpy import ndarray

RADKERNEL = 2.0
RADKERNEL2 = 4.0
CNORMK3D = 1.0 / np.pi

NPTS = 100
MAXCOLTABLE = 1000
DQ2TABLE = RADKERNEL2 / MAXCOLTABLE
DDQ2TABLE = 1.0 / DQ2TABLE

IVERBOSE = -1


@numba.njit
def w_cubic(q2: float):
    """Cubic spline kernel.

    Parameters
    ----------
    q2

    Returns
    -------
    w
    """
    w = 0.0
    if q2 < 1.0:
        q = np.sqrt(q2)
        w = 1.0 - 1.5 * q2 + 0.75 * q2 * q
    elif q2 < 4:
        q = np.sqrt(q2)
        w = 0.25 * (2.0 - q) ** 3
    return w


@numba.njit
def setup_integratedkernel():
    """Set up integrated kernel.

    Tabulates the integral through the cubic spline kernel tabulated in
    (r/h)**2 so that sqrt is not necessary.

    Returns
    -------
    coltable
    """
    coltable = np.zeros(MAXCOLTABLE + 1)

    for idx in range(MAXCOLTABLE):
        # Tabulate for (cylindrical) r**2 between 0 and RADKERNEL**2
        rxy2 = idx * DQ2TABLE

        # Integrate z between 0 and sqrt(RADKERNEL^2 - rxy^2)
        deltaz = np.sqrt(RADKERNEL2 - rxy2)
        dz = deltaz / (NPTS - 1)
        coldens = 0
        for j in range(NPTS):
            z = j * dz
            q2 = rxy2 + z * z
            wkern = w_cubic(q2)
            if j in (0, NPTS - 1):
                coldens = coldens + 0.5 * wkern * dz
            else:
                coldens = coldens + wkern * dz
        coltable[idx] = 2.0 * coldens * CNORMK3D
    coltable[MAXCOLTABLE] = 0.0

    return coltable


@numba.njit
def wfromtable(q2, coltable):
    """Interpolate from integrated kernel table values to give w(q).

    Parameters
    ----------
    q2
    coltable

    Returns
    -------
    w
    """
    # Find nearest index in table
    index = int(q2 * DDQ2TABLE)
    index1 = min(index, MAXCOLTABLE)

    # Find increment along from this index
    dxx = q2 - index * DQ2TABLE

    # Find gradient
    dwdx = (coltable[index1] - coltable[index]) * DDQ2TABLE

    # Compute value of integrated kernel
    return coltable[index] + dwdx * dxx


@numba.njit
def interpolate_projection(
    x: ndarray,
    y: ndarray,
    hh: ndarray,
    weight: ndarray,
    dat: ndarray,
    itype: ndarray,
    npart: int,
    xmin: float,
    ymin: float,
    npixx: int,
    npixy: int,
    pixwidthx: float,
    pixwidthy: float,
    normalise: bool,
):
    """Interpolate particles to grid via projection.

    Parameters
    ----------
    x
        The particle x positions.
    y
        The particle y positions.
    hh
        The particle smoothing length.
    weight
        The particle weight.
    dat
        The scalar data to interpolate.
    itype
        The particle type.
    npart
        The number of particles.
    xmin
        The minimum x position.
    ymin
        The minimum y position.
    npixx
        The number of pixels in the x direction.
    npixy
        The number of pixels in the y direction.
    normalise
        Whether to normalize.

    Return
    ------
    datsmooth
        The data smoothed to a pixel grid.
    """
    coltable = setup_integratedkernel()

    datsmooth = np.zeros((npixx, npixy))
    datnorm = np.zeros((npixx, npixy))
    dx2i = np.zeros(npixx)
    xpix = np.zeros(npixx)
    term = 0.0

    xminpix = xmin - 0.5 * pixwidthx
    yminpix = ymin - 0.5 * pixwidthy
    xmax = xmin + npixx * pixwidthx
    ymax = ymin + npixy * pixwidthy

    # Use a minimum smoothing length on the grid to make sure that particles
    # contribute to at least one pixel
    hmin = 0.5 * max(pixwidthx, pixwidthy)

    xpix = xminpix + np.arange(1, npixx + 1) * pixwidthx
    nsubgrid = 0
    nok = 0
    hminall = 1e10

    # Loop over particles
    for idx in range(npart):

        # Skip particles with itype < 0
        if itype[idx] < 0:
            continue

        # Set h related quantities
        hi = hh[idx]
        horigi = hi
        if not hi > 0.0:
            continue

        # Radius of the smoothing kernel
        radkern = RADKERNEL * hi

        # Cycle as soon as we know the particle does not contribute
        xi = x[idx]
        xpixmin = xi - radkern
        if xpixmin > xmax:
            continue
        xpixmax = xi + radkern
        if xpixmax < xmin:
            continue

        yi = y[idx]
        ypixmin = yi - radkern
        if ypixmin > ymax:
            continue
        ypixmax = yi + radkern
        if ypixmax < ymin:
            continue

        # Take resolution length as max of h and 1/2 pixel width
        if hi < hmin:
            hminall = min(hi, hminall)
            nsubgrid = nsubgrid + 1
            hsmooth = hmin
        else:
            hsmooth = hi
            nok = nok + 1
        radkern = RADKERNEL * hsmooth

        # Set kernel related quantities
        hi1 = 1.0 / hsmooth
        hi21 = hi1 * hi1
        termnorm = weight[idx] * horigi
        term = termnorm * dat[idx]

        # Loop over pixels, adding the contribution from this particle
        # copy by quarters if all pixels within domain
        ipixmin = int((xi - radkern - xmin) / pixwidthx)
        ipixmax = int((xi + radkern - xmin) / pixwidthx) + 1
        jpixmin = int((yi - radkern - ymin) / pixwidthy)
        jpixmax = int((yi + radkern - ymin) / pixwidthy) + 1

        # Make sure they only contribute to pixels in the image
        # (note that this optimises much better than using min/max)
        if ipixmin < 0:
            ipixmin = 0
        if jpixmin < 0:
            jpixmin = 0
        if ipixmax > npixx:
            ipixmax = npixx
        if jpixmax > npixy:
            jpixmax = npixy

        # Precalculate an array of dx2 for this particle (optimisation)
        for ipix in range(ipixmin, ipixmax):
            dx2i[ipix] = ((xpix[ipix] - xi) ** 2) * hi21

        for jpix in range(jpixmin, jpixmax):
            ypix = yminpix + jpix * pixwidthy
            dy = ypix - yi
            dy2 = dy * dy * hi21
            for ipix in range(ipixmin, ipixmax):
                # dx2 pre-calculated; dy2 pre-multiplied by hi21
                q2 = dx2i[ipix] + dy2
                # SPH kernel - integral through cubic spline
                # interpolate from a pre-calculated table
                if q2 < RADKERNEL2:
                    wab = wfromtable(q2, coltable)
                    # Calculate data value at this pixel using the summation interpolant
                    datsmooth[ipix, jpix] = datsmooth[ipix, jpix] + term * wab
                    if normalise:
                        datnorm[ipix, jpix] = datnorm[ipix, jpix] + termnorm * wab

    # Normalise dat array
    if normalise:
        # Normalise everywhere (required if not using SPH weighting)
        for idxi in range(npixx):
            for idxj in range(npixy):
                if datnorm[idxi, idxj] > 0.0:
                    datsmooth[idxi, idxj] /= datnorm[idxi, idxj]

    # Warn about subgrid interpolation
    if nsubgrid > 1:
        nfull = int((xmax - xmin) / (hminall)) + 1
        if nsubgrid > 0.1 * nok and IVERBOSE > -1:
            print('Warning: pixel size > 2h for ', nsubgrid, ' particles')
            print('need ', nfull, ' pixels for full resolution')

    # Return datsmooth
    return datsmooth.T


@numba.njit
def interpolate_cross_section(
    x: ndarray,
    y: ndarray,
    z: ndarray,
    hh: ndarray,
    weight: ndarray,
    dat: ndarray,
    itype: ndarray,
    npart: int,
    xmin: float,
    ymin: float,
    zslice: float,
    npixx: int,
    npixy: int,
    pixwidthx: float,
    pixwidthy: float,
    normalise: bool,
):
    """Interpolate particles to grid via cross section.

    Parameters
    ----------
    x
        The particle x positions.
    y
        The particle y positions.
    z
        The particle z positions.
    hh
        The particle smoothing length.
    weight
        The particle weight.
    dat
        The scalar data to interpolate.
    itype
        The particle type.
    npart
        The number of particles.
    xmin
        The minimum x position.
    ymin
        The minimum y position.
    zslice
        Cross section location.
    npixx
        The number of pixels in the x direction.
    npixy
        The number of pixels in the y direction.
    normalise
        Whether to normalize.

    Return
    ------
    datsmooth
        The data smoothed to a pixel grid.
    """
    datsmooth = np.zeros((npixx, npixy))
    datnorm = np.zeros((npixx, npixy))
    dx2i = np.zeros(npixx)
    const = CNORMK3D

    # Loop over particles
    for idx in range(npart):

        # Skip particles with itype < 0
        if itype[idx] < 0:
            continue

        # Set h related quantities
        hi = hh[idx]
        if not hi > 0.0:
            continue
        hi1 = 1.0 / hi
        hi21 = hi1 * hi1
        radkern = RADKERNEL * hi

        # For each particle, work out distance from the cross section slice
        dz = zslice - z[idx]
        dz2 = dz ** 2 * hi21

        # If this is < 2h then add the particle's contribution to the pixels
        # otherwise skip all this and start on the next particle
        if dz2 < RADKERNEL2:

            xi = x[idx]
            yi = y[idx]
            termnorm = const * weight[idx]
            term = termnorm * dat[idx]

            # Loop over pixels, adding the contribution from this particle
            # copy by quarters if all pixels within domain
            ipixmin = int((xi - radkern - xmin) / pixwidthx)
            ipixmax = int((xi + radkern - xmin) / pixwidthx) + 1
            jpixmin = int((yi - radkern - ymin) / pixwidthy)
            jpixmax = int((yi + radkern - ymin) / pixwidthy) + 1

            # Make sure they only contribute to pixels in the image
            # (note that this optimises much better than using min/max)
            if ipixmin < 0:
                ipixmin = 0
            if jpixmin < 0:
                jpixmin = 0
            if ipixmax > npixx:
                ipixmax = npixx
            if jpixmax > npixy:
                jpixmax = npixy

            # Precalculate an array of dx2 for this particle (optimisation)
            for ipix in range(ipixmin, ipixmax):
                dx2i[ipix] = ((xmin + (ipix - 0.5) * pixwidthx - xi) ** 2) * hi21 + dz2

            # Loop over pixels, adding the contribution from this particle
            for jpix in range(jpixmin, jpixmax):
                ypix = ymin + (jpix - 0.5) * pixwidthy
                dy = ypix - yi
                dy2 = dy * dy * hi21
                for ipix in range(ipixmin, ipixmax):
                    q2 = dx2i[ipix] + dy2
                    # SPH kernel - cubic spline
                    if q2 < RADKERNEL2:
                        wab = w_cubic(q2)
                        # Calculate data value at this pixel using the summation
                        # interpolant
                        datsmooth[ipix, jpix] = datsmooth[ipix, jpix] + term * wab
                        if normalise:
                            datnorm[ipix, jpix] = datnorm[ipix, jpix] + termnorm * wab

    # Normalise dat array
    if normalise:
        # Normalise everywhere (required if not using SPH weighting)
        for idxi in range(npixx):
            for idxj in range(npixy):
                if datnorm[idxi, idxj] > 0.0:
                    datsmooth[idxi, idxj] /= datnorm[idxi, idxj]

    # Return datsmooth
    return datsmooth.T
