'''
splash.py

Daniel Mentiplay, 2019.
'''

import numpy as np

try:
    from ._splash import interpolate3d_projection, interpolate3d_proj_vec, \
                         interpolate3d_fastxsec, interpolate3d_xsec_vec,   \
                         interp3d_proj_opacity
except ImportError:
    raise Exception('Cannot import Splash interpolation routines. See ' + \
                    'documentation.')

def scalar_interpolation(positions, smoothing_length, weights, scalar_data,
                         particle_mass, horizontal_range, vertical_range, npix,
                         cross_section, zslice, opacity, normalize, zobserver,
                         dscreen, accelerate):
    '''
    Interpolate a scalar quantity to pixels by projection.

    TODO: add docs
    '''

    if npix is None:
        npixx = 512
        npixy = 512
    else:
        npixx = npix[0]
        npixy = npix[1]

    if cross_section is None:
        cross_section = False

    if cross_section and zslice is None:
        zslice = 0.

    if opacity is None:
        opacity = False

    if zobserver is None:
        zobserver = 1e10

    if dscreen is None:
        dscreen = 1e10

    if normalize is None:
        normalize = False

    if accelerate is None:
        accelerate = False

    npart = len(smoothing_length)
    itype = np.ones(npart, dtype=np.int32)

    xmin      = horizontal_range[0]
    ymin      = vertical_range[0]
    xmax      = horizontal_range[1]
    ymax      = vertical_range[1]
    pixwidthx = (xmax - xmin) / npixx
    pixwidthy = (ymax - ymin) / npixy

    # Splash routines expect single precision
    x      = np.array(positions[0, :],  dtype=np.single)
    y      = np.array(positions[1, :],  dtype=np.single)
    z      = np.array(positions[2, :],  dtype=np.single)
    pmass  = np.array(particle_mass,    dtype=np.single)
    hh     = np.array(smoothing_length, dtype=np.single)
    dat    = np.array(scalar_data,      dtype=np.single)
    weight = np.array(weights,          dtype=np.single)

    normalise     = normalize
    useaccelerate = accelerate

################################################################################
# TODO: temporary; testing phase
    npmass = npart
    zorig = z
    pixwidth = pixwidthx
    dscreenfromobserver = dscreen
    rkappa = np.pi * hh.mean()**2 / pmass[0]
    zcut = zobserver
################################################################################

    if cross_section:
        datsmooth = \
            interpolate3d_fastxsec(x=x,
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
                                   normalise=normalise)
    else:
        if opacity:
            datsmooth = \
                interp3d_proj_opacity(x=x,
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
                                      zcut=zcut)
        else:
            datsmooth = \
                interpolate3d_projection(x=x,
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
                                         useaccelerate=useaccelerate)

    # TODO: check whether we need to transpose: Fortran vs C array ordering
    datsmooth = datsmooth.T

    smoothed_scalar = np.array(datsmooth)

    return smoothed_scalar

def vector_interpolation(positions, smoothing_length, weights, vector_data,
                         horizontal_range, vertical_range, npix, cross_section,
                         zslice, normalize, zobserver, dscreen):
    '''
    Interpolate a vector quantity to pixels by projection.

    TODO: add docs
    '''

    if npix is None:
        npixx = 512
        npixy = 512
    else:
        npixx = npix[0]
        npixy = npix[1]

    if cross_section is None:
        cross_section = False

    if cross_section and zslice is None:
        zslice = 0.

    if zobserver is None:
        zobserver = 1e10

    if dscreen is None:
        dscreen = 1e10

    if normalize is None:
        normalize = False

    npart = len(smoothing_length)
    itype = np.ones(npart, dtype=np.int32)

    xmin      = horizontal_range[0]
    ymin      = vertical_range[0]
    xmax      = horizontal_range[1]
    ymax      = vertical_range[1]
    pixwidthx = (xmax - xmin) / npixx
    pixwidthy = (ymax - ymin) / npixy

    # Splash routines expect single precision
    x      = np.array(positions[0, :],   dtype=np.single)
    y      = np.array(positions[1, :],   dtype=np.single)
    z      = np.array(positions[2, :],   dtype=np.single)
    hh     = np.array(smoothing_length,  dtype=np.single)
    vecx   = np.array(vector_data[0, :], dtype=np.single)
    vecy   = np.array(vector_data[1, :], dtype=np.single)
    weight = np.array(weights,           dtype=np.single)

    normalise = normalize

    if cross_section:
        vecsmoothx, vecsmoothy = interpolate3d_xsec_vec(x=x,
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
                                                        normalise=normalise)
    else:
        vecsmoothx, vecsmoothy = interpolate3d_proj_vec(x=x,
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
                                                        dscreen=dscreen)

    smoothed_vector = np.stack((np.array(vecsmoothx), np.array(vecsmoothy)))

    return smoothed_vector
