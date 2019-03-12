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

def scalar_interpolation(horizontal_data, vertical_data, depth_data,
                         smoothing_length, weights, render_data, particle_mass,
                         horizontal_range, vertical_range, cross_section,
                         zslice, opacity, normalize, zobserver, dscreen,
                         accelerate):
    '''
    Interpolate a scalar quantity to pixels by projection.

    TODO: add docs
    '''

    # TODO: set number of pixels by options
    npixx = 512
    npixy = 512

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

    itype = np.ones_like(horizontal_data, dtype=np.int32)
    npart = len(smoothing_length)

    xmax      = horizontal_range[1]
    ymax      = vertical_range[1]
    xmin      = horizontal_range[0]
    ymin      = vertical_range[0]
    pixwidthx = (xmax - xmin) / npixx
    pixwidthy = (ymax - ymin) / npixy

    # interpolate3d_projection expects single precision reals
    x   = np.array(horizontal_data, dtype=np.single)
    y   = np.array(vertical_data,   dtype=np.single)
    z   = np.array(depth_data,      dtype=np.single)
    dat = np.array(render_data,     dtype=np.single)

    hh            = smoothing_length
    weight        = weights
    normalise     = normalize
    useaccelerate = accelerate

################################################################################
# TODO: temporary; testing phase
    pmass = np.array(particle_mass, dtype=np.single)
    zorig = z
    pixwidth = pixwidthx
    dscreenfromobserver = dscreen
    rkappa = np.pi * hh.mean()**2 / pmass[0]
    zcut = zobserver
################################################################################

    if cross_section:
        image_data = interpolate3d_fastxsec(x=x,
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
            image_data = interp3d_proj_opacity(x=x,
                                               y=y,
                                               z=z,
                                               hh=hh,
                                               pmass=pmass,
                                               npmass=npart,
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
            image_data = interpolate3d_projection(x=x,
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
    image_data = image_data.T

    image_data = np.array(image_data)

    return image_data

def vector_interpolation(horizontal_data, vertical_data, depth_data,
                         smoothing_length, weights, vecx, vecy,
                         horizontal_range, vertical_range, cross_section,
                         zslice, normalize, zobserver, dscreen):
    '''
    Interpolate a vector quantity to pixels by projection.

    TODO: add docs
    '''

    # TODO: set number of pixels by options
    npixx = 512
    npixy = 512

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

    itype = np.ones_like(horizontal_data, dtype=np.int32)
    npart = len(smoothing_length)

    xmax      = horizontal_range[1]
    ymax      = vertical_range[1]
    xmin      = horizontal_range[0]
    ymin      = vertical_range[0]
    pixwidthx = (xmax - xmin) / npixx
    pixwidthy = (ymax - ymin) / npixy

    # interpolate3d_projection expects single precision reals
    x    = np.array(horizontal_data, dtype=np.single)
    y    = np.array(vertical_data,   dtype=np.single)
    z    = np.array(depth_data,      dtype=np.single)
    vecx = np.array(vecx,            dtype=np.single)
    vecy = np.array(vecy,            dtype=np.single)

    hh            = smoothing_length
    weight        = weights
    normalise     = normalize

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

    vecsmoothx = np.array(vecsmoothx)
    vecsmoothy = np.array(vecsmoothy)

    return vecsmoothx, vecsmoothy
