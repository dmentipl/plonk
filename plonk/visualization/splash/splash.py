'''
splash.py

Daniel Mentiplay, 2019.
'''

import numpy as np

try:
    from ._splash import interpolate3d_projection # pylint: disable-msg=no-name-in-module
except ImportError:
    raise Exception('Set LD_LIBRARY_PATH or DYLD_LIBRARY_PATH to contain ' + \
                    'plonk/plonk/visualization/splash')

def interpolate_to_pixelgrid(horizontal_data, vertical_data, depth_data,
                             smoothing_length, weights, render_data,
                             horizontal_range, vertical_range, normalize,
                             zobserver, dscreen, accelerate):
    '''
    Interpolate a quantity to a pixel grid by projection.

    TODO: add docs
    '''

    # TODO: set number of pixels based on smoothing length
    npixx = 512
    npixy = 512

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

    # TODO: check whether we need to transpose
    image_data = image_data.T

    image_data = np.array(image_data)

    return image_data
