'''
image.py

Daniel Mentiplay, 2019.
'''

import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np

from ..particles import density_from_smoothing_length
from ..visualization.splash import splash

interpolate2d = splash.interpolate2d
set_interpolation_weights = splash.set_interpolation_weights

# ---------------------------------------------------------------------------- #

cmap            = 'inferno'
densityWeighted = False
npixx           = 1000
npixy           = 1000

# ---------------------------------------------------------------------------- #

class Image:
    '''
    Rendered image.

    Arguments:
        dump : Dump object
    '''

    def __init__(self, dump):

        self._get_data_from_dump(dump)

    def _get_data_from_dump(self, dump):

        self.hfact = dump.parameters.numerical['hfact']
        self.fullDump = dump.fullDump

        self.gas = dump.gas
        self.dust = dump.dust

        nDustSmall = dump.parameters.dust['nDustSmall']
        nDustLarge = dump.parameters.dust['nDustLarge']
        nDustTypes = nDustSmall + nDustLarge
        containsDust = bool(nDustTypes > 0)

        interpolationWeights = dict()

        interpolationWeights['gas'] = self._interpolation_weights(dump.gas)

        if containsDust:
            interpolationWeights['dust'] = list()
            for idx in range(nDustLarge):
                interpolationWeights['dust'].append(
                    self._interpolation_weights(dump.dust[idx]) )

        self.interpolationWeights = interpolationWeights

    def _interpolation_weights(self, fluid):

        nParticles      = fluid.number
        massParticle    = fluid.mass
        smoothingLength = fluid.smoothingLength

        interpolationWeights = np.zeros(nParticles)

        if densityWeighted:
            interpolationWeights[:] = massParticle / smoothingLength**3
        else:
            interpolationWeights[:] = 1 / self.hfact**3

        return interpolationWeights

    def _interpolate_to_pixelgrid(self):

        # TODO: choose fluid type
        # TODO: choose xData, yData, render
        xData = self.gas.position[:, 0]
        yData = self.gas.position[:, 1]

        # TODO: choose these options
        normalise = False
        exact     = False
        periodicx = False
        periodicy = False

        smoothingLength = self.gas.smoothingLength
        render = density_from_smoothing_length(smoothingLength, self.gas.mass,
                                               self.hfact)

        itype = np.ones_like(xData)
        interpolationWeights = self.interpolationWeights['gas']
        npart = len(smoothingLength)

        maxx      = np.max(xData)
        maxy      = np.max(yData)
        xmin      = - maxx
        ymin      = - maxy
        pixwidthx = 2*maxx / npixx
        pixwidthy = 2*maxy / npixy

        imageData = np.zeros((npixx, npixy), dtype=np.float32, order='F')

        interpolate2d(x=xData,
                      y=yData,
                      hh=smoothingLength,
                      weight=interpolationWeights,
                      dat=render,
                      itype=itype,
                      npart=npart,
                      xmin=xmin,
                      ymin=ymin,
                      datsmooth=imageData,
                      npixx=npixx,
                      npixy=npixy,
                      pixwidthx=pixwidthx,
                      pixwidthy=pixwidthy,
                      normalise=normalise,
                      exact=exact,
                      periodicx=periodicx,
                      periodicy=periodicy)

        return imageData

    def plot(self):
        '''
        Make image.
        '''

        # TODO: add options
        # TODO: choose what to plot
        # TODO: add docs
        # TODO: check if need to interpolate again

        imageData = self._interpolate_to_pixelgrid()

        vmin = None
        vmax = None
        fpeak = None
        scale = 'log'
        _scale = 'lin'

        if vmax is None:
            vmax = imageData.max()

        if fpeak is not None:
            vmax = imageData.max() * fpeak

        if vmin is None:
            vmin = imageData.min()

        if scale is None:
            scale = _scale

        if scale == 'log':
            norm = colors.SymLogNorm(1e-1*vmax, clip=True)
        elif scale == 'lin':
            norm = colors.Normalize(vmin=vmin, vmax=vmax, clip=True)
        else:
            raise ValueError("Unknown color scale: " + scale)

        plt.imshow(imageData, norm=norm, origin='lower', cmap=cmap)
