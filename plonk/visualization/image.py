'''
image.py

Daniel Mentiplay, 2019.
'''

import matplotlib.pyplot as plt
import matplotlib.colors as colors
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np

from ..particles import density_from_smoothing_length
from ..visualization.splash import splash

interpolate2d = splash.interpolate2d
set_interpolation_weights = splash.set_interpolation_weights

# ---------------------------------------------------------------------------- #

cmap            = 'gist_heat'
densityWeighted = False
npixx           = 512
npixy           = 512

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

    def plot(self, horizontalAxis, verticalAxis, render, ax=None, scale='lin',
             vmin=None, vmax=None, fpeak=None, limit=-1, limits=None,
             title=None, colorbar=False):
        '''
        Make image.
        '''

        # TODO: choose fluid type: gas, dust1, dust2, ...
        # TODO: choose xData, yData, renderData

        if horizontalAxis == 'x':
            ix = 0
        elif horizontalAxis == 'y':
            ix = 1
        elif horizontalAxis == 'z':
            ix = 2
        else:
            raise ValueError('horizontalAxis should be "x", "y", or "z"')

        if verticalAxis == 'x':
            iy = 0
        elif verticalAxis == 'y':
            iy = 1
        elif verticalAxis == 'z':
            iy = 2
        else:
            raise ValueError('verticalAxis should be "x", "y", or "z"')

        xData = self.gas.position[:, ix]
        yData = self.gas.position[:, iy]

        if render in ['rho', 'dens', 'density']:
            renderData = density_from_smoothing_length(self.gas.smoothingLength,
                                                       self.gas.mass,
                                                       self.hfact)
        else:
            raise ValueError('render unknown')

        # TODO: set these according to fluid choice
        interpolationWeights = self.interpolationWeights['gas']
        smoothingLength = self.gas.smoothingLength

        # TODO: set x and y range (check units)
        xRange = [-500, 500]
        yRange = [-500, 500]

        # TODO: choose these options
        normalise = False
        exact     = False
        periodicx = False
        periodicy = False

        imageData = _interpolate_to_pixelgrid(
            xData, yData, renderData, interpolationWeights, smoothingLength,
            xRange, yRange, normalise, exact, periodicx, periodicy)

        # TODO: add options
        # TODO: choose what to plot
        # TODO: add docs
        # TODO: check if need to interpolate again

        extent = xRange + yRange

        if ax is None:
            ax = plt.gca()

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

        img = ax.imshow(imageData, norm=norm, extent=extent, origin='lower',
                        cmap=cmap)

        # TODO: set labels correctly
        xlabel = r'x [au]'
        ylabel = r'y [au]'

        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)

        if limit > 0:
            limits = [-limit,limit,-limit,limit]

        if limits is not None:
            ax.set_xlim(limits[0],limits[1])
            ax.set_ylim(limits[2],limits[3])

        if title is not None:
            ax.set_title(title)

        cb = None
        if colorbar:
            divider = make_axes_locatable(ax)
            cax = divider.append_axes("right", size="5%", pad=0.05)
            cb = plt.colorbar(img, cax=cax)

        return img, cb

    def _convert_units(self):
        '''
        Convert units.
        '''

        # TODO: write this

    def _calculate_quantity(self):
        '''
        Calculate an extra quantity.
        '''

        # TODO: write this

def _interpolate_to_pixelgrid(xData, yData, renderData,
                              interpolationWeights, smoothingLength, xRange,
                              yRange, normalise=False, exact=False,
                              periodicx=False, periodicy=False):

    itype = np.ones_like(xData)
    npart = len(smoothingLength)

    xmax      = xRange[1]
    ymax      = yRange[1]
    xmin      = xRange[0]
    ymin      = yRange[0]
    pixwidthx = (xmax - xmin) / npixx
    pixwidthy = (ymax - ymin) / npixy

    imageData = np.zeros((npixx, npixy), dtype=np.float32, order='F')

    interpolate2d(x=xData,
                  y=yData,
                  hh=smoothingLength,
                  weight=interpolationWeights,
                  dat=renderData,
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

