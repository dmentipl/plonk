'''
image.py

Daniel Mentiplay, 2019.
'''

import matplotlib.pyplot as plt
import matplotlib.colors as colors
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np

from ..ParticleData import density_from_smoothing_length
from ..visualization.interpolate3D_projection import Projections3D

class Image:
    '''
    Rendered image.

    Arguments:
        dump : Dump object
    '''

    def __init__(self, dump):

        self._get_data_from_dump(dump)

    def _get_data_from_dump(self, dump):

        self.hfact = dump.Parameters['hfact']

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

        # TODO: density weighting as option
        densityWeighted = False

        nParticles      = fluid.number
        massParticle    = fluid.mass
        smoothingLength = fluid.smoothingLength

        interpolationWeights = np.zeros(len(smoothingLength))

        if densityWeighted:
            interpolationWeights[:] = massParticle / smoothingLength**2
        else:
            interpolationWeights[:] = 1 / self.hfact**2

        return interpolationWeights

    def plot(self,
             horizontalAxisLabel,
             verticalAxisLabel,
             renderLabel,
             horizontalRange=None,
             verticalRange=None,
             imageRange=-1,
             renderScaleLabel='lin',
             renderMin=None,
             renderMax=None,
             renderFractionMax=None,
             title=None,
             ax=None,
             colorbar=False,
             colormap=None):
        '''
        Make image.
        '''

        _scale           = 'lin'
        _cmap            = 'gist_heat'

        # TODO: add options
        # TODO: choose what to plot
        # TODO: add docs
        # TODO: check if need to interpolate again
        # TODO: choose fluid type: gas, dust1, dust2, ...

        # TODO: set these according to fluid choice
        interpolationWeights = self.interpolationWeights['gas']
        smoothingLength = self.gas.smoothingLength

        # TODO: choose these options
        normalise = False
        exact     = False
        periodicx = False
        periodicy = False

        if horizontalAxisLabel == 'x':
            iHorizontal = 0
            horizontalAxisLabel = r'x [au]'
        elif horizontalAxisLabel == 'y':
            iHorizontal = 1
            horizontalAxisLabel = r'y [au]'
        elif horizontalAxisLabel == 'z':
            iHorizontal = 2
            horizontalAxisLabel = r'z [au]'
        else:
            raise ValueError('horizontalAxisLabel should be "x", "y", or "z"')

        if verticalAxisLabel == 'x':
            iVertical = 0
            verticalAxisLabel = r'x [au]'
        elif verticalAxisLabel == 'y':
            iVertical = 1
            verticalAxisLabel = r'y [au]'
        elif verticalAxisLabel == 'z':
            iVertical = 2
            verticalAxisLabel = r'z [au]'
        else:
            raise ValueError('verticalAxisLabel should be "x", "y", or "z"')

        horizontalData = self.gas.position[:, iHorizontal]
        verticalData   = self.gas.position[:, iVertical]

        if renderLabel in ['rho', 'dens', 'density']:
            renderData = density_from_smoothing_length(self.gas.smoothingLength,
                                                       self.gas.mass,
                                                       self.hfact)
        else:
            raise ValueError('renderLabel unknown')

        if imageRange > 0:
            if horizontalRange is not None or verticalRange is not None:
                raise ValueError( 'Cannot set imageRange and horizontalRange ' \
                                + '(or verticalRange at the same time' )
            horizontalRange = [-imageRange, imageRange]
            verticalRange   = [-imageRange, imageRange]

        if horizontalRange is None:
            horizontalRange = [horizontalData.min(), horizontalData.max()]

        if verticalRange is None:
            verticalRange = [verticalData.min(), verticalData.max()]

        depthData = self.gas.position[:, 2]

        imageData = _interpolate_to_pixelgrid(
            horizontalData, verticalData, depthData, renderData,
            interpolationWeights, smoothingLength, horizontalRange,
            verticalRange, normalise, exact, periodicx, periodicy )

        extent = horizontalRange + verticalRange

        cmap = _cmap
        if colormap is not None:
            cmap = colormap

        if renderMax is None:
            vmax = imageData.max()

        if renderFractionMax is not None:
            vmax = imageData.max() * renderFractionMax

        if renderMin is None:
            vmin = imageData.min()

        if renderScaleLabel is None:
            renderScaleLabel = _scale

        if renderScaleLabel == 'log':
            norm = colors.SymLogNorm(1e-1*vmax, clip=True)
        elif renderScaleLabel == 'lin':
            norm = colors.Normalize(vmin=vmin, vmax=vmax, clip=True)
        else:
            raise ValueError("Unknown color renderScaleLabel: " \
                + renderScaleLabel)

        if ax is None:
            ax = plt.gca()

        img = ax.imshow(imageData, norm=norm, origin='lower', extent=extent,
                        cmap=cmap)

        ax.set_xlabel(horizontalAxisLabel)
        ax.set_ylabel(verticalAxisLabel)

        ax.set_xlim(horizontalRange[0], horizontalRange[1])
        ax.set_ylim(verticalRange[0], verticalRange[1])

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

def _interpolate_to_pixelgrid(horizontalData, verticalData, depthData,
                              renderData, interpolationWeights, smoothingLength,
                              horizontalRange, verticalRange, normalise=False,
                              zobserver=100., dscreen=100., useaccelerate=False):

    # TODO: set number of pixels based on smoothing length
    npixx = 512
    npixy = 512

    itype = np.ones_like(horizontalData)
    npart = len(smoothingLength)

    xmax      = horizontalRange[1]
    ymax      = verticalRange[1]
    xmin      = horizontalRange[0]
    ymin      = verticalRange[0]
    pixwidthx = (xmax - xmin) / npixx
    pixwidthy = (ymax - ymin) / npixy

    imageData = np.zeros((npixx, npixy), dtype=np.float32, order='F')

    Projections3D.interpolate3d_projection(x=horizontalData,
                                           y=verticalData,
                                           z=depthData,
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
                                           zobserver=zobserver,
                                           dscreen=dscreen,
                                           useaccelerate=useaccelerate)

    imageData = imageData.T

    return imageData
