'''
image.py

Daniel Mentiplay, 2019.
'''

import matplotlib.pyplot as plt
import matplotlib.colors as colors
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np

from ..visualization.interpolate3D_projection import Projections3D

class Image:
    '''
    Rendered image.

    Arguments:
        dump : Dump object
    '''

    def __init__(self, dump):

        self.ParticleData = dump.ParticleData
        self.SinkData     = dump.SinkData
        self.Parameters   = dump.Parameters
        self.Units        = dump.Units

        self._default_plot_options()

    def _interpolation_weights(self, densityWeighted):

        if densityWeighted:
            self.ParticleData['interpolationWeights'] = \
                self.ParticleData['m'] / self.ParticleData['h']**2
        else:
            self.ParticleData['interpolationWeights'] = \
                1 / self.Parameters['hfact']

    def _default_plot_options(self):
        '''
        Set default plot options.
        '''

        PlotOptions = {
            'densityWeighted': False,
            'normalise':       False,
            'zobserver':       100.,
            'dscreen':         100.,
            'useaccelerate':   False
            }

        self._interpolation_weights(PlotOptions['densityWeighted'])

        self.PlotOptions = PlotOptions

    def set_plot_options(self, densityWeighted=None, normalise=None,
                         zobserver=None, dscreen=None, useaccelerate=None):
        '''
        Set plot options.
        '''

        if densityWeighted is not None:
            self.PlotOptions['densityWeighted'] = densityWeighted
            self._interpolation_weights(densityWeighted)

        if normalise is not None:
            self.PlotOptions['normalise'] = normalise

        if zobserver is not None:
            self.PlotOptions['zobserver'] = zobserver

        if dscreen is not None:
            self.PlotOptions['dscreen'] = dscreen

        if useaccelerate is not None:
            self.PlotOptions['useaccelerate'] = useaccelerate

    def get_plot_options(self, option):
        '''
        Get plot options.
        '''

        if option in ['densityWeighted', 'normalise', 'zobserver', 'dscreen',
                      'useaccelerate']:
            return self.PlotOptions[option]

        print(f'No option: {option}')
        return None

    def plot(self,
             horizontalAxis,
             verticalAxis,
             render,
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

        horizontalData = self.ParticleData['x']
        verticalData   = self.ParticleData['y']
        depthData      = self.ParticleData['z']
        renderData     = self.ParticleData[render]

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

        imageData = _interpolate_to_pixelgrid(
            horizontalData, verticalData, depthData, renderData,
            self.ParticleData['interpolationWeights'], self.ParticleData['h'],
            horizontalRange, verticalRange, self.PlotOptions['normalise'],
            self.PlotOptions['zobserver'], self.PlotOptions['dscreen'],
            self.PlotOptions['useaccelerate'] )

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

        ax.set_xlabel('x')
        ax.set_ylabel('y')

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
