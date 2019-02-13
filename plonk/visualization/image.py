'''
image.py

Daniel Mentiplay, 2019.
'''

import matplotlib.pyplot as plt
import matplotlib.colors as colors
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np

from ..ParticleData import iGas, iDust
from .splash import splash

options = ['accelerate',
           'colorbar',
           'colorscale',
           'densityWeighted',
           'dscreen',
           'normalize',
           'zobserver']

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

        self._axis     = None
        self._image    = None
        self._colorbar = None

    def _default_plot_options(self):
        '''
        Set default plot options.
        '''

        PlotOptions = {
            'accelerate':      False,
            'colorbar':        True,
            'colormap':        'gist_heat',
            'colorscale':      'linear',
            'densityWeighted': False,
            'dscreen':         None,
            'normalize':       False,
            'zobserver':       None
            }

        self.PlotOptions = PlotOptions

    def print_plot_options(self):
        '''
        Print plot options.
        '''

        print('Current plot options:\n')
        for opt in self.PlotOptions:
            print(f'{opt:20}:  {self.PlotOptions[opt]}')

    def set_plot_options(self, accelerate=None, colorbar=None, colorscale=None,
                         densityWeighted=None, dscreen=None, normalize=None,
                         zobserver=None):
        '''
        Set plot options.
        '''

        args = list(locals().values())[1:]
        if all([arg is None for arg in args]):
            self.print_plot_options()

        if accelerate is not None:
            self.PlotOptions['accelerate'] = accelerate

        if colorbar is not None:
            self.PlotOptions['colorbar'] = colorbar

        if colorscale is not None:
            self.PlotOptions['colorscale'] = colorscale

        if densityWeighted is not None:
            self.PlotOptions['densityWeighted'] = densityWeighted

        if dscreen is not None:
            self.PlotOptions['dscreen'] = dscreen

        if normalize is not None:
            self.PlotOptions['normalize'] = normalize

        if zobserver is not None:
            self.PlotOptions['zobserver'] = zobserver

    def get_plot_option(self, option):
        '''
        Get plot options.
        '''

        if option in options:
            return self.PlotOptions[option]

        print(f'No option: {option}')
        return None

    def plot(self,
             horizontalAxis=None,
             verticalAxis=None,
             render=None,
             particleTypes=None,
             itype=None,
             horizontalRange=None,
             verticalRange=None,
             imageRange=-1,
             renderScale=None,
             renderMin=None,
             renderMax=None,
             renderFractionMax=None,
             title=None,
             ax=None,
             colorbar=None,
             colormap=None):
        '''
        Make image.
        '''

        # TODO: add options
        # TODO: choose what to plot
        # TODO: add docs
        # TODO: check if need to interpolate again
        # TODO: choose fluid type: gas, dust1, dust2, ...

        if horizontalAxis is None:
            horizontalAxis = 'x'

        if verticalAxis is None:
            verticalAxis = 'y'

        if render is None:
            render = 'rho'

        itypes = list()

        if itype is not None and particleTypes is not None:
            raise ValueError('Cannot set itype and particleTypes together')

        if particleTypes is None and itype is None:
            itypes = [iGas]

        if itype is not None:
            itypes = [itype]

        if particleTypes is not None:

            if 'gas' in particleTypes:
                itypes.append(iGas)

            # TODO: can only plot one type at the moment
            if 'dust' in particleTypes:
                for i in range(self.Parameters['ndustlarge']):
                    itypes.append(iDust + i)

        if len(itypes) > 1:
            raise ValueError('plotting multiple types at once is not working')

        print(f'Plotting {render} on [{horizontalAxis}, {verticalAxis}] window')

        xyz = set(['x', 'y', 'z'])
        depthAxis = xyz.difference([horizontalAxis, verticalAxis]).pop()

        pd = self.ParticleData.loc[self.ParticleData['itype'].isin(itypes)]

        horizontalData  = pd[horizontalAxis]
        verticalData    = pd[verticalAxis]
        depthData       = pd[depthAxis]
        renderData      = pd[render]
        smoothingLength = pd['h']

        weights = _interpolation_weights(
            self.PlotOptions['densityWeighted'], pd,
            self.Parameters['hfact']
            )

        if imageRange > 0:
            if horizontalRange is not None or verticalRange is not None:
                raise ValueError( 'Cannot set imageRange and horizontalRange ' \
                                + '(or verticalRange) at the same time' )
            horizontalRange = [-imageRange, imageRange]
            verticalRange   = [-imageRange, imageRange]

        if horizontalRange is None:
            horizontalRange = [horizontalData.min(), horizontalData.max()]

        if verticalRange is None:
            verticalRange = [verticalData.min(), verticalData.max()]

        imageData = self._interpolate_to_pixelgrid(
            horizontalData, verticalData, depthData, smoothingLength,
            weights, renderData, horizontalRange, verticalRange )

        extent = horizontalRange + verticalRange

        cmap = self.PlotOptions['colormap']
        if colormap is not None:
            cmap = colormap

        if renderMax is None:
            vmax = imageData.max()

        if renderFractionMax is not None:
            vmax = imageData.max() * renderFractionMax

        if renderMin is None:
            vmin = imageData.min()

        if renderScale is None:
            renderScale = self.PlotOptions['colorscale']

        if renderScale == 'log':
            norm = colors.SymLogNorm(1e-1*vmax, clip=True)
        elif renderScale == 'linear':
            norm = colors.Normalize(vmin=vmin, vmax=vmax, clip=True)
        else:
            raise ValueError("Unknown color renderScale: " \
                + renderScale)

        if ax is None:
            plt.clf()
            ax = plt.gca()

        img = ax.imshow(imageData, norm=norm, origin='lower', extent=extent,
                        cmap=cmap)

        ax.set_xlabel(horizontalAxis)
        ax.set_ylabel(verticalAxis)

        ax.set_xlim(horizontalRange[0], horizontalRange[1])
        ax.set_ylim(verticalRange[0], verticalRange[1])

        if title is not None:
            ax.set_title(title)

        cb = None
        if colorbar is not None:
            colorbar_ = colorbar
        else:
            colorbar_ = self.PlotOptions['colorbar']
        if colorbar_:
            divider = make_axes_locatable(ax)
            cax = divider.append_axes("right", size="5%", pad=0.05)
            cb = plt.colorbar(img, cax=cax)

        self._axis = ax
        self._image = img
        self._colorbar = cb

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

    def set_colorbar(self, vmin=None, vmax=None):
        '''
        Set colorbar limits.
        '''

        if vmin is None:
            vmin = self._colorbar.vmin

        if vmax is None:
            vmax = self._colorbar.vmax

        self._image.set_clim([vmin, vmax])

    def _interpolate_to_pixelgrid(self, horizontalData, verticalData, depthData,
                                  smoothingLength, weights, renderData,
                                  horizontalRange, verticalRange):

        # TODO: set number of pixels based on smoothing length
        npixx = 512
        npixy = 512

        normalize = self.PlotOptions['normalize']
        zobserver = self.PlotOptions['zobserver']
        dscreen = self.PlotOptions['dscreen']
        accelerate = self.PlotOptions['accelerate']

        if zobserver is None:
            zobserver = 1e10

        if dscreen is None:
            dscreen = 1e10

        if normalize is None:
            normalize = False

        if accelerate is None:
            accelerate = False

        itype = np.ones_like(horizontalData)
        npart = len(smoothingLength)

        xmax      = horizontalRange[1]
        ymax      = verticalRange[1]
        xmin      = horizontalRange[0]
        ymin      = verticalRange[0]
        pixwidthx = (xmax - xmin) / npixx
        pixwidthy = (ymax - ymin) / npixy

        imageData = splash.interpolate3d_projection(x=horizontalData,
                                                    y=verticalData,
                                                    z=depthData,
                                                    hh=smoothingLength,
                                                    weight=weights,
                                                    dat=renderData,
                                                    itype=itype,
                                                    npart=npart,
                                                    xmin=xmin,
                                                    ymin=ymin,
                                                    npixx=npixx,
                                                    npixy=npixy,
                                                    pixwidthx=pixwidthx,
                                                    pixwidthy=pixwidthy,
                                                    normalise=normalize,
                                                    zobserver=zobserver,
                                                    dscreen=dscreen,
                                                    useaccelerate=accelerate)

        imageData = imageData.T

        return imageData

def _interpolation_weights(densityWeighted, ParticleData, hfact):

    if densityWeighted:
        return np.array(ParticleData['m'] / ParticleData['h']**2)

    return np.full_like(ParticleData['h'], 1/hfact)
