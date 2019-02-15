'''
image.py

Daniel Mentiplay, 2019.
'''

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np

from ..particles import I_GAS, I_DUST
from .splash import splash

options = ['accelerate',
           'colorbar',
           'colorscale',
           'density_weighted',
           'dscreen',
           'fontfamily',
           'fontsize',
           'normalize',
           'zobserver']

class Image:
    '''
    Rendered image.

    Arguments:
        dump : Dump object
    '''

    def __init__(self, dump):

        self.particles  = dump.particles
        self.sinks      = dump.sinks
        self.parameters = dump.parameters
        self.units      = dump.units

        self._default_plot_options()

        self._axis     = None
        self._image    = None
        self._colorbar = None

    def _default_plot_options(self):
        '''
        Set default plot options.
        '''

        plot_options = {
            'accelerate':      False,
            'colorbar':        True,
            'colormap':        'gist_heat',
            'colorscale':      'linear',
            'density_weighted': False,
            'dscreen':         None,
            'fontfamily':      'sans-serif',
            'fontsize':        12,
            'normalize':       False,
            'zobserver':       None
            }

        self.plot_options = plot_options

    def print_plot_options(self):
        '''
        Print plot options.
        '''

        print('Current plot options:\n')
        for opt in self.plot_options:
            print(f'{opt:20}:  {self.plot_options[opt]}')

    def set_plot_options(self,
                         accelerate=None,
                         colorbar=None,
                         colorscale=None,
                         density_weighted=None,
                         dscreen=None,
                         fontfamily=None,
                         fontsize=None,
                         normalize=None,
                         zobserver=None):
        '''
        Set plot options.
        '''

        args = list(locals().values())[1:]
        if all([arg is None for arg in args]):
            self.print_plot_options()

        if accelerate is not None:
            self.plot_options['accelerate'] = accelerate

        if colorbar is not None:
            self.plot_options['colorbar'] = colorbar

        if colorscale is not None:
            self.plot_options['colorscale'] = colorscale

        if density_weighted is not None:
            self.plot_options['density_weighted'] = density_weighted

        if dscreen is not None:
            self.plot_options['dscreen'] = dscreen

        if fontfamily is not None:
            self.plot_options['fontfamily'] = fontfamily

        if fontsize is not None:
            self.plot_options['fontsize'] = fontsize

        if normalize is not None:
            self.plot_options['normalize'] = normalize

        if zobserver is not None:
            self.plot_options['zobserver'] = zobserver

    def get_plot_option(self, option):
        '''
        Get plot options.
        '''

        if option in options:
            return self.plot_options[option]

        print(f'No option: {option}')
        return None

    def plot(self,
             horizontal_axis=None,
             vertical_axis=None,
             render=None,
             particle_types=None,
             itype=None,
             horizontal_range=None,
             vertical_range=None,
             image_range=-1,
             render_scale=None,
             render_min=None,
             render_max=None,
             render_fraction_max=None,
             title=None,
             ax=None,
             newfig=None,
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

        if horizontal_axis is None:
            horizontal_axis = 'x'

        if vertical_axis is None:
            vertical_axis = 'y'

        if render is None:
            render = 'rho'

        itypes = list()

        if itype is not None and particle_types is not None:
            raise ValueError('Cannot set itype and particle_types together')

        if particle_types is None and itype is None:
            itypes = [I_GAS]

        if itype is not None:
            itypes = [itype]

        if particle_types is not None:

            if 'gas' in particle_types:
                itypes.append(I_GAS)

            # TODO: can only plot one type at the moment
            if 'dust' in particle_types:
                for i in range(self.parameters['ndustlarge']):
                    itypes.append(I_DUST + i)

        if len(itypes) > 1:
            raise ValueError('plotting multiple types at once is not working')

        print(f'Rendering {render} on [{horizontal_axis}, {vertical_axis}] window')

        xyz = set(['x', 'y', 'z'])
        depth_axis = xyz.difference([horizontal_axis, vertical_axis]).pop()

        pd = self.particles.loc[self.particles['itype'].isin(itypes)]

        horizontal_data  = pd[horizontal_axis]
        vertical_data    = pd[vertical_axis]
        depth_data       = pd[depth_axis]
        render_data      = pd[render]
        smoothing_length = pd['h']

        weights = _interpolation_weights(
            self.plot_options['density_weighted'], pd,
            self.parameters['hfact']
            )

        if image_range > 0:
            if horizontal_range is not None or vertical_range is not None:
                raise ValueError( 'Cannot set image_range and horizontal_range ' \
                                + '(or vertical_range) at the same time' )
            horizontal_range = [-image_range, image_range]
            vertical_range   = [-image_range, image_range]

        if horizontal_range is None:
            horizontal_range = [horizontal_data.min(), horizontal_data.max()]

        if vertical_range is None:
            vertical_range = [vertical_data.min(), vertical_data.max()]

        image_data = self._interpolate_to_pixelgrid(
            horizontal_data, vertical_data, depth_data, smoothing_length,
            weights, render_data, horizontal_range, vertical_range )

        extent = horizontal_range + vertical_range

        cmap = self.plot_options['colormap']
        if colormap is not None:
            cmap = colormap

        if render_max is None:
            vmax = image_data.max()

        if render_fraction_max is not None:
            vmax = image_data.max() * render_fraction_max

        if render_min is None:
            vmin = image_data.min()

        if render_scale is None:
            render_scale = self.plot_options['colorscale']

        if render_scale == 'log':
            norm = colors.Sym_log_norm(1e-1*vmax, clip=True)
        elif render_scale == 'linear':
            norm = colors.Normalize(vmin=vmin, vmax=vmax, clip=True)
        else:
            raise ValueError("Unknown color render_scale: " \
                + render_scale)

        mpl.rcParams['font.family'] = self.plot_options['fontfamily']
        mpl.rcParams['font.size']   = self.plot_options['fontsize']

        if ax is None:
            if newfig:
                plt.figure()
            plt.clf()
            ax = plt.gca()

        img = ax.imshow(image_data, norm=norm, origin='lower', extent=extent,
                        cmap=cmap)

        ax.set_xlabel(horizontal_axis)
        ax.set_ylabel(vertical_axis)

        ax.set_xlim(horizontal_range[0], horizontal_range[1])
        ax.set_ylim(vertical_range[0], vertical_range[1])

        if title is not None:
            ax.set_title(title)

        render_label = r'$\int$ '+ f'{render}' + ' dz'

        cb = None
        if colorbar is not None:
            colorbar_ = colorbar
        else:
            colorbar_ = self.plot_options['colorbar']
        if colorbar_:
            divider = make_axes_locatable(ax)
            cax = divider.append_axes("right", size="5%", pad=0.05)
            cb = plt.colorbar(img, cax=cax)
            if render_label:
                cb.set_label(render_label)

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

    def _interpolate_to_pixelgrid(self, horizontal_data, vertical_data, depth_data,
                                  smoothing_length, weights, render_data,
                                  horizontal_range, vertical_range):

        # TODO: set number of pixels based on smoothing length
        npixx = 512
        npixy = 512

        normalize  = self.plot_options['normalize']
        zobserver  = self.plot_options['zobserver']
        dscreen    = self.plot_options['dscreen']
        accelerate = self.plot_options['accelerate']

        if zobserver is None:
            zobserver = 1e10

        if dscreen is None:
            dscreen = 1e10

        if normalize is None:
            normalize = False

        if accelerate is None:
            accelerate = False

        itype = np.ones_like(horizontal_data)
        npart = len(smoothing_length)

        xmax      = horizontal_range[1]
        ymax      = vertical_range[1]
        xmin      = horizontal_range[0]
        ymin      = vertical_range[0]
        pixwidthx = (xmax - xmin) / npixx
        pixwidthy = (ymax - ymin) / npixy

        image_data = splash.interpolate3d_projection(x=horizontal_data,
                                                    y=vertical_data,
                                                    z=depth_data,
                                                    hh=smoothing_length,
                                                    weight=weights,
                                                    dat=render_data,
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

        image_data = image_data.T

        return image_data

def _interpolation_weights(density_weighted, particles, hfact):

    if density_weighted:
        return np.array(particles['m'] / particles['h']**2)

    return np.full_like(particles['h'], 1/hfact)
