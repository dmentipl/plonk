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
from .splash.splash import scalar_interpolation, vector_interpolation
from ..utils import normalize_vector, rotate_vector_arbitrary_axis

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
        self._quiver   = None

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
             render=None,
             vector=None,
             particle_types=None,
             itype=None,
             horizontal_range=None,
             vertical_range=None,
             rotation_axis=None,
             rotation_angle=None,
             position_angle=None,
             inclination=None,
             number_pixels=None,
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
        # TODO: add docs
        # TODO: choose fluid type: gas, dust1, dust2, ...
        # TODO: add more checks on input

        #--- Render and vector field

        if render is None and vector is None:
            render = 'rho'

        if render is not None and render not in self.particles:
            raise ValueError(f'{render} not available for rendering')

        if vector is not None and vector not in ['v', 'velocity']:
            raise ValueError(f'{vector} not available for vector field overlay')

        #--- Rotate frame

        rotate = False

        if (rotation_axis is None or rotation_angle is None) and \
           (position_angle is None and inclination is None):
            raise ValueError('Cannot set rotation_axis/rotation_angle and ' + \
                             ' position_angle/inclination at the same time')

        if rotation_axis is not None:
            if rotation_angle is not None:
                rotate = True
                rotation_axis = normalize_vector(rotation_axis)
            else:
                raise ValueError('Must specify rotation_angle')
            if isinstance(rotation_axis, list):
                rotation_axis = np.array(rotation_axis)
            rotation_axis = normalize_vector(rotation_axis)

        if rotation_angle is not None and rotation_axis is None:
            raise ValueError('Must specify rotation_axis')

        if position_angle is not None:
            if inclination is not None:
                rotate = True
                rotation_angle = inclination
                rotation_axis = np.array([np.cos(position_angle),
                                          np.sin(position_angle), 0])
            else:
                raise ValueError('Must specify inclination')

        if inclination is not None and position_angle is None:
            raise ValueError('Must specify position_angle')

        #--- Particle type

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

        #--- Number of pixels

        if number_pixels is None:
            npix = [512, 512]
        else:
            npix = number_pixels

        #--- Get options

        normalize  = self.plot_options['normalize']
        zobserver  = self.plot_options['zobserver']
        dscreen    = self.plot_options['dscreen']
        accelerate = self.plot_options['accelerate']
################################################################################
# TODO: temporary; testing phase
        cross_section = False
        opacity = False
        zslice = None
################################################################################

        #--- Dataframe subsets

        pd = self.particles.loc[self.particles['itype'].isin(itypes)]

        positions        = np.array(pd[['x', 'y', 'z']])
        velocities       = np.array(pd[['vx', 'vy', 'vz']])
        smoothing_length = np.array(pd['h'])
        particle_mass    = np.array(pd['m'])

        if render:
            render_data = np.array(pd[render])

        if vector:
            vector_data = velocities

        #--- Interpolation weights

        weights = _interpolation_weights(
            self.plot_options['density_weighted'], pd,
            self.parameters['hfact']
            )

        #--- Rotate frame

        if rotate:
            print(f'Rotating {rotation_angle*180/np.pi:.0f} deg around ' + \
                  f'[{rotation_axis[0]:.2f}, {rotation_axis[1]:.2f}, {rotation_axis[2]:.2f}]')
            positions, velocities = _rotate_frame(positions,
                                                  velocities,
                                                  rotation_axis,
                                                  rotation_angle)

        #--- Image window range

        if image_range > 0:
            if horizontal_range is not None or vertical_range is not None:
                raise ValueError( 'Cannot set image_range and horizontal_range ' \
                                + '(or vertical_range) at the same time' )
            horizontal_range = [-image_range, image_range]
            vertical_range   = [-image_range, image_range]

        if horizontal_range is None and vertical_range is None:
            range_min = min(positions[:, 0].min(), positions[:, 1].min())
            range_max = max(positions[:, 0].max(), positions[:, 1].max())
            horizontal_range = [range_min, range_max]
            vertical_range   = [range_min, range_max]

        if horizontal_range is None:
            horizontal_range = [positions[:, 0].min(), positions[:, 0].max()]

        if vertical_range is None:
            vertical_range = [positions[:, 1].min(), positions[:, 1].max()]

        extent = horizontal_range + vertical_range

        #--- Interpolate scalar data

        if render:

            print(f'Rendering {render}')

            image_data = scalar_interpolation(
                positions, smoothing_length, weights, render_data, particle_mass,
                horizontal_range, vertical_range, npix, cross_section, zslice,
                opacity, normalize, zobserver, dscreen, accelerate )

        #--- Interpolate vector data

        if vector:

            print(f'Vector field {vector}')

            vector_data = vector_interpolation(
                positions, smoothing_length, weights, vector_data,
                horizontal_range, vertical_range, npix, cross_section, zslice,
                normalize, zobserver, dscreen )

            xvector_data = vector_data[0]
            yvector_data = vector_data[1]

        #--- Render settings

        if render:

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

        #--- Font settings

        mpl.rcParams['font.family'] = self.plot_options['fontfamily']
        mpl.rcParams['font.size']   = self.plot_options['fontsize']

        #--- Figure and axis handles

        if ax is None:
            if newfig:
                plt.figure()
            plt.clf()
            ax = plt.gca()

        #--- Make rendered image

        if render:

            img = ax.imshow(image_data, norm=norm, origin='lower', extent=extent,
                            cmap=cmap)

        #--- Make vector field

        if vector:

            X, Y = np.meshgrid(np.linspace(*horizontal_range, len(xvector_data)),
                               np.linspace(*vertical_range,   len(yvector_data)))

################################################################################
# TODO: temporary; testing phase
# set stride such that number of vector arrows is ~ 15-25 x 15-25
            stride = 25

            vector_color = 'black'
            if render:
                vector_color = 'white'
################################################################################

            q = ax.quiver(X[::stride, ::stride],
                          Y[::stride, ::stride],
                          xvector_data[::stride, ::stride],
                          yvector_data[::stride, ::stride],
                          color=vector_color)
            ax.set_aspect('equal', 'box')

        #--- Axis labels/limits, title

        if not rotate:
            ax.set_xlabel('x')
            ax.set_ylabel('y')

        ax.set_xlim(horizontal_range[0], horizontal_range[1])
        ax.set_ylim(vertical_range[0], vertical_range[1])

        if title is not None:
            ax.set_title(title)

        # TODO: make render_label respond to settings/options
        render_label = r'$\int$ '+ f'{render}' + ' dz'

        #--- Colorbar

        if render:

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

        #--- Handles for Image object

        self._axis = ax
        if render:
            self._image = img
            self._colorbar = cb
        if vector:
            self._quiver = q

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

def _interpolation_weights(density_weighted, particles, hfact):
    '''
    Calculate interpolation weights.
    '''

    if density_weighted:
        return np.array(particles['m'] / particles['h']**2)

    return np.full_like(particles['h'], 1/hfact)

def _rotate_frame(positions, velocities, axis, theta):
    '''
    Rotate around axis.
    '''

    return rotate_vector_arbitrary_axis(positions, axis, theta), \
           rotate_vector_arbitrary_axis(velocities, axis, theta)
