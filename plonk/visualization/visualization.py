"""
This module contains the Visualization class.

This class contains methods for visualizing smoothed particle
hydrodynamics simulations.
"""

from collections import namedtuple

import matplotlib.colors as colors
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable

from ..core.particles import I_GAS
from ..core.utils import normalize_vector, rotate_vector_arbitrary_axis
from .interpolation import scalar_interpolation, vector_interpolation

PlotOptions = namedtuple(
    'PlotOptions',
    [
        'accelerate',
        'colorbar',
        'colormap',
        'colorscale',
        'cross_section',
        'density_weighted',
        'distance_to_screen',
        'font_family',
        'font_size',
        'normalize',
        'npixx',
        'npixy',
        'opacity',
        'slice_position',
        'stream',
        'stride',
        'vector_color',
        'z_observer',
    ],
)

_DEFAULT_OPTS = PlotOptions(
    accelerate=False,
    colorbar=True,
    colormap='gist_heat',
    colorscale='linear',
    cross_section=False,
    density_weighted=False,
    distance_to_screen=None,
    font_family='sans-serif',
    font_size=12,
    normalize=False,
    npixx=512,
    npixy=512,
    opacity=False,
    slice_position=0.0,
    stream=False,
    stride=25,
    vector_color='black',
    z_observer=None,
)


class Visualization:
    """
    Visualize a dump as a particle plot, a rendered plot, or a
    vector plot.

    Parameters
    ----------
    dump : Dump object
        The plonk.Dump object to visualize.

    render : str, default ``None``
        Scalar quantity to render.

    vector : str, default ``None``
        Vector quantity to be represented as arrows or stream
        function.

    Other Parameters
    ----------------
    particle_type : int, default ``None``
        Particle type to plot, represented as an integer type.

    xrange : list or tuple of float (len=2), default ``None``
        The range of values for the horizontal (x) axis.
    yrange : list or tuple of float (len=2), default ``None``
        The range of values for the vertical (y) axis.
    extent : list or tuple of float (len=4), default ``None``
        Specify the x and y image range as [xmin, xmax, ymin, ymax].

    rotation_axis : list of float (len=3), default ``None``
        A 3-dimensional vector specifying an axis around which to
        rotate the reference frame. Must also specify
        ``rotation_angle``.
    rotation_angle : float, default ``None``
        An angle (radians) to rotate the reference frame around a
        ``rotation_axis``.
    position_angle : float, default ``None``
        An angle (radians) East of North specifying an axis around
        which to incline the reference frame as specified by
        ``inclination``.
    inclination : float, default ``None``
        An angle (radians) of inclination specified relative to a
        ``position_angle``.

    render_scale : str, default 'linear'
        Render scale options include: 'log', 'linear'.
    render_min : float, default None
        Minimum value of the rendered quantity.
    render_max : float, default None
        Maximum value of the rendered quantity.
    render_fraction_max : float, default None
        Maximum value of the rendered quantity specified as a
        fraction of the maximum in the data.

    stream : bool, default False
        If true, plot vector plots as stream functions.
    stride : int, default 25
        Striding through vector interpolation for vector plots.
    vector_color : str, default None
        Color of vector field as arrows or streamfunction.

    axis : Axes, default None
        Matplotlib Axes object to plot to.
    colorbar : bool, default True
        If true, plot colorbar.
    colormap : str, default 'gist_heat'
        Specify the colormap.
    new_figure : bool, default False
        If true, create new figure using Matplotlib.
    title : str, default None
        Plot title.

    accelerate : bool, default False
        Use accelerated interpolation.
    cross_section : bool, default False
        If true, plot a cross section rather than column density
        with slice position specified by ``slice_position``.
    density_weighted : bool, default False
        Use density weighted interpolation.
    distance_to_screen : float, default 0.0
        Distance to screen.
    normalize : bool, default False
        Use normalized interpolation.
    number_pixels : list of float (len=2), default (512, 512)
        The number of pixels in the horizontal and vertical
        directions, like [npixx, npixy], for interpolation. This
        determines the resolution of the image.
    opacity : bool, default False
        Use opacity rendering.
    slice_position : float, 0.0
        Position of the cross sectional slice.
    z_observer : float, default 0.0
        Z position of observer.

    Examples
    --------
    Rendering density.
    >>> viz = plonk.Visualization(dump, render='rho')

    Plotting velocity vectors.
    >>> viz = plonk.Visualization(dump, vector='velocity')
    """

    def __init__(self, dump, render=None, vector=None, **kwargs):

        # TODO: add options
        # TODO: add docs
        # TODO: choose fluid type: gas, dust1, dust2, ...
        # TODO: add more checks on input
        # TODO: physical units
        # TODO: calculated quantities

        self._particles = dump.particles
        self._sinks = dump.sinks
        self._header = dump.header

        render_options = {
            key: value
            for key, value in kwargs.items()
            if key
            in [
                'render_scale',
                'render_min',
                'render_max',
                'render_fraction_max',
                'colormap',
            ]
        }

        vector_options = {
            key: value
            for key, value in kwargs.items()
            if key in ['stream', 'stride', 'vector_color']
        }

        figure_options = {
            key: value
            for key, value in kwargs.items()
            if key in ['title', 'colorbar']
        }

        interpolation_options = {
            key: value
            for key, value in kwargs.items()
            if key
            in [
                'accelerate',
                'cross_section',
                'density_weighted',
                'distance_to_screen',
                'normalize',
                'number_pixels',
                'opacity',
                'slice_position',
                'z_observer',
            ]
        }

        self.axis = kwargs.pop('axis', None)
        new_figure = kwargs.pop('new_figure', None)
        if self.axis is None:
            if new_figure:
                plt.figure()
            plt.clf()
            self.axis = plt.gca()

        self._render = None
        self._vector = None
        self._plot_particles = False
        self._plot_render = False
        self._plot_vector = False
        if render is None:
            self._plot_particles = True
        if render is not None:
            self._render = render
            self._plot_render = True
        if vector is not None:
            self._vector = vector
            self._plot_vector = True

        particle_type = kwargs.pop('particle_type', None)
        if particle_type is None:
            particle_type = I_GAS
        else:
            if not isinstance(particle_type, int):
                raise ValueError('particle_type must be int')
        self._particle_mask = self._particles.itype[()] == particle_type
        self._positions = self._particles.xyz[()][self._particle_mask]
        self._smoothing_length = self._particles.h[()][self._particle_mask]
        self._particle_mass = self._header['massoftype'][particle_type - 1]
        try:
            self._velocities = self._particles.vxyz[()]
        except Exception:
            self._velocities = None

        density_weighted = kwargs.pop('density_weighted', None)
        self._weights = _interpolation_weights(
            density_weighted,
            self._smoothing_length,
            self._particle_mass,
            self._header['hfact'],
        )

        _rotation_options = {
            key: value
            for key, value in kwargs.items()
            if key
            in [
                'rotation_axis',
                'rotation_angle',
                'position_angle',
                'inclination',
            ]
        }
        self._frame_rotation(_rotation_options)

        _xy_range = {
            key: value
            for key, value in kwargs.items()
            if key in ['xrange', 'yrange', 'extent']
        }
        self._set_image_size(_xy_range)

        if self._plot_render:
            image, colorbar = self._render_image(
                interpolation_options, render_options, figure_options
            )

        if self._plot_vector:
            self._vector_image(vector_options, interpolation_options)

        if self._plot_particles:
            self._particle_scatter_plot()

        self._set_axis_title(figure_options)

    def _set_image_size(self, _xy_range):
        """Set image size."""

        self._extent = None

        _xrange = _xy_range.pop('xrange', None)
        _yrange = _xy_range.pop('yrange', None)
        _extent = _xy_range.pop('extent', None)

        if _extent is not None:
            if _xrange is not None or _yrange is not None:
                raise ValueError(
                    'Cannot set extent and xrange/yrange at the same time'
                )
            if len(_extent) == 4:
                self._extent = _extent

        if _xrange is None and _yrange is None:
            _x = min(self._positions[:, 0].min(), self._positions[:, 1].min())
            _y = max(self._positions[:, 0].max(), self._positions[:, 1].max())
        if _xrange is None:
            _x = [self._positions[:, 0].min(), self._positions[:, 0].max()]
        if _yrange is None:
            _y = [self._positions[:, 0].min(), self._positions[:, 0].max()]

        if self._extent is None:
            self._extent = _x + _y

    def _render_image(
        self, interpolation_options, render_options, figure_options
    ):

        accelerate = interpolation_options.pop(
            'accelerate', _DEFAULT_OPTS.accelerate
        )
        cross_section = interpolation_options.pop(
            'cross_section', _DEFAULT_OPTS.cross_section
        )
        distance_to_screen = interpolation_options.pop(
            'distance_to_screen', _DEFAULT_OPTS.distance_to_screen
        )
        normalize = interpolation_options.pop(
            'normalize', _DEFAULT_OPTS.normalize
        )
        number_pixels = interpolation_options.pop(
            'number_pixels', [_DEFAULT_OPTS.npixx, _DEFAULT_OPTS.npixy]
        )
        opacity = interpolation_options.pop('opacity', _DEFAULT_OPTS.opacity)
        slice_position = interpolation_options.pop(
            'slice_position', _DEFAULT_OPTS.slice_position
        )
        z_observer = interpolation_options.pop(
            'z_observer', _DEFAULT_OPTS.z_observer
        )

        if self._render in ['rho', 'dens', 'density']:
            render_data = self._particles.rho
        elif self._render == 'x':
            render_data = self._particles.xyz[()][self._particle_mask][:, 0]
        elif self._render == 'y':
            render_data = self._particles.xyz[()][self._particle_mask][:, 1]
        elif self._render == 'z':
            render_data = self._particles.xyz[()][self._particle_mask][:, 2]
        elif self._render == 'vx':
            render_data = self._particles.vxyz[()][self._particle_mask][:, 0]
        elif self._render == 'vy':
            render_data = self._particles.vxyz[()][self._particle_mask][:, 1]
        elif self._render == 'vz':
            render_data = self._particles.vxyz[()][self._particle_mask][:, 2]
        elif self._render in ['v', 'velocity']:
            render_data = self._particles.extra_quantity('velocity magnitude')[
                self._particle_mask
            ]
        else:
            try:
                render_data = self._particles.arrays[self._render][
                    self._particle_mask
                ]
            except Exception:
                raise ValueError(
                    f'Cannot determine quantity to render: {self._render}'
                )
            if render_data.ndim != 1:
                raise ValueError(f'{self._render} is not 1-dimensional')

        print(f'Rendering {self._render} using Splash')

        image_data = scalar_interpolation(
            self._positions,
            self._smoothing_length,
            self._weights,
            render_data,
            self._particle_mass,
            self._extent[0:2],
            self._extent[2:],
            number_pixels,
            cross_section,
            slice_position,
            opacity,
            normalize,
            z_observer,
            distance_to_screen,
            accelerate,
        )

        image, colorbar = self._render_image_matplotlib(
            image_data, render_options, figure_options
        )

        return image, colorbar

    def _render_image_matplotlib(
        self, image_data, render_options, figure_options
    ):

        render_scale = render_options.pop(
            'render_scale', _DEFAULT_OPTS.colorscale
        )
        render_min = render_options.pop('render_min', None)
        render_max = render_options.pop('render_max', None)
        render_fraction_max = render_options.pop('render_fraction_max', None)
        cmap = render_options.pop('colormap', _DEFAULT_OPTS.colormap)
        colorbar_ = figure_options.pop('colorbar', _DEFAULT_OPTS.colorbar)

        if render_max is None:
            vmax = image_data.max()
        else:
            vmax = render_max

        if render_fraction_max is not None:
            vmax = image_data.max() * render_fraction_max

        if render_min is None:
            vmin = image_data.min()
        else:
            vmin = render_min

        if render_scale == 'log':
            norm = colors.SymLogNorm(1e-1 * vmax, clip=True)
        elif render_scale == 'linear':
            norm = colors.Normalize(vmin=vmin, vmax=vmax, clip=True)
        else:
            raise ValueError("Unknown color render_scale: " + render_scale)

        image = self.axis.imshow(
            image_data,
            norm=norm,
            origin='lower',
            extent=self._extent,
            cmap=cmap,
        )

        # TODO: make render_label respond to settings/options
        render_label = r'$\int$ ' + f'{self._render}' + ' dz'

        colorbar = None
        if colorbar_:
            divider = make_axes_locatable(self.axis)
            cax = divider.append_axes("right", size="5%", pad=0.05)
            colorbar = plt.colorbar(image, cax=cax)
            if render_label:
                colorbar.set_label(render_label)

        return image, colorbar

    def _vector_image(self, vector_options, interpolation_options):

        stream = vector_options.pop('stream', _DEFAULT_OPTS.stream)
        stride = vector_options.pop('stride', _DEFAULT_OPTS.stride)

        cross_section = interpolation_options.pop(
            'cross_section', _DEFAULT_OPTS.cross_section
        )
        distance_to_screen = interpolation_options.pop(
            'distance_to_screen', _DEFAULT_OPTS.distance_to_screen
        )
        normalize = interpolation_options.pop(
            'normalize', _DEFAULT_OPTS.normalize
        )
        number_pixels = interpolation_options.pop(
            'number_pixels', [_DEFAULT_OPTS.npixx, _DEFAULT_OPTS.npixy]
        )
        slice_position = interpolation_options.pop(
            'slice_position', _DEFAULT_OPTS.slice_position
        )
        z_observer = interpolation_options.pop(
            'z_observer', _DEFAULT_OPTS.z_observer
        )

        if self._vector in ['v', 'vel', 'velocity']:
            try:
                vector_data = self._particles.vxyz[()][self._particle_mask]
            except Exception:
                raise ValueError('Velocity not available in dump')
        else:
            try:
                vector_data = self._particles.arrays[self._vector][
                    self._particle_mask
                ]
                if vector_data.ndim != 2 and vector_data.shape[1] != 3:
                    raise ValueError(
                        f'{self._vector} does not have appropriate dimensions'
                    )
            except Exception:
                raise ValueError(
                    f'Cannot determine vector quantity to plot: {self._vector}'
                )

        print(f'Plotting vector field {self._vector} using Splash')

        _xrange, _yrange = self._extent[0:2], self._extent[2:]

        vector_data = vector_interpolation(
            self._positions,
            self._smoothing_length,
            self._weights,
            vector_data,
            _xrange,
            _yrange,
            number_pixels,
            cross_section,
            slice_position,
            normalize,
            z_observer,
            distance_to_screen,
        )

        xvector_data = vector_data[0]
        yvector_data = vector_data[1]

        X, Y = np.meshgrid(
            np.linspace(*_xrange, len(xvector_data)),
            np.linspace(*_yrange, len(yvector_data)),
        )

        vector_color = _DEFAULT_OPTS.vector_color
        if self._render:
            vector_color = 'white'

        if stream:
            self.axis.streamplot(
                X[::stride, ::stride],
                Y[::stride, ::stride],
                xvector_data[::stride, ::stride],
                yvector_data[::stride, ::stride],
                color=vector_color,
            )
        else:
            self.axis.quiver(
                X[::stride, ::stride],
                Y[::stride, ::stride],
                xvector_data[::stride, ::stride],
                yvector_data[::stride, ::stride],
                color=vector_color,
            )

        self.axis.set_aspect('equal', 'box')

    def _particle_scatter_plot(self):
        marker_size = 0.01
        self.axis.scatter(
            self._positions[:, 0], self._positions[:, 1], s=marker_size, c='k'
        )
        self.axis.set_aspect('equal', 'box')

    def _frame_rotation(self, _rotation_options):
        rotation_axis = _rotation_options.get('rotation_axis', None)
        rotation_angle = _rotation_options.get('rotation_angle', None)
        position_angle = _rotation_options.get('position_angle', None)
        inclination = _rotation_options.get('inclination', None)
        self._rotate_frame = False

        if (rotation_axis is not None or rotation_angle is not None) and (
            position_angle is not None and inclination is not None
        ):
            raise ValueError(
                'Cannot set rotation_axis/rotation_angle and '
                + ' position_angle/inclination at the same time'
            )

        if rotation_axis is not None:
            if rotation_angle is not None:
                self._rotate_frame = True
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
                self._rotate_frame = True
                rotation_angle = inclination
                rotation_axis = np.array(
                    [np.cos(position_angle), np.sin(position_angle), 0]
                )
            else:
                raise ValueError('Must specify inclination')

        if inclination is not None and position_angle is None:
            raise ValueError('Must specify position_angle')

        if rotation_axis is not None and rotation_angle is not None:
            self._rotation_axis = rotation_axis
            self._rotation_angle = rotation_angle
            print(
                f'Rotating {rotation_angle*180/np.pi:.0f} deg around '
                f'[{rotation_axis[0]:.2f},'
                f' {rotation_axis[1]:.2f},'
                f' {rotation_axis[2]:.2f}]'
            )
            self._positions = rotate_vector_arbitrary_axis(
                self._positions, rotation_axis, rotation_angle
            )
            self._velocities = rotate_vector_arbitrary_axis(
                self._velocities, rotation_axis, rotation_angle
            )

    def _set_axis_title(self, figure_options):

        title = figure_options.pop('title', None)

        if not self._rotate_frame:
            self.axis.set_xlabel('x')
            self.axis.set_ylabel('y')

        self.axis.set_xlim(self._extent[0], self._extent[1])
        self.axis.set_ylim(self._extent[2], self._extent[3])

        if title is not None:
            self.axis.set_title(title)


def plot_options(**kwargs):
    """
    Create a dictionary with plot options.

    Returns
    -------
    """

    options = dict(_DEFAULT_OPTS._asdict())
    for key, item in kwargs.items():
        if key in options:
            options[key] = item

    return options


def _convert_units():
    """Convert units."""
    # TODO: write this


def _calculate_quantity():
    """Calculate an extra quantity."""
    # TODO: write this


def set_colorbar(cb, vmin=None, vmax=None):
    """Set colorbar limits."""
    if vmin is None:
        vmin = cb._colorbar.vmin
    if vmax is None:
        vmax = cb._colorbar.vmax
    cb._image.set_clim([vmin, vmax])


def _interpolation_weights(density_weighted, smoothing_length, mass, hfact):
    """Calculate interpolation weights."""
    if density_weighted:
        return np.array(mass / smoothing_length ** 2)
    return np.full_like(smoothing_length, 1 / hfact)
