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
        'cross_section',
        'density_weighted',
        'distance_to_screen',
        'font_family',
        'font_size',
        'normalize',
        'number_pixels',
        'opacity',
        'render_scale',
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
    cross_section=False,
    density_weighted=False,
    distance_to_screen=None,
    font_family='sans-serif',
    font_size=12,
    normalize=False,
    number_pixels=(512, 512),
    opacity=False,
    render_scale='linear',
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
    particle_type : int, default I_GAS
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
    figure : Figure, default None
        Matplotlib Figure to add Axes to.
    font_family : str, default 'sans-serif'
        Font family style for axes title, label, etc.
    font_size : int, default 12
        Font size for axes title, label, etc.
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

        # TODO: physical units
        # TODO: calculated extra quantities

        self._initialized = False

        self._particles = dump.particles
        self._sinks = dump.sinks
        self._header = dump.header
        self._positions = self._particles.xyz
        self._smoothing_length = self._particles.h
        try:
            self._velocities = self._particles.vxyz
        except Exception:
            self._velocities = None

        self._figure_options = {
            key: value
            for key, value in kwargs.items()
            if key in ['colorbar', 'colormap', 'title']
        }

        self._interpolation_options = {
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

        self.axis = kwargs.get('axis', None)
        self.figure = kwargs.get('figure', None)
        if self.axis is None and self.figure is None:
            plt.clf()
            self.figure = plt.gcf()
        if self.axis is None:
            self.axis = self.figure.gca()

        self._render = None
        self._vector = None
        self._plot_particles = False
        self._plot_render = False
        self._plot_vector = False
        if render is None and vector is None:
            self._plot_particles = True
        if render is not None:
            self._render = render
            self._plot_render = True
        if vector is not None:
            self._vector = vector
            self._plot_vector = True

        particle_type = kwargs.get('particle_type', None)
        self.set_particle_type(particle_type)

        density_weighted = kwargs.get('density_weighted', None)
        self._weights = _interpolation_weights(
            density_weighted,
            self._smoothing_length[()][self._particle_mask],
            self._particle_mass,
            self._header['hfact'],
        )

        self._rotation_options = {
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
        self._frame_rotation()

        self._xy_range = {
            key: value
            for key, value in kwargs.items()
            if key in ['xrange', 'yrange', 'extent']
        }
        self.set_image_size()

        self._render_options = {
            key: value
            for key, value in kwargs.items()
            if key
            in [
                'render_scale',
                'render_min',
                'render_max',
                'render_fraction_max',
            ]
        }

        self._vector_options = {
            key: value
            for key, value in kwargs.items()
            if key in ['stream', 'stride', 'vector_color']
        }

        self._make_plot()

        self._initialized = True

    def _make_plot(self):

        if self._plot_render:
            self._render_image()

        if self._plot_vector:
            self._vector_image()

        if self._plot_particles:
            self._particle_scatter_plot()

        self._set_axis()
        self.set_axis_labels()
        self.set_title()

    def set_render_scale(self, render_scale=None):
        """
        Set render scale.

        Parameters
        ----------
        render_scale : str
            A string representing the render color scale, e.g.
            'linear' or 'log'.
        """

        if render_scale is None:
            render_scale = 'linear'

        if render_scale == 'log':
            norm = colors.SymLogNorm(1e-1 * self._vmax, clip=True)
        elif render_scale == 'linear':
            norm = colors.Normalize(vmin=self._vmin, vmax=self._vmax, clip=True)
        else:
            raise ValueError("Unknown color render_scale: " + render_scale)

        self._norm = norm
        self._render_scale = render_scale

        if self._initialized:
            self.image.set_norm(norm)
            self._make_colorbar()

    def set_render_range(self, vmin=None, vmax=None):
        """
        Set render range for colorbar.

        Parameters
        ----------
        vmin : float
            Minimum for the render colorbar.
        vmax : float
            Maximum for the render colorbar.
        """
        if vmin is not None and vmax is not None:
            self.image.set_clim(vmin=vmin, vmax=vmax)
            self._vmin, self._vmax = vmin, vmax
        if vmin is not None:
            self.image.set_clim(vmin=vmin)
            self._vmin = vmin
        if vmax is not None:
            self.image.set_clim(vmax=vmax)
            self._vmax = vmax

    def set_particle_type(self, particle_type):
        """
        Set particle type to visualize.

        Parameters
        ----------
        particle_type : int
            Integer representing the particle type.
        """

        if hasattr(self, '_particle_type'):
            if particle_type == self._particle_type:
                return

        if particle_type is None:
            particle_type = I_GAS
        else:
            if not isinstance(particle_type, int):
                raise ValueError('particle_type must be int')

        self._particle_type = particle_type
        self._particle_mask = self._particles.itype[()] == particle_type
        self._particle_mass = self._header['massoftype'][particle_type - 1]

        if self._initialized:
            self._make_plot()

    def set_image_size(self, extent=None, size=None):
        """
        Set image size.

        Parameters
        ----------
        extent : list or numpy.ndarray
            Extent is the image size: [xmin, xmax, ymin, ymax].
        size : float
            Extent specified by a single value:
            [-size, size, -size, size].
        """

        if extent is not None:
            if hasattr(self, '_extent'):
                if np.all(extent == self._extent):
                    return
            self._extent = extent
            if self._initialized:
                self._make_plot()
            return
        if size is not None:
            if hasattr(self, '_extent'):
                if np.all(size == np.abs(self._extent)):
                    return
            self._extent = [-size, size, -size, size]
            if self._initialized:
                self._make_plot()
            return

        self._extent = None

        _xrange = self._xy_range.get('xrange', None)
        _yrange = self._xy_range.get('yrange', None)
        _extent = self._xy_range.get('extent', None)

        if _extent is not None:
            if _xrange is not None or _yrange is not None:
                raise ValueError(
                    'Cannot set extent and xrange/yrange at the same time'
                )
            if len(_extent) == 4:
                self._extent = _extent

        if self._extent is None:
            if _xrange is None and _yrange is None:
                _min = self._positions[()][self._particle_mask][:, 0:2].min(
                    axis=0
                )
                _max = self._positions[()][self._particle_mask][:, 0:2].max(
                    axis=0
                )
                self._extent = (_min[0], _max[0], _min[1], _max[1])
            else:
                if _xrange is None:
                    _x = (
                        self._positions[()][self._particle_mask][:, 0].min(),
                        self._positions[()][self._particle_mask][:, 0].max(),
                    )
                if _yrange is None:
                    _y = (
                        self._positions[()][self._particle_mask][:, 0].min(),
                        self._positions[()][self._particle_mask][:, 0].max(),
                    )
                if self._extent is None:
                    self._extent = _x + _y

    def _render_image(self):

        if self._render in ['rho', 'dens', 'density']:
            render_data = self._particles.rho[()][self._particle_mask]
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
                render_data = self._particles.arrays[self._render][()][
                    self._particle_mask
                ]
            except Exception:
                raise ValueError(
                    f'Cannot determine quantity to render: {self._render}'
                )
            if render_data.ndim != 1:
                raise ValueError(f'{self._render} is not 1-dimensional')

        accelerate = self._interpolation_options.get(
            'accelerate', _DEFAULT_OPTS.accelerate
        )
        cross_section = self._interpolation_options.get(
            'cross_section', _DEFAULT_OPTS.cross_section
        )
        distance_to_screen = self._interpolation_options.get(
            'distance_to_screen', _DEFAULT_OPTS.distance_to_screen
        )
        normalize = self._interpolation_options.get(
            'normalize', _DEFAULT_OPTS.normalize
        )
        number_pixels = self._interpolation_options.get(
            'number_pixels', _DEFAULT_OPTS.number_pixels
        )
        opacity = self._interpolation_options.get(
            'opacity', _DEFAULT_OPTS.opacity
        )
        slice_position = self._interpolation_options.get(
            'slice_position', _DEFAULT_OPTS.slice_position
        )
        z_observer = self._interpolation_options.get(
            'z_observer', _DEFAULT_OPTS.z_observer
        )

        print(f'Rendering {self._render} using Splash')

        image_data = scalar_interpolation(
            self._positions[()][self._particle_mask],
            self._smoothing_length[()][self._particle_mask],
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

        self._render_image_matplotlib(image_data)

    def _render_image_matplotlib(self, image_data):

        render_scale = self._render_options.get(
            'render_scale', _DEFAULT_OPTS.render_scale
        )
        render_min = self._render_options.get('render_min', None)
        render_max = self._render_options.get('render_max', None)
        render_fraction_max = self._render_options.get(
            'render_fraction_max', None
        )
        cmap = self._render_options.get('colormap', _DEFAULT_OPTS.colormap)
        colorbar_ = self._figure_options.get('colorbar', _DEFAULT_OPTS.colorbar)

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
        self._vmin, self._vmax = vmin, vmax

        self.set_render_scale(render_scale)

        self.image = self.axis.imshow(
            image_data,
            norm=self._norm,
            origin='lower',
            extent=self._extent,
            cmap=cmap,
        )

        if not hasattr(self, 'colorbar'):
            if colorbar_:
                self._make_colorbar()
            else:
                self.colorbar = None
        else:
            self._make_colorbar()

    def _set_render_label(self):
        self._render_label = r'$\int$ ' + f'{self._render}' + ' dz'
        if self._render_scale == 'log':
            self._render_label = ' '.join(('log', self._render_label))

    def _make_colorbar(self):
        if hasattr(self, 'colorbar'):
            self.colorbar.remove()
        divider = make_axes_locatable(self.axis)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        self.colorbar = self.figure.colorbar(self.image, cax=cax)
        self._set_render_label()
        if self._render_label is not None:
            self.colorbar.set_label(self._render_label)

    def _vector_image(self):

        stream = self._vector_options.get('stream', _DEFAULT_OPTS.stream)
        stride = self._vector_options.get('stride', _DEFAULT_OPTS.stride)
        vector_color = self._vector_options.get('stride', _DEFAULT_OPTS.stride)

        cross_section = self._interpolation_options.get(
            'cross_section', _DEFAULT_OPTS.cross_section
        )
        distance_to_screen = self._interpolation_options.get(
            'distance_to_screen', _DEFAULT_OPTS.distance_to_screen
        )
        normalize = self._interpolation_options.get(
            'normalize', _DEFAULT_OPTS.normalize
        )
        number_pixels = self._interpolation_options.get(
            'number_pixels', _DEFAULT_OPTS.number_pixels
        )
        slice_position = self._interpolation_options.get(
            'slice_position', _DEFAULT_OPTS.slice_position
        )
        z_observer = self._interpolation_options.get(
            'z_observer', _DEFAULT_OPTS.z_observer
        )

        if self._vector in ['v', 'vel', 'velocity']:
            try:
                vector_data = self._particles.vxyz[()][self._particle_mask]
            except Exception:
                raise ValueError('Velocity not available in dump')
        else:
            try:
                vector_data = self._particles.arrays[self._vector][()][
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
            self._positions[()][self._particle_mask],
            self._smoothing_length[()][self._particle_mask],
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

        if self._render:
            vector_color = 'white'

        self.stream = None
        self.quiver = None
        if stream:
            self.stream = self.axis.streamplot(
                X[::stride, ::stride],
                Y[::stride, ::stride],
                xvector_data[::stride, ::stride],
                yvector_data[::stride, ::stride],
                color=vector_color,
            )
        else:
            self.quiver = self.axis.quiver(
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
            self._positions[()][self._particle_mask][:, 0],
            self._positions[()][self._particle_mask][:, 1],
            s=marker_size,
            c='k',
        )
        self.axis.set_aspect('equal', 'box')

    def _frame_rotation(self):
        rotation_axis = self._rotation_options.get('rotation_axis', None)
        rotation_angle = self._rotation_options.get('rotation_angle', None)
        position_angle = self._rotation_options.get('position_angle', None)
        inclination = self._rotation_options.get('inclination', None)
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

    def _set_axis(self):

        self.axis.set_xlim(self._extent[0], self._extent[1])
        self.axis.set_ylim(self._extent[2], self._extent[3])

    def set_axis_labels(self, xlabel=None, ylabel=None):
        """
        Set axis labels.

        Parameters
        ----------
        xlabel : str
            Label for the x axis.
        ylabel : str
            Label for the y axis.
        """

        if not self._rotate_frame:
            self.axis.set_xlabel('x')
            self.axis.set_ylabel('y')

        if xlabel is not None:
            self.axis.set_xlabel(xlabel)
        if ylabel is not None:
            self.axis.set_ylabel(ylabel)

    def set_title(self, title=None):
        """
        Set title.

        Parameters
        ----------
        title : str
            Figure title.
        """

        if title is None:
            title = self._figure_options.get('title', None)

        if title is not None:
            self.axis.set_title(title)


def plot_options(**kwargs):
    """
    Create a dictionary with plot options.

    Parameters
    ----------
    **kwargs
        Any values in PlotOptions namedtuple.

    Returns
    -------
    dict
        A dictionary of visualization options.
    """

    options = dict(_DEFAULT_OPTS._asdict())
    for key, item in kwargs.items():
        if key in options:
            options[key] = item

    return options


def _interpolation_weights(density_weighted, smoothing_length, mass, hfact):
    """Calculate interpolation weights."""
    if density_weighted:
        return np.array(mass / smoothing_length ** 2)
    return np.full_like(smoothing_length, 1 / hfact)
