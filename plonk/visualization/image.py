"""
image.py

Daniel Mentiplay, 2019.
"""

import matplotlib.colors as colors
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable

from ..particles import I_GAS
from ..utils import normalize_vector, rotate_vector_arbitrary_axis
from .interpolation import scalar_interpolation, vector_interpolation

_accelerate = False
_colorbar = True
_colormap = 'gist_heat'
_colorscale = 'linear'
_cross_section = False
_density_weighted = False
_dscreen = None
_fontfamily = 'sans-serif'
_fontsize = 12
_normalize = False
_npixx = 512
_npixy = 512
_opacity = False
_slice_thickness = 0.0
_stride = 25
_vector_color = 'black'
_zobserver = None


class Visualisation:
    """
    Plonk visualisation class.
    """

    def __init__(self):
        pass

    def plot(self, dump, render=None, vector=None, **kwargs):
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

        horizontal_range : list of float (len=2), default ``None``
            The range of values for the horizontal (x) axis.
        vertical_range : list of float (len=2), default ``None``
            The range of values for the vertical (y) axis.
        image_range : float, default ``None``
            Specify the horizontal and vertical image range to both
            equal ``image_range``.

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

        stream : bool, default ``None``
            If true, plot vector plots as stream functions.
        cross_section : bool, default ``None``
            If true, plot a cross section rather than column density
            with slice thickness specified by ``slice_thickness``.
        slice_thickness : float, default ``None``
            Thickness of a cross sectional slice.
        number_pixels : list of float (len=2), default ``None``
            The number of pixels in the horizontal and vertical
            directions, like [npixx, npixy], for interpolation. This
            determines the resolution of the image.

        render_scale : str, default ``None``
            Render scale options include: 'log', 'linear'.
        render_min : float, default ``None``
            Minimum value of the rendered quantity.
        render_max : float, default ``None``
            Maximum value of the rendered quantity.
        render_fraction_max : float, default ``None``
            Maximum value of the rendered quantity specified as a
            fraction of the maximum in the data.

        title : str, default ``None``
            Plot title.
        axis : Axes, default ``None``
            Matplotlib Axes object to plot to.
        new_figure : bool, default ``None``
            If true, create new figure using Matplotlib.
        colorbar : bool, default ``None``
            If true, plot colorbar.
        colormap : str, default ``None``
            Specify the colormap.

        Returns
        -------

        Examples
        --------
        Rendering density.
        >>> plonk.plot(dump, render='rho')

        Plotting velocity vectors.
        >>> plonk.plot(dump, vector='velocity')
        """

        # TODO: add options
        # TODO: add docs
        # TODO: choose fluid type: gas, dust1, dust2, ...
        # TODO: add more checks on input
        # TODO: physical units
        # TODO: calculated quantities

        range_options = {
            key: value
            for key, value in kwargs.items()
            if key in ['horizontal_range', 'vertical_range', 'image_range']
        }

        rotation_options = {
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
            key: value for key, value in kwargs.items() if key in ['stream']
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
                'dscreen',
                'normalize',
                'number_pixels',
                'opacity',
                'slice_thickness',
                'zobserver',
            ]
        }

        axis = kwargs.pop('axis', None)
        new_figure = kwargs.pop('new_figure', None)

        if axis is None:
            if new_figure:
                plt.figure()
            plt.clf()
            axis = plt.gca()

        particle_plot = False
        if render is None and vector is None:
            particle_plot = True

        particle_type = kwargs.pop('particle_type', None)

        if particle_type is None:
            particle_type = I_GAS
        else:
            if not isinstance(particle_type, int):
                raise ValueError('particle_type must be int')

        particle_mask = dump.particles['itype'] == particle_type

        positions = dump.particles['xyz'][particle_mask]
        smoothing_length = dump.particles['h'][particle_mask]
        particle_mass = dump.header['massoftype'][particle_type - 1]
        try:
            velocities = dump.particles['vxyz']
        except Exception:
            velocities = None

        density_weighted = kwargs.pop('new_figure', None)
        weights = _interpolation_weights(
            density_weighted,
            smoothing_length,
            dump.header['massoftype'],
            dump.header['hfact'],
        )

        positions, velocities = self._rotate_frame(
            positions, velocities, rotation_options
        )

        horizontal_range, vertical_range = self._set_image_size(
            positions, range_options
        )

        if render:
            self._render_image(
                dump,
                render,
                positions,
                weights,
                particle_mask,
                particle_mass,
                horizontal_range,
                vertical_range,
                interpolation_options,
                render_options,
                figure_options,
                axis,
            )

        if vector:
            self._vector_image(dump, vector, weights, vector_options)

        if particle_plot:
            self._plot_particles(positions, axis)

        self._set_axis_title(
            axis, horizontal_range, vertical_range, figure_options
        )

    def _set_image_size(self, positions, range_options):
        """Set image size."""

        horizontal_range = range_options.pop('horizontal_range', None)
        vertical_range = range_options.pop('vertical_range', None)
        image_range = range_options.pop('image_range', -1)

        if image_range > 0:
            if horizontal_range is not None or vertical_range is not None:
                raise ValueError(
                    'Cannot set image_range and horizontal_range '
                    + '(or vertical_range) at the same time'
                )
            horizontal_range = [-image_range, image_range]
            vertical_range = [-image_range, image_range]

        if horizontal_range is None and vertical_range is None:
            range_min = min(positions[:, 0].min(), positions[:, 1].min())
            range_max = max(positions[:, 0].max(), positions[:, 1].max())
            horizontal_range = [range_min, range_max]
            vertical_range = [range_min, range_max]

        if horizontal_range is None:
            horizontal_range = [positions[:, 0].min(), positions[:, 0].max()]

        if vertical_range is None:
            vertical_range = [positions[:, 1].min(), positions[:, 1].max()]

        return horizontal_range, vertical_range

    def _render_image(
        self,
        dump,
        render,
        positions,
        weights,
        particle_mask,
        particle_mass,
        horizontal_range,
        vertical_range,
        interpolation_options,
        render_options,
        figure_options,
        axis,
    ):

        accelerate = interpolation_options.pop('accelerate', _accelerate)
        cross_section = interpolation_options.pop(
            'cross_section', _cross_section
        )
        dscreen = interpolation_options.pop('dscreen', _dscreen)
        normalize = interpolation_options.pop('normalize', _normalize)
        number_pixels = interpolation_options.pop(
            'number_pixels', [_npixx, _npixy]
        )
        opacity = interpolation_options.pop('opacity', _opacity)
        slice_thickness = interpolation_options.pop(
            'slice_thickness', _slice_thickness
        )
        zobserver = interpolation_options.pop('zobserver', _zobserver)

        if render == 'rho':
            render_data = dump.density_from_smoothing_length()
        elif render == 'x':
            render_data = dump.particles['xyz'][particle_mask][:, 0]
        elif render == 'y':
            render_data = dump.particles['xyz'][particle_mask][:, 1]
        elif render == 'z':
            render_data = dump.particles['xyz'][particle_mask][:, 2]
        elif render == 'vx':
            render_data = dump.particles['vxyz'][particle_mask][:, 0]
        elif render == 'vy':
            render_data = dump.particles['vxyz'][particle_mask][:, 1]
        elif render == 'vz':
            render_data = dump.particles['vxyz'][particle_mask][:, 2]
        else:
            try:
                render_data = dump.particles[render][particle_mask]
                if render_data.ndim != 1:
                    raise ValueError(f'{render} is not 1-dimensional')
            except Exception:
                raise ValueError(
                    f'Cannot determine quantity to render: {render}'
                )

        print(f'Rendering {render} using Splash')

        image_data = scalar_interpolation(
            positions,
            dump.particles['h'][particle_mask],
            weights,
            render_data,
            particle_mass,
            horizontal_range,
            vertical_range,
            number_pixels,
            cross_section,
            slice_thickness,
            opacity,
            normalize,
            zobserver,
            dscreen,
            accelerate,
        )

        self._render_image_matplotlib(
            image_data,
            render,
            horizontal_range,
            vertical_range,
            render_options,
            figure_options,
            axis,
        )

    def _render_image_matplotlib(
        self,
        image_data,
        render,
        horizontal_range,
        vertical_range,
        render_options,
        figure_options,
        axis,
    ):

        render_scale = render_options.pop('render_scale', None)
        render_min = render_options.pop('render_min', None)
        render_max = render_options.pop('render_max', None)
        render_fraction_max = render_options.pop('render_fraction_max', None)
        colormap = render_options.pop('colormap', _colormap)
        colorbar = figure_options.pop('colorbar', _colorbar)

        cmap = _colormap
        if colormap is not None:
            cmap = colormap

        if render_max is None:
            vmax = image_data.max()
        else:
            vmax = render_max

        if render_fraction_max is not None:
            vmax = image_data.max() * render_fraction_max

        if render_min is None:
            vmin = image_data.min()

        if render_scale is None:
            render_scale = _colorscale

        if render_scale == 'log':
            norm = colors.Sym_log_norm(1e-1 * vmax, clip=True)
        elif render_scale == 'linear':
            norm = colors.Normalize(vmin=vmin, vmax=vmax, clip=True)
        else:
            raise ValueError("Unknown color render_scale: " + render_scale)

        extent = horizontal_range + vertical_range

        img = axis.imshow(
            image_data, norm=norm, origin='lower', extent=extent, cmap=cmap
        )

        # TODO: make render_label respond to settings/options
        render_label = r'$\int$ ' + f'{render}' + ' dz'

        cb = None
        if colorbar is not None:
            colorbar_ = colorbar
        else:
            colorbar_ = _colorbar
        if colorbar_:
            divider = make_axes_locatable(axis)
            cax = divider.append_axes("right", size="5%", pad=0.05)
            cb = plt.colorbar(img, cax=cax)
            if render_label:
                cb.set_label(render_label)

    def _vector_image(self, dump, vector, weights, vector_options):

        stream = vector_options.pop('stream', None)

        if stream is None:
            stream = False

        vector_data = velocities

        print(f'Plotting vector field {vector} using Splash')

        vector_data = vector_interpolation(
            positions,
            smoothing_length,
            weights,
            vector_data,
            horizontal_range,
            vertical_range,
            npix,
            cross_section,
            slice_thickness,
            normalize,
            zobserver,
            dscreen,
        )

        xvector_data = vector_data[0]
        yvector_data = vector_data[1]

        X, Y = np.meshgrid(
            np.linspace(*horizontal_range, len(xvector_data)),
            np.linspace(*vertical_range, len(yvector_data)),
        )

        stride = _stride
        vector_color = _vector_color
        if render:
            vector_color = 'white'

        if stream:
            axis.streamplot(
                X[::stride, ::stride],
                Y[::stride, ::stride],
                xvector_data[::stride, ::stride],
                yvector_data[::stride, ::stride],
                color=vector_color,
            )
        else:
            axis.quiver(
                X[::stride, ::stride],
                Y[::stride, ::stride],
                xvector_data[::stride, ::stride],
                yvector_data[::stride, ::stride],
                color=vector_color,
            )

        axis.set_aspect('equal', 'box')

    def _plot_particles(self, positions, axis):
        print('Plotting particles')
        marker_size = 0.01
        axis.scatter(positions[:, 0], positions[:, 1], s=marker_size, c='k')
        axis.set_aspect('equal', 'box')

    def _rotate_frame(self, positions, velocities, rotation_options):

        rotation_axis = rotation_options.pop('rotation_axis', None)
        rotation_angle = rotation_options.pop('rotation_angle', None)
        position_angle = rotation_options.pop('position_angle', None)
        inclination = rotation_options.pop('inclination', None)

        self._rotate = False

        if (rotation_axis is not None or rotation_angle is not None) and (
            position_angle is not None and inclination is not None
        ):
            raise ValueError(
                'Cannot set rotation_axis/rotation_angle and '
                + ' position_angle/inclination at the same time'
            )

        if rotation_axis is not None:
            if rotation_angle is not None:
                self._rotate = True
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
                self._rotate = True
                rotation_angle = inclination
                rotation_axis = np.array(
                    [np.cos(position_angle), np.sin(position_angle), 0]
                )
            else:
                raise ValueError('Must specify inclination')

        if inclination is not None and position_angle is None:
            raise ValueError('Must specify position_angle')

        if rotation_axis is not None and rotation_angle is not None:

            print(
                f'Rotating {rotation_angle*180/np.pi:.0f} deg around '
                f'[{rotation_axis[0]:.2f},'
                f' {rotation_axis[1]:.2f},'
                f' {rotation_axis[2]:.2f}]'
            )

            positions = rotate_vector_arbitrary_axis(
                positions, rotation_axis, rotation_angle
            )

            velocities = rotate_vector_arbitrary_axis(
                velocities, rotation_axis, rotation_angle
            )

        return positions, velocities

    def _set_axis_title(
        self, axis, horizontal_range, vertical_range, figure_options
    ):

        title = figure_options.pop('title', None)

        if not self._rotate:
            axis.set_xlabel('x')
            axis.set_ylabel('y')

        axis.set_xlim(horizontal_range[0], horizontal_range[1])
        axis.set_ylim(vertical_range[0], vertical_range[1])

        if title is not None:
            axis.set_title(title)


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
