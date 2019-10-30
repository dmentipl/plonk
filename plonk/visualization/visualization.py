"""The Visualization class.

This class contains methods for visualizing smoothed particle
hydrodynamics simulation data.
"""
from typing import Any, Dict, Optional, Tuple

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable
from numpy import ndarray

from .interpolation import scalar_interpolation, vector_interpolation


class Visualization:
    """Visualize scalar and vector smoothed particle hydrodynamics data.

    Visualize SPH data as a particle plot, a rendered image, a vector
    plot, or a combination. The x, y coordinates are the particle
    Cartesian coordinates, corresponding to the x and y axes of the
    plot.

    Pass in options for scalar plots, vector plots, and for
    interpolation via the 'scalar_options', 'vector_options', and
    'interpolation_options' dictionaries.

    Parameters
    ----------
    scalar_data, optional
        The 1d array (N,) of scalar data to visualize.
    vector_data, optional
        The 2d array (N, 2) of vector data to visualize.
    x_coordinate
        The x-position on the particles, where x is the required plot
        x-axis.
    y_coordinate
        The y-position on the particles, where y is the required plot
        y-axis.
    z_coordinate, optional
        The z-position on the particles, where z is the depth-axis for
        the required plot.
    extent
        The range in the x- and y-direction as (xmin, xmax, ymin, ymax).
    particle_mass
        The particle mass for each particle.
    smoothing_length
        The smoothing length for each particle.
    axis, optional
        A matplotlib axis handle.
    scalar_options, optional
        A dictionary of options for scalar plots.
    vector_options, optional
        A dictionary of options for vector plots.
    interpolation_options, optional
        A dictionary of options for interpolation.
    """

    def __init__(
        self,
        *,
        scalar_data: Optional[ndarray] = None,
        vector_data: Optional[ndarray] = None,
        x_coordinate: ndarray,
        y_coordinate: ndarray,
        z_coordinate: Optional[ndarray] = None,
        extent: Tuple[float, float, float, float],
        particle_mass: ndarray,
        smoothing_length: ndarray,
        axis: Optional[Any] = None,
        scalar_options: Dict[str, Any] = None,
        vector_options: Dict[str, Any] = None,
        interpolation_options: Dict[str, Any] = None,
    ):

        self.fig: Any = None
        self.axis: Any = None
        self.scalar: Dict[str, Any] = {
            'image': None,
            'contours': None,
            'colorbar': None,
            'interpolated_data': None,
        }
        self.vector: Dict[str, Any] = {
            'quiver': None,
            'streamplot': None,
            'interpolated_data': None,
        }

        if axis is None:
            self.fig, self.axis = plt.subplots()
        else:
            self.fig = axis.get_figure()
            self.axis = axis

        self.axis.set_xlim(*extent[:2])
        self.axis.set_ylim(*extent[2:])

        if scalar_data is None and vector_data is None:
            self._particle_plot(
                x_coordinate=x_coordinate,
                y_coordinate=y_coordinate,
                smoothing_length=smoothing_length,
                extent=extent,
                axis=self.axis,
            )

        if scalar_data is not None:
            _scalar_plot = self._scalar_plot(
                scalar_data=scalar_data,
                x_coordinate=x_coordinate,
                y_coordinate=y_coordinate,
                z_coordinate=z_coordinate,
                particle_mass=particle_mass,
                smoothing_length=smoothing_length,
                extent=extent,
                fig=self.fig,
                axis=self.axis,
                plot_options=scalar_options,
                interpolation_options=interpolation_options,
            )
            self.scalar['image'] = _scalar_plot[0]
            self.scalar['contours'] = _scalar_plot[1]
            self.scalar['colorbar'] = _scalar_plot[2]
            self.scalar['interpolated_data'] = _scalar_plot[3]

        if vector_data is not None:
            _vector_plot = self._vector_plot(
                vector_data=vector_data,
                x_coordinate=x_coordinate,
                y_coordinate=y_coordinate,
                z_coordinate=z_coordinate,
                particle_mass=particle_mass,
                smoothing_length=smoothing_length,
                extent=extent,
                axis=self.axis,
                plot_options=vector_options,
                interpolation_options=interpolation_options,
            )
            self.vector['quiver'] = _vector_plot[0]
            self.vector['streamplot'] = _vector_plot[1]
            self.vector['interpolated_data'] = _vector_plot[2]

        self.axis.set_aspect('equal')

    def _particle_plot(
        self,
        x_coordinate: ndarray,
        y_coordinate: ndarray,
        smoothing_length: ndarray,
        extent: Tuple[float, float, float, float],
        axis: Any,
    ):

        axis.plot(
            x_coordinate[smoothing_length > 0],
            y_coordinate[smoothing_length > 0],
            'k.',
            markersize=0.5,
        )
        axis.set_xlim(extent[:2])
        axis.set_ylim(extent[2:])
        return

    def _scalar_plot(
        self,
        *,
        scalar_data: ndarray,
        x_coordinate: ndarray,
        y_coordinate: ndarray,
        z_coordinate: Optional[ndarray] = None,
        particle_mass: ndarray,
        smoothing_length: ndarray,
        extent: Tuple[float, float, float, float],
        axis: Any,
        fig: Any,
        plot_options: Dict[str, Any] = None,
        interpolation_options: Dict[str, Any] = None,
    ):

        if plot_options is None:
            plot_options = {}

        plot_render = plot_options.pop('plot_render', True)
        plot_contour = plot_options.pop('plot_contour', False)
        colors = plot_options.pop('contour_color', 'black')
        fmt = plot_options.pop('contour_format', '%.2g')
        _norm = plot_options.pop('norm')
        cmap = plot_options.pop('cmap')

        if len(plot_options) > 0:
            raise ValueError(f'plot_options: {list(plot_options.keys())} not available')

        if interpolation_options is None:
            interpolation_options = {}

        interpolated_data = scalar_interpolation(
            data=scalar_data,
            x_position=x_coordinate,
            y_position=y_coordinate,
            z_position=z_coordinate,
            extent=extent,
            smoothing_length=smoothing_length,
            particle_mass=particle_mass,
            **interpolation_options,
        )

        scalar_image = None
        scalar_colorbar = None
        if plot_render:
            norm = mpl.colors.Normalize()
            if _norm is not None:
                if _norm.lower() in ('linear', 'lin'):
                    norm = mpl.colors.Normalize()
                elif _norm.lower() in ('logarithic', 'logarithm', 'log', 'log10'):
                    norm = mpl.colors.LogNorm()
                else:
                    raise ValueError('Cannot determine normalization for colorbar')

            if cmap is None:
                cmap = 'gist_heat'

            scalar_image = axis.imshow(
                interpolated_data, origin='lower', norm=norm, extent=extent, cmap=cmap
            )

            divider = make_axes_locatable(axis)
            cax = divider.append_axes("right", size="5%", pad=0.05)
            scalar_colorbar = fig.colorbar(scalar_image, cax)

        scalar_contour = None
        if plot_contour:
            n_interp_x, n_interp_y = interpolated_data.shape
            X, Y = np.meshgrid(
                np.linspace(*extent[:2], n_interp_x),
                np.linspace(*extent[2:], n_interp_y),
            )

            scalar_contour = axis.contour(X, Y, interpolated_data, colors=colors)
            scalar_contour.clabel(inline=True, fmt=fmt, fontsize=8)

        return scalar_image, scalar_contour, scalar_colorbar, interpolated_data

    def _vector_plot(
        self,
        *,
        vector_data: ndarray,
        x_coordinate: ndarray,
        y_coordinate: ndarray,
        z_coordinate: Optional[ndarray] = None,
        particle_mass: ndarray,
        smoothing_length: ndarray,
        extent: Tuple[float, float, float, float],
        axis: Any,
        plot_options: Dict[str, Any] = None,
        interpolation_options: Dict[str, Any] = None,
    ):

        if plot_options is None:
            plot_options = {}

        plot_stream = plot_options.pop('plot_stream', False)
        color = plot_options.pop('vector_color', 'black')
        scale = plot_options.pop('vector_scale', None)
        scale_units = plot_options.pop('vector_scale_units', None)
        normalize_vectors = plot_options.pop('normalize_vectors', None)
        number_of_arrows = plot_options.pop('number_of_arrows', (25, 25))

        if len(plot_options) > 0:
            raise ValueError(f'plot_options: {list(plot_options.keys())} not available')

        if interpolation_options is None:
            interpolation_options = {}

        interpolated_data = vector_interpolation(
            x_data=vector_data[:, 0],
            y_data=vector_data[:, 1],
            x_position=x_coordinate,
            y_position=y_coordinate,
            z_position=z_coordinate,
            extent=extent,
            smoothing_length=smoothing_length,
            particle_mass=particle_mass,
            **interpolation_options,
        )

        n_interp_x, n_interp_y = interpolated_data[0].shape
        X, Y = np.meshgrid(
            np.linspace(*extent[:2], n_interp_x), np.linspace(*extent[2:], n_interp_y)
        )
        U, V = interpolated_data[0], interpolated_data[1]

        quiver, streamplot = None, None
        if not plot_stream:
            n_x, n_y = number_of_arrows[0], number_of_arrows[1]
            stride_x = int(n_interp_x / n_x)
            stride_y = int(n_interp_y / n_y)
            X = X[::stride_y, ::stride_x]
            Y = Y[::stride_y, ::stride_x]
            U = U[::stride_y, ::stride_x]
            V = V[::stride_y, ::stride_x]
            if normalize_vectors:
                norm = np.hypot(U, V)
                U /= norm
                V /= norm
            quiver = axis.quiver(
                X, Y, U, V, scale=scale, scale_units=scale_units, color=color
            )

        else:
            streamplot = axis.streamplot(X, Y, U, V, color=color)

        return quiver, streamplot, interpolated_data

    def __repr__(self):
        """Dunder repr method."""
        return '<plonk.Visualization>'
