"""The Visualization class.

This class contains methods for visualizing smoothed particle
hydrodynamics simulation data.
"""
from typing import Any, Dict, Optional, Tuple

import matplotlib.colors as colors
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable
from numpy import ndarray

from .interpolation import scalar_interpolation, vector_interpolation

_HFACT = 1.2


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
        The 2d array (2, N) of vector data to visualize.
    x_coordinate
        The x-position on the particles, where x is the required plot
        x-axis.
    y_coordinate
        The y-position on the particles, where y is the required plot
        y-axis.
    x_range
        The range of x-axis values, like (x_min, x_max).
    y_range
        The range of y-axis values, like (y_min, y_max).
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
        A dictionary of options for interpolation. These get passed to
        the Splash interpolation routines.
    """

    def __init__(
        self,
        *,
        scalar_data: Optional[ndarray] = None,
        vector_data: Optional[ndarray] = None,
        x_coordinate: ndarray,
        y_coordinate: ndarray,
        x_range: Tuple[float, float],
        y_range: Tuple[float, float],
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
            'arrows': None,
            'streamlines': None,
            'interpolated_data': None,
        }

        if axis is None:
            self.fig, self.axis = plt.subplots()
        else:
            self.fig = axis.get_figure()
            self.axis = axis

        self.axis.set_xlim(*x_range)
        self.axis.set_ylim(*y_range)
        self.axis.set_aspect('equal')

        if scalar_data is None and vector_data is None:
            self._particle_plot(
                x_coordinate,
                y_coordinate,
                smoothing_length,
                x_range,
                y_range,
                self.axis,
            )

        if scalar_data is not None:
            _scalar_plot = self._scalar_plot(
                x_coordinate=x_coordinate,
                y_coordinate=y_coordinate,
                scalar_data=scalar_data,
                particle_mass=particle_mass,
                smoothing_length=smoothing_length,
                x_range=x_range,
                y_range=y_range,
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
                x_coordinate=x_coordinate,
                y_coordinate=y_coordinate,
                vector_data=vector_data,
                particle_mass=particle_mass,
                smoothing_length=smoothing_length,
                x_range=x_range,
                y_range=y_range,
                axis=self.axis,
                plot_options=vector_options,
                interpolation_options=interpolation_options,
            )
            self.vector['arrows'] = _vector_plot[0]
            self.vector['streamlines'] = _vector_plot[1]
            self.vector['interpolated_data'] = _vector_plot[2]

    def _particle_plot(
        self,
        x_coordinate: ndarray,
        y_coordinate: ndarray,
        smoothing_length: ndarray,
        x_range: Tuple[float, float],
        y_range: Tuple[float, float],
        axis: Any,
    ):

        axis.plot(
            x_coordinate[smoothing_length > 0],
            y_coordinate[smoothing_length > 0],
            'k.',
            markersize=0.5,
        )
        return

    def _scalar_plot(
        self,
        x_coordinate: ndarray,
        y_coordinate: ndarray,
        scalar_data: ndarray,
        particle_mass: ndarray,
        smoothing_length: ndarray,
        x_range: Tuple[float, float],
        y_range: Tuple[float, float],
        axis: Any,
        fig: Any,
        plot_options: Dict[str, Any] = None,
        interpolation_options: Dict[str, Any] = None,
    ):

        if plot_options is None:
            plot_options = {}

        plot_render = plot_options.get('plot_render', True)
        plot_contour = plot_options.get('plot_contour', False)

        if interpolation_options is None:
            interpolation_options = {}

        interpolated_data = scalar_interpolation(
            x_coordinate=x_coordinate,
            y_coordinate=y_coordinate,
            z_coordinate=np.zeros(x_coordinate.shape),
            smoothing_length=smoothing_length,
            smoothing_length_factor=_HFACT,
            scalar_data=scalar_data,
            particle_mass=particle_mass,
            x_range=x_range,
            y_range=y_range,
            **interpolation_options,
        )

        scalar_image = None
        scalar_colorbar = None
        if plot_render:
            norm = colors.Normalize()
            _norm = plot_options.get('norm')
            if _norm is not None:
                if _norm.lower() in ('linear', 'lin'):
                    norm = colors.Normalize()
                elif _norm.lower() in ('logarithic', 'logarithm', 'log', 'log10'):
                    norm = colors.LogNorm()
                else:
                    raise ValueError('Cannot determine normalization for colorbar')

            cmap = plot_options.get('cmap')
            if cmap is None:
                cmap = 'gist_heat'

            scalar_image = axis.imshow(
                interpolated_data,
                origin='lower',
                norm=norm,
                extent=x_range + y_range,
                cmap=cmap,
            )

            divider = make_axes_locatable(axis)
            cax = divider.append_axes("right", size="5%", pad=0.05)
            scalar_colorbar = fig.colorbar(scalar_image, cax)

        scalar_contour = None
        if plot_contour:
            contour_color = plot_options.get('contour_color', 'black')
            contour_format = plot_options.get('contour_format', '%.2g')
            n_interp_x, n_interp_y = interpolated_data.shape
            X, Y = np.meshgrid(
                np.linspace(*x_range, n_interp_x), np.linspace(*y_range, n_interp_y)
            )

            scalar_contour = axis.contour(X, Y, interpolated_data, colors=contour_color)
            scalar_contour.clabel(inline=True, fmt=contour_format, fontsize=8)

        return scalar_image, scalar_contour, scalar_colorbar, interpolated_data

    def _vector_plot(
        self,
        x_coordinate: ndarray,
        y_coordinate: ndarray,
        vector_data: ndarray,
        particle_mass: ndarray,
        smoothing_length: ndarray,
        x_range: Tuple[float, float],
        y_range: Tuple[float, float],
        axis: Any,
        plot_options: Dict[str, Any] = None,
        interpolation_options: Dict[str, Any] = None,
    ):

        if plot_options is None:
            plot_options = {}

        plot_stream = plot_options.get('plot_stream', False)
        vector_color = plot_options.get('vector_color', 'black')

        if interpolation_options is None:
            interpolation_options = {}

        interpolated_data = vector_interpolation(
            x_coordinate=x_coordinate,
            y_coordinate=y_coordinate,
            z_coordinate=np.zeros(x_coordinate.shape),
            smoothing_length=smoothing_length,
            smoothing_length_factor=_HFACT,
            x_vector_data=vector_data[:, 0],
            y_vector_data=vector_data[:, 1],
            particle_mass=particle_mass,
            x_range=x_range,
            y_range=y_range,
            **interpolation_options,
        )

        n_interp_x, n_interp_y = interpolated_data[0].shape
        X, Y = np.meshgrid(
            np.linspace(*x_range, n_interp_x), np.linspace(*y_range, n_interp_y)
        )
        U, V = interpolated_data[0], interpolated_data[1]

        arrows, streamlines = None, None
        if not plot_stream:
            number_of_arrows = plot_options.get('number_of_arrows', (25, 25))
            n_x, n_y = number_of_arrows[0], number_of_arrows[1]
            stride_x = int(n_interp_x / n_x)
            stride_y = int(n_interp_y / n_y)
            X = X[::stride_x, ::stride_y]
            Y = Y[::stride_x, ::stride_y]
            U = U[::stride_x, ::stride_y]
            V = V[::stride_x, ::stride_y]
            arrows = axis.quiver(X, Y, U, V, color=vector_color)

        else:
            streamlines = axis.streamplot(X, Y, U, V, color=vector_color)

        return arrows, streamlines, interpolated_data

    def __repr__(self):
        """Dunder repr method."""
        return '<plonk.Visualization>'
