"""The Visualization class.

This class contains methods for visualizing smoothed particle
hydrodynamics simulation data.
"""

from __future__ import annotations

from typing import Any, Dict, Optional, Tuple, Union

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable
from numpy import ndarray

from .interpolation import scalar_interpolation, vector_interpolation
from ..snap.snap import Snap, SubSnap

SnapLike = Union[Snap, SubSnap]


class Visualization:
    """Visualize scalar and vector smoothed particle hydrodynamics data.

    Visualize SPH data as:
        - a particle plot,
        - a rendered image,
        - a vector plot,
        - or a combination.
    """

    def __init__(self, snap):
        self.snap = snap
        self.fig: Any = None
        self.axis: Any = None
        self.lines: Any = None
        self.image: Any = None
        self.colorbar: Any = None
        self.contour: Any = None
        self.quiver: Any = None
        self.streamplot: Any = None
        self.extent: Tuple[float, float, float, float] = None
        self.data: Dict[str, ndarray] = {
            'render': None,
            'contour': None,
            'arrow': None,
            'stream': None,
        }

    def plot(
        self,
        *,
        data: Optional[Union[str, ndarray]] = None,
        x: Union[str, ndarray] = 'x',
        y: Union[str, ndarray] = 'y',
        z: Union[str, ndarray] = 'z',
        kind: Optional[str] = None,
        interp: str = 'projection',
        z_slice: float = 0.0,
        extent: Tuple[float, float, float, float],
        axis: Optional[Any] = None,
        **kwargs,
    ) -> Visualization:
        """Plot scalar and vector smoothed particle hydrodynamics data.

        Visualize SPH data as a particle plot, a rendered image, a vector
        plot, or a combination. The x, y coordinates are the particle
        Cartesian coordinates, corresponding to the x and y axes of the
        plot.

        Pass in options for scalar plots, vector plots, and for
        interpolation via the 'scalar_options', 'vector_options', and
        'interpolation_options' dictionaries.

        Parameters
        ----------
        data
            The data to visualize. Can be a string to pass to Snap, or
            a 1d array (N,) of scalar data, or a 2d array (N, 3) of
            vector data. Default is None.
        x
            The x-coordinate for the visualization. Can be a string to
            pass to Snap, or a 1d array (N,). Default is 'x'.
        y
            The y-coordinate for the visualization. Can be a string to
            pass to Snap, or a 1d array (N,).
        z
            The z-coordinate for the visualization. Can be a string to
            pass to Snap, or a 1d array (N,). This is only required for
            cross-section plots.
        kind
            The type of plot.
            - 'particle' : particle plot (default if data is None)
            - 'render' : rendered image (default for scalar data)
            - 'contour' : contour plot (scalar data)
            - 'arrow' : quiver (arrow) plot (default for vector data)
            - 'stream' : stream plot (vector data)
        interp
            The interpolation type.
            - 'projection' : 2d interpolation via projection to xy-plane
            - 'cross_section' : 3d interpolation via cross-section in
              z-direction
        z_slice
            The z-coordinate value of the cross-section slice. Default
            is 0.0.
        extent
            The range in the x and y-coord as (xmin, xmax, ymin, ymax).
        axis
            A matplotlib axis handle.
        """
        if axis is None:
            self.fig, self.axis = plt.subplots()
        else:
            self.fig = axis.get_figure()
            self.axis = axis

        data, x, y, z, kind = _check_input(
            snap=self.snap, data=data, x=x, y=y, z=z, kind=kind
        )

        interpolation_kwargs = ('number_of_pixels', 'density_weighted')
        _kwargs = {
            key: val for key, val in kwargs.items() if key in interpolation_kwargs
        }
        for key in _kwargs:
            kwargs.pop(key)
        interpolated_data = _interpolate(
            snap=self.snap,
            data=data,
            x=x,
            y=y,
            z=z,
            interp=interp,
            z_slice=z_slice,
            extent=extent,
            **_kwargs,
        )

        if kind == 'particle':
            self.lines = particle_plot(
                snap=self.snap, x=x, y=y, extent=extent, axis=self.axis, **kwargs,
            )

        elif kind == 'render':
            show_colorbar = kwargs.pop('show_colorbar', True)
            self.image, self.data['render'] = render_plot(
                data=interpolated_data, extent=extent, axis=self.axis, **kwargs,
            )
            if show_colorbar:
                divider = make_axes_locatable(self.axis)
                cax = divider.append_axes("right", size="5%", pad=0.05)
                self.colorbar = self.fig.colorbar(self.image, cax)

        elif kind == 'contour':
            self.contour, self.data['contour'] = contour_plot(
                data=interpolated_data, extent=extent, axis=self.axis, **kwargs,
            )

        elif kind == 'arrow':
            self.quiver, self.data['arrow'] = arrow_plot(
                data=interpolated_data, extent=extent, axis=self.axis, **kwargs,
            )

        elif kind == 'stream':
            self.streamplot, self.data['stream'] = stream_plot(
                data=interpolated_data, extent=extent, axis=self.axis, **kwargs,
            )

        self.axis.set_xlim(*extent[:2])
        self.axis.set_ylim(*extent[2:])
        self.axis.set_aspect('equal')

        return self

    def __repr__(self):
        """Dunder repr method."""
        return '<plonk.Visualization>'


def _interpolate(
    *,
    snap: SnapLike,
    data: ndarray,
    x: ndarray,
    y: ndarray,
    z: Optional[ndarray] = None,
    interp: 'str',
    z_slice: Optional[float] = None,
    extent: Tuple[float, float, float, float],
    **kwargs,
):
    if interp == 'projection':
        cross_section = None
    elif interp == 'cross_section':
        if z_slice is None:
            z_slice = 0.0
        cross_section = z_slice

    if data.ndim == 1:
        interpolated_data = scalar_interpolation(
            data=data,
            x_coordinate=x,
            y_coordinate=y,
            z_coordinate=z,
            extent=extent,
            smoothing_length=snap['smooth'],
            particle_mass=snap['mass'],
            hfact=snap.properties['hfact'],
            cross_section=cross_section,
            **kwargs,
        )

    elif data.ndim == 2:
        interpolated_data = vector_interpolation(
            x_data=data[:, 0],
            y_data=data[:, 1],
            x_coordinate=x,
            y_coordinate=y,
            z_coordinate=z,
            extent=extent,
            smoothing_length=snap['smooth'],
            particle_mass=snap['mass'],
            hfact=snap.properties['hfact'],
            cross_section=cross_section,
            **kwargs,
        )

    else:
        raise ValueError('data.ndim > 2: cannot determine data')

    return interpolated_data


def particle_plot(
    *,
    snap: SnapLike,
    x: ndarray,
    y: ndarray,
    extent: Tuple[float, float, float, float],
    axis: Any,
    **kwargs,
):
    """Plot the particles.

    Parameters
    ----------
    """
    h: ndarray = snap['smooth']
    mask = (
        (h > 0) & (x > extent[0]) & (x < extent[1]) & (y > extent[2]) & (y < extent[3])
    )
    fmt = kwargs.get('fmt', 'k.')
    lines = axis.plot(x[mask], y[mask], fmt, **kwargs)
    return lines


def render_plot(
    *, data: ndarray, extent: Tuple[float, float, float, float], axis: Any, **kwargs,
):
    """Plot scalar data as a rendered image.

    Visualize scalar SPH data as a rendered image with either projection
    or cross section interpolation.

    Parameters
    ----------
    extent
        The range in the x- and y-direction as (xmin, xmax, ymin, ymax).
    axis
        A matplotlib axis handle.
    interpolation_kwargs
        A dictionary of kwargs to pass to interpolation function.
    """
    try:
        norm = kwargs.pop('norm')
    except KeyError:
        norm = 'linear'
    if norm.lower() in ('linear', 'lin'):
        norm = mpl.colors.Normalize()
    elif norm.lower() in ('logarithic', 'logarithm', 'log', 'log10'):
        norm = mpl.colors.LogNorm()
    else:
        raise ValueError('Cannot determine normalization for colorbar')

    image = axis.imshow(data, origin='lower', extent=extent, norm=norm, **kwargs)

    return image, data


def contour_plot(
    *, data: ndarray, extent: Tuple[float, float, float, float], axis: Any, **kwargs,
):
    """Plot scalar data as a contour plot.

    Visualize scalar SPH data as a contour plot with either projection
    or cross section interpolation.

    Parameters
    ----------
    """
    n_interp_x, n_interp_y = data.shape
    X, Y = np.meshgrid(
        np.linspace(*extent[:2], n_interp_x), np.linspace(*extent[2:], n_interp_y),
    )

    contour = axis.contour(X, Y, data, **kwargs)

    return contour, data


def arrow_plot(
    *, data: ndarray, extent: Tuple[float, float, float, float], axis: Any, **kwargs,
):
    """Plot vector data as a quiver plot.

    Parameters
    ----------
    """
    n_interp_x, n_interp_y = data[0].shape
    X, Y = np.meshgrid(
        np.linspace(*extent[:2], n_interp_x), np.linspace(*extent[2:], n_interp_y)
    )
    U, V = data[0], data[1]

    number_of_arrows = kwargs.pop('number_of_arrows', (25, 25))
    normalize_vectors = kwargs.pop('normalize_vectors', False)

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

    quiver = axis.quiver(X, Y, U, V, **kwargs)

    return quiver, data


def stream_plot(
    *, data: ndarray, extent: Tuple[float, float, float, float], axis: Any, **kwargs,
):
    """Plot vector data as a stream plot.

    Parameters
    ----------
    """
    n_interp_x, n_interp_y = data[0].shape
    X, Y = np.meshgrid(
        np.linspace(*extent[:2], n_interp_x), np.linspace(*extent[2:], n_interp_y)
    )
    U, V = data[0], data[1]

    streamplot = axis.streamplot(X, Y, U, V, **kwargs)

    return streamplot, data


def _get_array_from_input(
    snap: SnapLike, inp: Union[str, ndarray], default: str = None
) -> ndarray:
    if isinstance(inp, str):
        return snap[inp]
    elif isinstance(inp, ndarray):
        return inp
    elif default is not None:
        return snap[default]
    raise ValueError('Cannot determine array to return')


def _check_input(*, snap, data, x, y, z, kind):

    try:
        data = _get_array_from_input(snap, data)
    except ValueError:
        data = None

    x = _get_array_from_input(snap, x)
    y = _get_array_from_input(snap, y)
    z = _get_array_from_input(snap, z)

    if data is not None:
        if data.ndim > 2:
            raise ValueError('Cannot interpret data')
        if kind in ('render', 'contour') and data.ndim != 1:
            raise ValueError('Data is wrong shape for render or contour')
        if kind in ('quiver', 'stream') and data.ndim != 2:
            raise ValueError('Data is wrong shape for quiver or streamplot')
        if kind is None:
            if data.ndim == 1:
                kind = 'render'
            elif data.ndim == 2:
                kind = 'quiver'
    else:
        if kind is None:
            kind = 'particle'
        else:
            raise ValueError(f'No data: can only do particle plot')

    return data, x, y, z, kind
