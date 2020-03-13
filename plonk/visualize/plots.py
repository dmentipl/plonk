"""Plot functions for visualization."""

from typing import Any

import matplotlib as mpl
import numpy as np
from numpy import ndarray

from ..snap import SnapLike
from .interpolation import Extent


def particle_plot(
    *, snap: SnapLike, x: ndarray, y: ndarray, extent: Extent, ax: Any, **kwargs,
):
    """Plot particles.

    Parameters
    ----------
    snap
        The Snap object to visualize.
    x
        The x-coordinates for the particle plot.
    y
        The y-coordinates for the particle plot.
    extent
        The range in the x and y-coord as (xmin, xmax, ymin, ymax).
    ax
        A matplotlib Axes handle.
    **kwargs
        Keyword arguments to pass to ax.plot method.

    Returns
    -------
    lines
        A list of matplotlib Line2D objects.
    """
    h: ndarray = snap['smooth']
    mask = (
        (h > 0) & (x > extent[0]) & (x < extent[1]) & (y > extent[2]) & (y < extent[3])
    )
    fmt = kwargs.get('fmt', 'k.')
    lines = ax.plot(x[mask], y[mask], fmt, **kwargs)
    return lines


def render_plot(
    *, interpolated_data: ndarray, extent: Extent, ax: Any, **kwargs,
):
    """Plot 1d interpolated data as a rendered image.

    Parameters
    ----------
    interpolated_data
        The data interpolated to a pixel grid.
    extent
        The range in the x and y-coord as (xmin, xmax, ymin, ymax).
    ax
        A matplotlib Axes handle.
    **kwargs
        Keyword arguments to pass to ax.imshow method.

    Returns
    -------
    image
        A matplotlib AxesImage object.
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

    image = ax.imshow(
        interpolated_data, origin='lower', extent=extent, norm=norm, **kwargs
    )

    return image


def contour_plot(
    *, interpolated_data: ndarray, extent: Extent, ax: Any, **kwargs,
):
    """Plot 1d interpolated data as a contour plot.

    Parameters
    ----------
    interpolated_data
        The data interpolated to a pixel grid.
    extent
        The range in the x and y-coord as (xmin, xmax, ymin, ymax).
    ax
        A matplotlib Axes handle.
    **kwargs
        Keyword arguments to pass to ax.imshow method.

    Returns
    -------
    contour
        A matplotlib QuadContourSet object.
    """
    n_interp_x, n_interp_y = interpolated_data.shape
    X, Y = np.meshgrid(
        np.linspace(*extent[:2], n_interp_x), np.linspace(*extent[2:], n_interp_y),
    )

    contour = ax.contour(X, Y, interpolated_data, **kwargs)

    return contour


def quiver_plot(
    *, interpolated_data: ndarray, extent: Extent, ax: Any, **kwargs,
):
    """Plot 2d interpolated data as a quiver plot.

    Parameters
    ----------
    interpolated_data
        The data interpolated to a pixel grid.
    extent
        The range in the x and y-coord as (xmin, xmax, ymin, ymax).
    ax
        A matplotlib Axes handle.
    **kwargs
        Keyword arguments to pass to ax.imshow method.

    Returns
    -------
    quiver
        A matplotlib Quiver object.
    """
    n_interp_x, n_interp_y = interpolated_data[0].shape
    X, Y = np.meshgrid(
        np.linspace(*extent[:2], n_interp_x), np.linspace(*extent[2:], n_interp_y)
    )
    U, V = interpolated_data[0], interpolated_data[1]

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

    quiver = ax.quiver(X, Y, U, V, **kwargs)

    return quiver


def stream_plot(
    *, interpolated_data: ndarray, extent: Extent, ax: Any, **kwargs,
):
    """Plot 2d interpolated data as a stream plot.

    Parameters
    ----------
    interpolated_data
        The data interpolated to a pixel grid.
    extent
        The range in the x and y-coord as (xmin, xmax, ymin, ymax).
    ax
        A matplotlib Axes handle.
    **kwargs
        Keyword arguments to pass to ax.imshow method.

    Returns
    -------
    streamplot
        A matplotlib StreamplotSet object.
    """
    n_interp_x, n_interp_y = interpolated_data[0].shape
    X, Y = np.meshgrid(
        np.linspace(*extent[:2], n_interp_x), np.linspace(*extent[2:], n_interp_y)
    )
    U, V = interpolated_data[0], interpolated_data[1]

    streamplot = ax.streamplot(X, Y, U, V, **kwargs)

    return streamplot
