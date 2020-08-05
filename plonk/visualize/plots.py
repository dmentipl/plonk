"""Matplotlib wrapper functions for visualization."""

from copy import copy
from typing import Any

import matplotlib as mpl
import numpy as np
from numpy import ndarray

from .._logging import logger
from .interpolation import Extent


def plot(*, x: ndarray, y: ndarray, ax: Any, **kwargs):
    """Plot particles.

    Parameters
    ----------
    x
        The x-coordinates for the particle plot.
    y
        The y-coordinates for the particle plot.
    ax
        A matplotlib Axes handle.
    **kwargs
        Keyword arguments to pass to ax.plot method.

    Returns
    -------
    lines
        A list of matplotlib Line2D objects.
    """
    _kwargs = copy(kwargs)
    linestyle = _kwargs.pop('linestyle', '')
    marker = _kwargs.pop('marker', '.')
    return ax.plot(x, y, linestyle=linestyle, marker=marker, **_kwargs)


def scatter(
    *,
    x: ndarray,
    y: ndarray,
    c: ndarray = None,
    s: ndarray = None,
    n_samples=10_000,
    random_seed=None,
    ax: Any,
    **kwargs,
):
    """Plot particles as scatter plot.

    Parameters
    ----------
    x
        The x-coordinates for the scatter plot.
    y
        The y-coordinates for the scatter plot.
    c
        The quantity to color the particles.
    s
        The quantity to set the particle size.
    n_samples
        The number of samples to take. Default is 10,000.
    random_seed
        The random seed for sampling. Default is None.
    ax
        A matplotlib Axes handle.
    **kwargs
        Keyword arguments to pass to ax.scatter method.

    Returns
    -------
    paths
        A matplotlib PathCollection object.
    """
    if c is None and s is None:
        raise ValueError('Should set size or color')
    if n_samples > 100_000:
        logger.warning('n_samples > 100,000: this may be slow')
    if random_seed is not None:
        np.random.seed(random_seed)

    rand = np.random.choice(len(x), n_samples)

    x = x[rand]
    y = y[rand]
    if c is not None:
        c = c[rand]
    if s is not None:
        s = s[rand]

    _kwargs = copy(kwargs)
    alpha = _kwargs.pop('alpha', 0.5)

    return ax.scatter(x, y, c=c, s=s, alpha=alpha, **_kwargs)


def imshow(*, interpolated_data: ndarray, extent: Extent, ax: Any, **kwargs):
    """Plot 1d interpolated data as an image.

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
    _kwargs = copy(kwargs)
    try:
        norm = _kwargs.pop('norm')
    except KeyError:
        norm = 'linear'
    if norm.lower() in ('linear', 'lin'):
        norm = mpl.colors.Normalize()
    elif norm.lower() in ('logarithic', 'logarithm', 'log', 'log10'):
        norm = mpl.colors.LogNorm()
    else:
        raise ValueError('Cannot determine normalization for colorbar')

    return ax.imshow(
        interpolated_data, origin='lower', extent=extent, norm=norm, **_kwargs
    )


def contour(*, interpolated_data: ndarray, extent: Extent, ax: Any, **kwargs):
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

    return ax.contour(X, Y, interpolated_data, **kwargs)


def quiver(*, interpolated_data: ndarray, extent: Extent, ax: Any, **kwargs):
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

    _kwargs = copy(kwargs)
    number_of_arrows = _kwargs.pop('number_of_arrows', (25, 25))
    normalize_vectors = _kwargs.pop('normalize_vectors', False)

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

    return ax.quiver(X, Y, U, V, **_kwargs)


def streamplot(*, interpolated_data: ndarray, extent: Extent, ax: Any, **kwargs):
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

    return ax.streamplot(X, Y, U, V, **kwargs)
