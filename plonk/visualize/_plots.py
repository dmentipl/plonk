"""Plot functions for visualization."""

from typing import Any, Tuple

import matplotlib as mpl
import numpy as np
from numpy import ndarray

from ..snap.snap import SnapLike


def _particle_plot(
    *,
    snap: SnapLike,
    x: ndarray,
    y: ndarray,
    extent: Tuple[float, float, float, float],
    axis: Any,
    **kwargs,
):
    h: ndarray = snap['smooth']
    mask = (
        (h > 0) & (x > extent[0]) & (x < extent[1]) & (y > extent[2]) & (y < extent[3])
    )
    fmt = kwargs.get('fmt', 'k.')
    lines = axis.plot(x[mask], y[mask], fmt, **kwargs)
    return lines


def _render_plot(
    *,
    interpolated_data: ndarray,
    extent: Tuple[float, float, float, float],
    axis: Any,
    **kwargs,
):
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

    image = axis.imshow(
        interpolated_data, origin='lower', extent=extent, norm=norm, **kwargs
    )

    return image


def _contour_plot(
    *,
    interpolated_data: ndarray,
    extent: Tuple[float, float, float, float],
    axis: Any,
    **kwargs,
):
    n_interp_x, n_interp_y = interpolated_data.shape
    X, Y = np.meshgrid(
        np.linspace(*extent[:2], n_interp_x), np.linspace(*extent[2:], n_interp_y),
    )

    contour = axis.contour(X, Y, interpolated_data, **kwargs)

    return contour


def _quiver_plot(
    *,
    interpolated_data: ndarray,
    extent: Tuple[float, float, float, float],
    axis: Any,
    **kwargs,
):
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

    quiver = axis.quiver(X, Y, U, V, **kwargs)

    return quiver


def _stream_plot(
    *,
    interpolated_data: ndarray,
    extent: Tuple[float, float, float, float],
    axis: Any,
    **kwargs,
):
    n_interp_x, n_interp_y = interpolated_data[0].shape
    X, Y = np.meshgrid(
        np.linspace(*extent[:2], n_interp_x), np.linspace(*extent[2:], n_interp_y)
    )
    U, V = interpolated_data[0], interpolated_data[1]

    streamplot = axis.streamplot(X, Y, U, V, **kwargs)

    return streamplot
