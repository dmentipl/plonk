"""Functions for visualization."""

from __future__ import annotations

from typing import TYPE_CHECKING, Any, List, Union

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import PatchCollection
from matplotlib.patches import Circle

from .._units import Quantity

if TYPE_CHECKING:
    from ..snap.snap import SnapLike


def plot_smoothing_length(
    snap: SnapLike,
    indices: List[int],
    fac: float = 1.0,
    units: Union[str, Quantity] = None,
    x: str = 'x',
    y: str = 'y',
    ax: Any = None,
    **kwargs
) -> Any:
    """Plot smoothing length around particle.

    Also works for plotting the accretion radius on a sink particle.

    Parameters
    ----------
    snap
        The Snap object.
    indices
        The particle indices.
    fac
        Set the circle to be a multiple "fac" of the smoothing length.
    units
        The distance units.
    ax
        The matplotlib Axes to plot on.
    x
        The "x" coordinate.
    y
        The "y" coordinate.
    **kwargs
        Keyword arguments to pass to matplotlib PatchCollection.

    Returns
    -------
    circles
        A list of matplotlib circles.
    """
    px: Quantity = snap[x][indices]
    py: Quantity = snap[y][indices]
    try:
        h: Quantity = snap['smoothing_length'][indices]
    except ValueError:
        # Sinks accretion radius is equivalent to smoothing length
        h = snap['accretion_radius'][indices]
    if units is None:
        arr: Quantity = snap[x]
        units = arr.units
    px = px.to(units).magnitude
    py = py.to(units).magnitude
    h = h.to(units).magnitude

    circles = [Circle((_px, _py), fac * _h) for _px, _py, _h in zip(px, py, h)]
    collection = PatchCollection(circles, **kwargs)
    if ax is None:
        fig, ax = plt.subplots()
        ax.set_aspect('equal')
    ax.add_collection(collection)
    _set_new_lim(ax, px, py, h, fac)

    return collection


def get_extent_from_percentile(
    snap: SnapLike,
    x: str,
    y: str,
    percentile: float = 99,
    x_center_on: float = None,
    y_center_on: float = None,
    edge_factor: float = None,
):
    """Get extent from percentile.

    Parameters
    ----------
    snap
        The Snap object.
    x
        The "x" coordinate.
    y
        The "y" coordinate.
    percentile : optional
        The percentile used in the calculation. Default is 99.
    x_center_on : optional
        Center on some x-value. Default is None.
    y_center_on : optional
        Center on some y-value. Default is None.
    edge_factor : optional
        Add extra spacing to extent. E.g. to add extra 5%, set this
        value to 0.05. Default is None.

    Returns
    -------
    tuple
        The extent of the box as (xmin, xmax, ymin, ymax).
    """
    pl, pr = (100 - percentile) / 2, percentile + (100 - percentile) / 2
    xlim = np.percentile(snap[x], [pl, pr])
    ylim = np.percentile(snap[y], [pl, pr])

    if x_center_on is not None:
        xlim += x_center_on - xlim.mean()
    if y_center_on is not None:
        ylim += y_center_on - ylim.mean()

    if edge_factor is not None:
        dx = xlim[1] - xlim[0]
        dy = ylim[1] - ylim[0]
        xlim += (-dx * edge_factor, dx * edge_factor)
        ylim += (-dy * edge_factor, dy * edge_factor)
        return (xlim[0], xlim[1], ylim[0], ylim[1])

    return (xlim[0], xlim[1], ylim[0], ylim[1])


def _set_new_lim(ax, x, y, h, fac):
    xlim, ylim = ax.get_xlim(), ax.get_ylim()
    x_min = x.min() - fac * h.max()
    x_max = x.max() + fac * h.max()
    y_min = y.min() - fac * h.max()
    y_max = y.max() + fac * h.max()
    _xlim = (min(x_min, xlim[0]), max(x_max, xlim[1]))
    _ylim = (min(y_min, ylim[0]), max(y_max, ylim[1]))
    ax.set_xlim(_xlim)
    ax.set_ylim(_ylim)
