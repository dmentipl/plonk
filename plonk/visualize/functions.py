"""Functions for visualization."""

from __future__ import annotations

from typing import TYPE_CHECKING, Any, List

import matplotlib.pyplot as plt
import numpy as np

from .._units import Quantity

if TYPE_CHECKING:
    from ..snap.snap import SnapLike


def plot_smoothing_length(
    snap: SnapLike,
    indices: List[int],
    fac: float = 2.0,
    units: Quantity = None,
    x: str = 'x',
    y: str = 'y',
    ax: Any = None,
    **kwargs
) -> Any:
    """Plot smoothing length around particle.

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
        Keyword arguments to pass to plt.Circle.

    Returns
    -------
    circles
        A list of matplotlib circles.
    """
    if ax is None:
        fig, ax = plt.subplots()
    if units is None:
        units = snap[x].units
    circles = list()
    for index in indices:
        px: Quantity = snap[x][index]
        py: Quantity = snap[y][index]
        h: Quantity = snap['smoothing_length'][index]
        px = px.to(units).magnitude
        py = py.to(units).magnitude
        h = h.to(units).magnitude
        pos = (px, py)
        circles.append(plt.Circle(pos, fac * h, **kwargs))
    for c in circles:
        ax.add_artist(c)
    return circles


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
