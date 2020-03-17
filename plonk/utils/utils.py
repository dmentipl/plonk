"""Utility functions."""

from typing import Any, Optional, Union

import numpy as np
from numpy import ndarray

from .. import units


def is_documented_by(original):
    """Wrap function to add docstring."""

    def wrapper(target):
        target.__doc__ = original.__doc__
        return target

    return wrapper


def time_string(
    snap,
    unit: Union[str, Any] = units('year'),
    unit_str: str = 'yr',
    float_format: str = '.0f',
):
    """Generate time stamp string.

    Parameters
    ----------
    snap
        The Snap object.
    unit
        The time unit. Can be a string to pass to Pint or a Pint unit.
        Default is plonk.units('year').
    unit_str
        The unit string to print. Default is 'yr'.
    float_format
        The format for the time float value. Default is '.0f'.

    Examples
    --------
    Generate a list of strings of snapshot time.

    >>> text = [time_string(snap) for snap in snaps]
    """
    time = snap.properties['time'] * snap.units['time'].to(unit).magnitude
    return f't = {time:{float_format}} {unit_str}'


def get_extent_from_percentile(
    x: ndarray,
    y: ndarray,
    percentile: float = 99,
    x_center_on: Optional[float] = None,
    y_center_on: Optional[float] = None,
    edge_factor: Optional[float] = None,
):
    """Get extent from percentile.

    Parameters
    ----------
    x
        The "x" coordinate.
    y
        The "x" coordinate.
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
    xlim = np.percentile(x, [pl, pr])
    ylim = np.percentile(y, [pl, pr])

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
