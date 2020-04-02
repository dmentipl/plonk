"""Utility functions."""

from typing import Optional

import numpy as np


def is_documented_by(original):
    """Wrap function to add docstring."""

    def wrapper(target):
        target.__doc__ = original.__doc__
        return target

    return wrapper


def time_string(
    snap, unit: str, unit_str: Optional[str] = None, float_format: str = '.0f',
):
    """Generate time stamp string.

    Parameters
    ----------
    snap
        The Snap object.
    unit
        The time unit as a string to pass to Pint. E.g. 'year'.
    unit_str
        The unit string to print. If None, use the value from 'unit'.
        Default is None.
    float_format
        The format for the time float value. Default is '.0f'.

    Examples
    --------
    Generate a list of strings of snapshot time, like
    ['0 yr', '10 yr', ...].

    >>> text = [time_string(snap, 'year', 'yr') for snap in snaps]

    Or, in terms of an orbital time, like '10 orbits'.

    >>> plonk.units.define('binary_orbit = 100 years')
    >>> time_string(snap, 'binary_orbit', 'orbits')
    """
    time = snap.properties['time'].to(unit).magnitude
    if unit_str is None:
        unit_str = unit
    return f't = {time:{float_format}} {unit_str}'


def get_extent_from_percentile(
    snap,
    x: str,
    y: str,
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
