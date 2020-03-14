"""Utility functions."""

from typing import Any, Union

import numpy as np

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


def get_extent_from_percentile(snap, percentile=99):
    """Get extent from percentile.

    Calculate a xy-plane square box such that some percentile of
    particles is contained within a sphere inscribed in the box.

    Parameters
    ----------
    snap
        The Snap object.
    percentile : optional
        The percentile used in the calculation. Default is 99.

    Returns
    -------
    tuple
        The extent of the box as (xmin, xmax, ymin, ymax).
    """
    r = np.sqrt(snap['x'] ** 2 + snap['y'] ** 2 + snap['z'] ** 2)
    size = np.percentile(r, percentile)
    return (-size, size, -size, size)
