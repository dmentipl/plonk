"""Utility functions."""

from typing import Any, Union

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
