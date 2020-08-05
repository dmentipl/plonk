"""Utility functions."""


def is_documented_by(original):
    """Wrap function to add docstring."""

    def wrapper(target):
        target.__doc__ = original.__doc__
        return target

    return wrapper


def time_string(
    snap, unit: str, unit_str: str = None, float_format: str = '.0f',
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
