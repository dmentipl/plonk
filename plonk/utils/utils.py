"""Utility functions."""


def is_documented_by(original):
    """Wrap function to add docstring."""

    def wrapper(target):
        target.__doc__ = original.__doc__
        return target

    return wrapper


def time_string(snap, unit_str: str = 'year', float_format: str = '.0f'):
    """Generate time stamp string.

    Parameters
    ----------
    snap
        The Snap object.
    unit_str
        The unit string to pass to Pint. Default is 'year'.
    float_format
        The format for the time float value. Default is '.0f'.

    Examples
    --------
    Generate a list of strings of snapshot time.

    >>> text = [time_string(snap) for snap in snaps]
    """
    time = snap.properties['time'] * snap.units['time'].to(unit_str).m
    return f't = {time:{float_format}} {unit_str}'
