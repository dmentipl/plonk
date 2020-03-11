"""The utils sub-package."""

from .geometry import cartesian_to_polar
from .math import cross, norm

__all__ = ['cartesian_to_polar', 'cross', 'norm']


def is_documented_by(original):
    """Wrap function to add docstring."""

    def wrapper(target):
        target.__doc__ = original.__doc__
        return target

    return wrapper
