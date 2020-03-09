"""The utils sub-package."""

from .geometry import coordinate_transform

__all__ = ['coordinate_transform']


def is_documented_by(original):
    """Wrap function to add docstring."""

    def wrapper(target):
        target.__doc__ = original.__doc__
        return target

    return wrapper
