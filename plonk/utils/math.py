"""NumPy functions with Pint support."""

import numpy as np

from .. import Quantity


def cross(x, y, **kwargs):
    """Cross product.

    Parameters
    ----------
    x, y
        The two arrays (N, 3) to take the cross product of. Can be
        ndarray or pint Quantity.
    **kwargs
        Keyword arguments to pass to np.cross.

    Returns
    -------
    ndarray
        The cross product of x and y.
    """
    if isinstance(x, Quantity):
        return np.cross(x.magnitude, y.magnitude, **kwargs) * x.units * y.units
    else:
        return np.cross(x, y, **kwargs)


def norm(x, **kwargs):
    """Cross product.

    Parameters
    ----------
    x
        The arrays (N, 3) to take the norm of. Can be ndarray or pint
        Quantity.
    **kwargs
        Keyword arguments to pass to np.linalg.norm.

    Returns
    -------
    ndarray
        The norm of x.
    """
    if isinstance(x, Quantity):
        return np.linalg.norm(x.magnitude, **kwargs) * x.units
    else:
        return np.linalg.norm(x, **kwargs)
