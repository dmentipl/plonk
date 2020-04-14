"""NumPy functions with Pint support."""

import dask.array as da
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
        if isinstance(x.magnitude, da.Array):
            result_x = x[:, 1] * y[:, 2] - x[:, 2] * y[:, 1]
            result_y = x[:, 2] * y[:, 0] - x[:, 0] * y[:, 2]
            result_z = x[:, 0] * y[:, 1] - x[:, 1] * y[:, 0]
            return da.stack([result_x, result_y, result_z]).T
    if isinstance(x, da.Array):
        result_x = x[:, 1] * y[:, 2] - x[:, 2] * y[:, 1]
        result_y = x[:, 2] * y[:, 0] - x[:, 0] * y[:, 2]
        result_z = x[:, 0] * y[:, 1] - x[:, 1] * y[:, 0]
        return da.stack([result_x, result_y, result_z]).T
    return np.cross(x, y, **kwargs)


def norm(x, **kwargs):
    """Norm of a vector.

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
    return np.linalg.norm(x, **kwargs)


def average(x, weights, **kwargs):
    """Average.

    Parameters
    ----------
    x
        The array (N,) to take the norm of. Can be ndarray or pint
        Quantity.
    weights
        The weights for averaging.
    **kwargs
        Keyword arguments to pass to np.average.

    Returns
    -------
    ndarray
        The average of x.
    """
    if isinstance(x, Quantity):
        return np.average(x.magnitude, weights=weights.magnitude, **kwargs) * x.units
    return np.average(x, weights=weights, **kwargs)
