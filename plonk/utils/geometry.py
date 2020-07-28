"""Geometry, and coordinate transformations."""

from typing import Tuple

import numpy as np
from numpy import ndarray
from scipy.interpolate import RectBivariateSpline

try:
    from skimage import transform
except ImportError:
    transform = None

from .._logging import logger


def cartesian_to_polar(
    interpolated_data_cartesian: ndarray,
    extent_cartesian: Tuple[float, float, float, float],
) -> Tuple[ndarray, Tuple[float, float, float, float]]:
    """Convert interpolated Cartesian pixel grid to polar coordinates.

    Parameters
    ----------
    interpolated_data_cartesian
        The interpolated data on a Cartesian grid.
    extent_cartesian
        The extent in Cartesian space as (xmin, xmax, ymin, ymax). It
        must be square.

    Returns
    -------
    interpolated_data_polar
        The interpolated data on a polar grid (R, phi).
    extent_polar
        The extent on a polar grid (0, Rmax, 0, 2π).
    """
    if transform is None:
        logger.error(
            'cartesian_to_polar requires skimage (scikit-image) which is unavailable\n'
            'try pip install skimage --or-- conda install skimage'
        )
    data, extent = interpolated_data_cartesian, extent_cartesian

    if not np.allclose(extent[1] - extent[0], extent[3] - extent[2]):
        raise ValueError('Bad polar plot: x and y have different scales')

    number_of_pixels = data.shape
    radius_pix = 0.5 * data.shape[0]

    data = transform.warp_polar(data, radius=radius_pix)

    radius = 0.5 * (extent[1] - extent[0])
    extent_polar = (0, radius, 0, 2 * np.pi)

    x_grid = np.linspace(*extent[:2], data.shape[0])
    y_grid = np.linspace(*extent[2:], data.shape[1])
    spl = RectBivariateSpline(x_grid, y_grid, data)
    x_regrid = np.linspace(extent[0], extent[1], number_of_pixels[0])
    y_regrid = np.linspace(extent[2], extent[3], number_of_pixels[1])
    interpolated_data_polar = spl(x_regrid, y_regrid)

    return interpolated_data_polar, extent_polar
