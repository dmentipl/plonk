"""Plot helper functions."""

from typing import Any, Dict, Optional, Tuple, Union

from numpy import ndarray

from ..snap.snap import Snap, SubSnap
from .interpolation import scalar_interpolation
from .visualization import Visualization

SnapLike = Union[Snap, SubSnap]

_number_of_pixels = (512, 512)


def plot(
    *,
    snap: SnapLike,
    scalar_data: Optional[Union[str, ndarray]] = None,
    vector_data: Optional[Union[str, ndarray]] = None,
    x_coordinate: Optional[Union[str, ndarray]] = None,
    y_coordinate: Optional[Union[str, ndarray]] = None,
    z_coordinate: Optional[Union[str, ndarray]] = None,
    extent: Optional[Tuple[float, float, float, float]] = None,
    axis: Optional[Any] = None,
    scalar_options: Optional[Dict[str, Any]] = None,
    vector_options: Optional[Dict[str, Any]] = None,
    interpolation_options: Optional[Dict[str, Any]] = None,
) -> Visualization:
    """Plot scalar and vector smoothed particle hydrodynamics data.

    Visualize SPH data as a particle plot, a rendered image, a vector
    plot, or a combination. The x, y coordinates are the particle
    Cartesian coordinates, corresponding to the x and y axes of the
    plot.

    Pass in options for scalar plots, vector plots, and for
    interpolation via the 'scalar_options', 'vector_options', and
    'interpolation_options' dictionaries.

    Parameters
    ----------
    snap
        The Snap or SubSnap to visualize.
    scalar_data : optional
        The scalar data to visualize. Can be 1d array (N,) or a string
        resolved by snap[scalar_data].
    vector_data : optional
        The vector data to visualize. Can be 2d array (N, :) or a string
        resolved by snap[scalar_data].
    x_coordinate : optional
        The x-position on the particles, where x is the required plot
        x-axis. Can be 1d array (N,) or a string resolved by
        snap[x_coordinate].
    y_coordinate : optional
        The y-position on the particles, where y is the required plot
        y-axis. Can be 1d array (N,) or a string resolved by
        snap[y_coordinate].
    z_coordinate : optional
        The z-position on the particles, where z is the required plot
        z-axis. Can be 1d array (N,) or a string resolved by
        snap[z_coordinate].
    extent : optional
        The range in the x- and y-direction as (xmin, xmax, ymin, ymax).
    axis
        A matplotlib axis handle.
    scalar_options
        A dictionary of options for scalar plots.
    vector_options
        A dictionary of options for vector plots.
    interpolation_options
        A dictionary of options for interpolation.

    Returns
    -------
    Visualization
        The Visualization object for the rendered image.
    """
    if isinstance(scalar_data, str):
        _scalar_data: ndarray = snap[scalar_data]
    elif isinstance(scalar_data, ndarray):
        _scalar_data = scalar_data
    elif scalar_data is None:
        _scalar_data = None
    else:
        raise ValueError('Cannot determine scalar_data')
    if isinstance(vector_data, str):
        _vector_data: ndarray = snap[vector_data]
    elif isinstance(vector_data, ndarray):
        _vector_data = vector_data
    elif vector_data is None:
        _vector_data = None
    else:
        raise ValueError('Cannot determine vector_data')
    if isinstance(x_coordinate, str):
        _x_coordinate: ndarray = snap[x_coordinate]
    elif isinstance(x_coordinate, ndarray):
        _x_coordinate = x_coordinate
    elif x_coordinate is None:
        _x_coordinate = snap['x']
    else:
        raise ValueError('Cannot determine x_coordinate')
    if isinstance(y_coordinate, str):
        _y_coordinate: ndarray = snap[y_coordinate]
    elif isinstance(y_coordinate, ndarray):
        _y_coordinate = y_coordinate
    elif y_coordinate is None:
        _y_coordinate = snap['y']
    else:
        raise ValueError('Cannot determine y_coordinate')
    if isinstance(z_coordinate, str):
        _z_coordinate: ndarray = snap[z_coordinate]
    elif isinstance(z_coordinate, ndarray):
        _z_coordinate = z_coordinate
    elif z_coordinate is None:
        _z_coordinate = snap['z']
    else:
        raise ValueError('Cannot determine z_coordinate')

    if isinstance(_scalar_data, ndarray):
        if _scalar_data.ndim != 1:
            raise ValueError('Scalar quantity to render must be 1-dimensional')
    if isinstance(_vector_data, ndarray):
        if _vector_data.ndim != 2:
            raise ValueError('Vector quantity to render must be 2-dimensional')
    if _x_coordinate.ndim != 1:
        raise ValueError('x-coordinate to render must be 1-dimensional')
    if _y_coordinate.ndim != 1:
        raise ValueError('y-coordinate to render must be 1-dimensional')
    if _z_coordinate.ndim != 1:
        raise ValueError('z-coordinate to render must be 1-dimensional')

    if scalar_options is None:
        scalar_options = {}

    polar = False
    if scalar_options.get('polar_coordinates'):
        polar = True

    if extent is None:
        if polar:
            # extent must be square for polar plots
            max_min_xy = max(_x_coordinate.min(), _y_coordinate.min())
            min_max_xy = min(_x_coordinate.max(), _y_coordinate.max())
            _extent = (max_min_xy, min_max_xy, max_min_xy, min_max_xy)
        else:
            min_xy = _x_coordinate.min(axis=0), _y_coordinate.min(axis=0)
            max_xy = _x_coordinate.max(axis=0), _y_coordinate.max(axis=0)
            _extent = (min_xy[0], max_xy[0], min_xy[1], max_xy[1])
    else:
        _extent = extent

    viz = Visualization().plot(
        scalar_data=_scalar_data,
        vector_data=_vector_data,
        x_coordinate=_x_coordinate,
        y_coordinate=_y_coordinate,
        z_coordinate=_z_coordinate,
        extent=_extent,
        particle_mass=snap['mass'],
        smoothing_length=snap['smooth'],
        hfact=snap.properties['hfact'],
        axis=axis,
        scalar_options=scalar_options,
        vector_options=vector_options,
        interpolation_options=interpolation_options,
    )

    if polar:
        viz.axis.set_aspect('auto')

    return viz


def render(
    snap: SnapLike,
    quantity: Union[str, ndarray],
    extent: Optional[Tuple[float, float, float, float]] = None,
    scalar_options: Optional[Dict[Any, Any]] = None,
    interpolation_options: Optional[Dict[Any, Any]] = None,
    axis: Optional[Any] = None,
) -> Visualization:
    """Produce a rendered image of a quantity on the snapshot.

    Parameters
    ----------
    snap
        The Snap or SubSnap containing the quantity.
    quantity
        The quantity to render, as a string or ndarray. If quantity is
        a string it must be accessible via snap[quantity]. If the
        quantity is an ndarray it must have the same number of particles
        as the snap.
    extent : optional
        The xy extent of the image as (xmin, xmax, ymin, ymax), by
        default None.
    scalar_options : optional
        Options passed to the scalar rendering function, by default
        None.
    interpolation_options : optional
        Options passed to the interpolation function, by default None.
    axis : optional
        The axis handle to plot on.

    Returns
    -------
    Visualization
        The Visualization object for the rendered image.
    """
    viz = plot(
        snap=snap,
        scalar_data=quantity,
        extent=extent,
        scalar_options=scalar_options,
        interpolation_options=interpolation_options,
        axis=axis,
    )

    return viz


def interpolate(
    snap: SnapLike,
    quantity: Union[str, ndarray],
    extent: Tuple[float, float, float, float],
    interpolation_options: Optional[Dict[Any, Any]] = None,
) -> ndarray:
    """Interpolate a quantity on the snapshot to a pixel grid.

    Parameters
    ----------
    snap
        The Snap or SubSnap containing the quantity.
    quantity
        The quantity to interpolate, as a string or ndarray. If quantity
        is a string it must be accessible via snap[quantity]. If the
        quantity is an ndarray it must have the same number of particles
        as the snap.
    extent
        The xy extent of the image as (xmin, xmax, ymin, ymax).
    interpolation_options : optional
        Options passed to the interpolation function, by default None.

    Returns
    -------
    ndarray
        The interpolated data as an ndarray.
    """
    if interpolation_options is None:
        interpolation_options = {}

    number_of_pixels = interpolation_options.pop('number_of_pixels', _number_of_pixels)
    cross_section = interpolation_options.get('cross_section')
    density_weighted = interpolation_options.get('density_weighted')

    if isinstance(quantity, str):
        try:
            data: ndarray = snap[quantity]
        except ValueError:
            raise ValueError('Cannot determine quantity to render.')
    elif isinstance(quantity, ndarray):
        if quantity.shape[0] != len(snap):
            raise ValueError(
                'Quantity to render must have same length as number of particles'
            )
        data = quantity
    if data.ndim > 1:
        raise ValueError('Quantity to render must be 1-dimensional')

    position: ndarray = snap['position']
    smoothing_length: ndarray = snap['smooth']
    particle_mass: ndarray = snap['mass']
    hfact = snap.properties['hfact']

    x_coordinate = position[:, 0]
    y_coordinate = position[:, 1]
    z_coordinate = None
    if cross_section is not None:
        z_coordinate = position[:, 2]

    interpolated_data = scalar_interpolation(
        data=data,
        x_coordinate=x_coordinate,
        y_coordinate=y_coordinate,
        z_coordinate=z_coordinate,
        extent=extent,
        particle_mass=particle_mass,
        smoothing_length=smoothing_length,
        hfact=hfact,
        number_of_pixels=number_of_pixels,
        cross_section=cross_section,
        density_weighted=density_weighted,
    )

    return interpolated_data
