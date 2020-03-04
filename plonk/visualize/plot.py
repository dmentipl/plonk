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
    scalar_data: Optional[ndarray] = None,
    vector_data: Optional[ndarray] = None,
    x_coordinate: ndarray,
    y_coordinate: ndarray,
    z_coordinate: Optional[ndarray] = None,
    extent: Tuple[float, float, float, float],
    particle_mass: ndarray,
    smoothing_length: ndarray,
    hfact: float,
    axis: Optional[Any] = None,
    scalar_options: Dict[str, Any] = None,
    vector_options: Dict[str, Any] = None,
    interpolation_options: Dict[str, Any] = None,
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
    scalar_data
        The 1d array (N,) of scalar data to visualize.
    vector_data
        The 2d array (N, 2) of vector data to visualize.
    x_coordinate
        The x-position on the particles, where x is the required plot
        x-axis.
    y_coordinate
        The y-position on the particles, where y is the required plot
        y-axis.
    z_coordinate
        The z-position on the particles, where z is the depth-axis for
        the required plot.
    extent
        The range in the x- and y-direction as (xmin, xmax, ymin, ymax).
    particle_mass
        The particle mass for each particle.
    smoothing_length
        The smoothing length for each particle.
    hfact
        The smoothing length factor.
    axis
        A matplotlib axis handle.
    scalar_options
        A dictionary of options for scalar plots.
    vector_options
        A dictionary of options for vector plots.
    interpolation_options
        A dictionary of options for interpolation.
    """
    return Visualization().plot(
        scalar_data=scalar_data,
        vector_data=vector_data,
        x_coordinate=x_coordinate,
        y_coordinate=y_coordinate,
        z_coordinate=z_coordinate,
        extent=extent,
        particle_mass=particle_mass,
        smoothing_length=smoothing_length,
        hfact=hfact,
        axis=axis,
        scalar_options=scalar_options,
        vector_options=vector_options,
        interpolation_options=interpolation_options,
    )


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
    extent
        The xy extent of the image as (xmin, xmax, ymin, ymax), by
        default None.
    scalar_options
        Options passed to the scalar rendering function, by default
        None.
    interpolation_options
        Options passed to the interpolation function, by default None.
    axis
        The axis handle to plot on.

    Returns
    -------
    Visualization
        The Visualization object for the rendered image.
    """
    if scalar_options is None:
        scalar_options = {}
    if interpolation_options is None:
        interpolation_options = {}

    polar = False
    if scalar_options.get('polar_coordinates'):
        polar = True

    need_z = False
    if interpolation_options.get('cross_section') is not None:
        need_z = True

    if isinstance(quantity, str):
        try:
            scalar_data: ndarray = snap[quantity]
        except ValueError:
            raise ValueError('Cannot determine quantity to render.')
    elif isinstance(quantity, ndarray):
        if quantity.shape[0] != len(snap):
            raise ValueError(
                'Quantity to render must have same length as number of particles'
            )
        scalar_data = quantity
    if scalar_data.ndim > 1:
        raise ValueError('Quantity to render must be 1-dimensional')

    position: ndarray = snap['position']
    smoothing_length: ndarray = snap['smooth']
    particle_mass: ndarray = snap['mass']
    hfact = snap.properties['hfact']

    if extent is None:
        minimum_xy = position[:, :2].min(axis=0)
        maximum_xy = position[:, :2].max(axis=0)
        extent = (minimum_xy[0], maximum_xy[0], minimum_xy[1], maximum_xy[1])
        if polar:
            # extent must be square for polar plots
            extent = (
                max(minimum_xy),
                min(maximum_xy),
                max(minimum_xy),
                min(maximum_xy),
            )

    x_coordinate = position[:, 0]
    y_coordinate = position[:, 1]
    z_coordinate = None
    if need_z:
        z_coordinate = position[:, 2]

    viz = plot(
        scalar_data=scalar_data,
        x_coordinate=x_coordinate,
        y_coordinate=y_coordinate,
        z_coordinate=z_coordinate,
        extent=extent,
        particle_mass=particle_mass,
        smoothing_length=smoothing_length,
        hfact=hfact,
        scalar_options=scalar_options,
        interpolation_options=interpolation_options,
        axis=axis,
    )

    if polar:
        viz.axis.set_aspect('auto')

    return viz


def interpolate(
    snap: SnapLike,
    quantity: Union[str, ndarray],
    extent: Tuple[float, float, float, float],
    scalar_options: Optional[Dict[Any, Any]] = None,
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
    scalar_options
        Options passed to the scalar rendering function, by default
        None.
    interpolation_options
        Options passed to the interpolation function, by default None.

    Returns
    -------
    ndarray
        The interpolated data as an ndarray.
    """
    if scalar_options is None:
        scalar_options = {}
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
