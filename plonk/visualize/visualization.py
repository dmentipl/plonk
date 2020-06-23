"""Visualize smoothed particle hydrodynamics data."""

from __future__ import annotations

from copy import copy
from typing import Any, Dict, Optional

import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from numpy import ndarray

from .. import Quantity, logger
from ..snap import SnapLike
from . import plots
from .functions import get_extent_from_percentile
from .interpolation import Extent, interpolate

_kind_to_function = {
    'image': plots.imshow,
    'contour': plots.contour,
    'quiver': plots.quiver,
    'streamplot': plots.streamplot,
}


def plot(
    *,
    snap: SnapLike,
    quantity: str,
    x: str = 'x',
    y: str = 'y',
    kind: Optional[str] = None,
    interp: str = 'projection',
    z_slice: float = 0.0,
    extent: Extent = (-1, -1, -1, -1),
    units: Dict[str, Any] = None,
    ax: Optional[Any] = None,
    colorbar_kwargs={},
    **kwargs,
) -> Any:
    """Visualize smoothed particle hydrodynamics data.

    Visualize SPH data by interpolation to a pixel grid. Including:
    an image or a contour plot for scalar data, a quiver (arrow)
    plot or stream plot for vector data.

    Parameters
    ----------
    snap
        The Snap (or SubSnap) object to visualize.
    quantity
        The quantity to visualize. Must be a string to pass to Snap.
    x
        The x-coordinate for the visualization. Must be a string to
        pass to Snap. Default is 'x'.
    y
        The y-coordinate for the visualization. Must be a string to
        pass to Snap. Default is 'y'.
    kind
        The type of plot.

        - 'image' : image plot (default for scalar quantities)
        - 'contour' : contour plot (scalar quantity)
        - 'quiver' : quiver plot (default for vector quantities)
        - 'streamplot' : stream plot (vector quantity)
    interp
        The interpolation type. Default is 'projection'.

        - 'projection' : 2d interpolation via projection to xy-plane
        - 'cross_section' : 3d interpolation via cross-section in
            z-direction
    z_slice
        The z-coordinate value of the cross-section slice. Default
        is 0.0.
    extent
        The range in the x and y-coord as (xmin, xmax, ymin, ymax).
        The default is to set the extent to a box of size such that
        99% of particles are contained within.
    units
        The units of the plot as a dictionary with keys 'quantity',
        'extent', 'projection'. The values are Pint Unit objects.
    ax
        A matplotlib Axes handle.
    colorbar_kwargs
        Keyword arguments to pass to matplotlib colorbar functions.
    **kwargs
        Additional keyword arguments to pass to interpolation and
        matplotlib functions.

    Returns
    -------
    ax
        The matplotlib Axes object.

    Notes
    -----
    Additional parameters passed as keyword arguments will be
    passed to lower level functions as required. E.g. Plonk uses
    matplotlib's imshow for a image plot, so additional arguments
    to imshow can be passed this way.

    See below for additional parameters for interpolation,
    colorbars, quiver plots, etc. All other keyword arguments are
    passed to the appropriate matplotlib function.

    Other Parameters
    ----------------
    number_of_pixels : tuple
        The number of pixels to interpolate particle quantities
        to as a tuple (nx, ny). Default is (512, 512).
    density_weighted : bool
        Whether to density weight the interpolation or not.
        Default is False.
    show_colorbar : bool
        Whether or not to display a colorbar. Default is True.
    number_of_arrows : tuple
        The number of arrows to display by sub-sampling the
        interpolated data. Default is (25, 25).
    normalize_vectors : bool
        Whether to normalize the arrows to all have the same
        length. Default is False.

    Examples
    --------
    Show an image of the surface density in xy-plane.

    >>> plonk.visualize.plot(
    ...     snap=snap,
    ...     quantity='density',
    ... )

    Quiver plot of velocity in xy-plane.

    >>> plonk.visualize.plot(
    ...     snap=snap,
    ...     quantity='velocity',
    ... )

    Set units for the plot.

    >>> units = {
    ...     'quantity': plonk.units('g/cm^3'),
    ...     'extent': plonk.units('au'),
    ...     'projection': plonk.units('cm'),
    ... }

    >>> plonk.visualize.plot(
    ...     snap=snap,
    ...     quantity='density',
    ...     units=units,
    ... )

    You can also use the utility function to generate the units
    dictionary.

    >>> units = plonk.visualize.str_to_units('g/cm^3', 'au', 'cm')
    """
    logger.debug(f'Visualizing "{quantity}" on snap: {snap.file_path.name}')
    _kwargs = copy(kwargs)

    if ax is None:
        fig, ax = plt.subplots()
    else:
        fig = ax.figure

    if extent == (-1, -1, -1, -1):
        extent = get_extent_from_percentile(snap=snap, x=x, y=y)
    if isinstance(extent[0], Quantity):
        extent = (
            (extent[0] / snap.units['length']).to_base_units().magnitude,
            (extent[1] / snap.units['length']).to_base_units().magnitude,
            (extent[2] / snap.units['length']).to_base_units().magnitude,
            (extent[3] / snap.units['length']).to_base_units().magnitude,
        )
    if isinstance(z_slice, Quantity):
        z_slice = (z_slice / snap.units['length']).to_base_units().magnitude

    interpolation_kwargs = ('number_of_pixels', 'density_weighted')
    __kwargs = {key: val for key, val in _kwargs.items() if key in interpolation_kwargs}
    for key in __kwargs:
        _kwargs.pop(key)
    interpolated_data = interpolate(
        snap=snap,
        quantity=quantity,
        x=x,
        y=y,
        interp=interp,
        z_slice=z_slice,
        extent=extent,
        **__kwargs,
    )
    if units is not None:
        interpolated_data, extent = _convert_units(
            snap=snap,
            quantity=quantity,
            interpolated_data=interpolated_data,
            extent=extent,
            units=units,
            interp=interp,
        )

    if kind is None:
        if interpolated_data.ndim == 2:
            kind = 'image'
        elif interpolated_data.ndim == 3:
            kind = 'quiver'

    show_colorbar = _kwargs.pop('show_colorbar', kind == 'image')

    if kind in ('image', 'contour', 'quiver', 'streamplot'):
        plot_object = _kind_to_function[kind](
            interpolated_data=interpolated_data, extent=extent, ax=ax, **_kwargs,
        )

    else:
        raise ValueError('Cannot determine plot type')

    if show_colorbar:
        divider = make_axes_locatable(ax)
        _kwargs = copy(colorbar_kwargs)
        position = _kwargs.pop('position', 'right')
        size = _kwargs.pop('size', '5%')
        pad = _kwargs.pop('pad', '2%')
        if position in ('top', 'bottom'):
            _kwargs.update({'orientation': 'horizontal'})
        cax = divider.append_axes(position=position, size=size, pad=pad)
        fig.colorbar(plot_object, cax, **_kwargs)

    ax.set_xlim(*extent[:2])
    ax.set_ylim(*extent[2:])

    ratio = (extent[1] - extent[0]) / (extent[3] - extent[2])
    if not max(ratio, 1 / ratio) > 10.0:
        ax.set_aspect('equal')

    return ax


def particle_plot(
    *,
    snap: SnapLike,
    x: str = 'x',
    y: str = 'y',
    c: Optional[str] = None,
    s: Optional[str] = None,
    xunit: Any = None,
    yunit: Any = None,
    cunit: Any = None,
    xscale: str = None,
    yscale: str = None,
    ax: Optional[Any] = None,
    colorbar_kwargs={},
    **kwargs,
) -> Any:
    """Particle plots.

    Visualize SPH data by plotting the particles, or a subset of
    the particles.

    Parameters
    ----------
    snap
        The Snap (or SubSnap) object to visualize.
    x
        The x-coordinate for the visualization. Must be a string to
        pass to Snap. Default is 'x'.
    y
        The y-coordinate for the visualization. Must be a string to
        pass to Snap. Default is 'y'.
    c
        The quantity to color the particles. Must be a string to
        pass to Snap.
    s
        The quantity to set the particle size. Must be a string to
        pass to Snap.
    xunit
        The units of the x-coordinate.
    yunit
        The units of the y-coordinate.
    cunit
        The units of the color.
    xscale
        The xscale to pass to the matplotlib Axes method set_xscale.
    yscale
        The yscale to pass to the matplotlib Axes method set_yscale.
    ax
        A matplotlib Axes handle.
    colorbar_kwargs
        Keyword arguments to pass to matplotlib colorbar functions.
    **kwargs
        Additional keyword arguments to pass to matplotlib
        functions.

    Returns
    -------
    ax
        The matplotlib Axes object.
    """
    if ax is None:
        fig, ax = plt.subplots()
    else:
        fig = ax.figure

    _kwargs = {
        'x': x,
        'y': y,
        'c': c,
        's': s,
        'xunit': xunit,
        'yunit': yunit,
        'cunit': cunit,
        'xscale': xscale,
        'yscale': yscale,
        'fig': fig,
        'ax': ax,
        'colorbar_kwargs': colorbar_kwargs,
        **kwargs,
    }
    if c is None and s is None:
        for subsnap in snap.subsnaps_as_list():
            _particle_plot(snap=subsnap, **_kwargs)
    else:
        _particle_plot(snap=snap, **_kwargs)

    return ax


def _particle_plot(
    *,
    snap,
    x='x',
    y='y',
    c=None,
    s=None,
    xunit=None,
    yunit=None,
    cunit=None,
    xscale=None,
    yscale=None,
    fig,
    ax,
    colorbar_kwargs={},
    **kwargs,
):
    _kwargs = copy(kwargs)

    _x: ndarray = snap[x]
    _y: ndarray = snap[y]
    _c: ndarray = snap[c] if c is not None else None
    _s: ndarray = snap[s] if s is not None else None

    if snap._physical_units:
        if xunit is not None:
            _x = _x.to(xunit).magnitude
        else:
            _x = _x.magnitude
        if yunit is not None:
            _y = _y.to(yunit).magnitude
        else:
            _y = _y.magnitude
        if cunit is not None:
            _c = _c.to(cunit).magnitude
        else:
            _c = _c.magnitude

    h: ndarray = snap['smoothing_length']
    mask = h > 0

    _x = _x[mask]
    _y = _y[mask]
    if _c is not None:
        _c = _c[mask]
    if _s is not None:
        _s = _s[mask]
        _s = 100 * _s / _s.max()
        if snap._physical_units:
            _s = _s.magnitude

    show_colorbar = _kwargs.pop('show_colorbar', _c is not None)

    if _s is None and _c is None:
        plots.plot(x=_x, y=_y, ax=ax, **_kwargs)

    else:
        plot_object = plots.scatter(x=_x, y=_y, c=_c, s=_s, ax=ax, **_kwargs)
        if show_colorbar:
            divider = make_axes_locatable(ax)
            _kwargs = copy(colorbar_kwargs)
            position = _kwargs.pop('position', 'right')
            size = _kwargs.pop('size', '5%')
            pad = _kwargs.pop('pad', '2%')
            if position in ('top', 'bottom'):
                _kwargs.update({'orientation': 'horizontal'})
            cax = divider.append_axes(position=position, size=size, pad=pad)
            fig.colorbar(plot_object, cax, **_kwargs)

    if xscale is not None:
        ax.set_xscale(xscale)
    if yscale is not None:
        ax.set_yscale(yscale)

    ratio = (_x.max() - _x.min()) / (_y.max() - _y.min())
    if not max(ratio, 1 / ratio) > 10.0:
        ax.set_aspect('equal')


def _convert_units(
    *,
    snap: SnapLike,
    quantity: str,
    interpolated_data: ndarray,
    extent: Extent,
    units: Dict[str, Any],
    interp: str,
):
    required_keys = {'extent', 'projection', 'quantity'}
    if not set(units) == required_keys:
        raise ValueError(f'units dictionary requires: {required_keys}')
    _quantity_unit = snap.get_array_unit(quantity)
    if interp == 'projection':
        data = (
            (interpolated_data * _quantity_unit * snap.units['length'])
            .to(units['quantity'] * units['projection'])
            .magnitude
        )
    elif interp == 'cross_section':
        data = (interpolated_data * _quantity_unit).to(units['quantity']).magnitude

    new_extent = tuple((extent * snap.units['length']).to(units['extent']).magnitude)

    return data, new_extent
