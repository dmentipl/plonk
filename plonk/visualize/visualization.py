"""Visualize smoothed particle hydrodynamics data."""

from __future__ import annotations

from copy import copy
from typing import TYPE_CHECKING, Any, Dict, Union

import matplotlib.pyplot as plt
import numpy as np
from numpy import ndarray

from .._logging import logger
from .._units import Quantity
from .._units import units as plonk_units
from . import plots
from .functions import get_extent_from_percentile
from .interpolation import Extent, interpolate

if TYPE_CHECKING:
    from ..snap.snap import SnapLike

_kind_to_function = {
    'image': plots.imshow,
    'contour': plots.contour,
    'quiver': plots.quiver,
    'streamplot': plots.streamplot,
}


def plot(
    snap: SnapLike,
    quantity: str,
    *,
    x: str = 'x',
    y: str = 'y',
    kind: str = None,
    interp: str = 'projection',
    z_slice: Union[Quantity, float] = None,
    extent: Extent = (-1, -1, -1, -1),
    units: Dict[str, str] = None,
    ax: Any = None,
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
        The z-coordinate value of the cross-section slice. Can be a
        float or quantity with units of length. Default is 0.0.
    extent
        The range in the x and y-coord as (xmin, xmax, ymin, ymax)
        where xmin, etc. can be floats or quantities with units of
        length. The default is to set the extent to a box of size such
        that 99% of particles are contained within.
    units
        The units of the plot as a dictionary with keys 'quantity',
        'extent', 'projection'. The values are strings representing
        units, e.g. 'g/cm^3'.
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

    >>> plonk.plot(snap=snap, quantity='density')

    Quiver plot of velocity in xy-plane.

    >>> plonk.plot(snap=snap, quantity='velocity')

    Set units for the plot.

    >>> units = {
    ...     'quantity': 'g/cm^3', 'extent': 'au', 'projection': 'cm',
    ... }

    >>> plonk.plot(snap=snap, quantity='density', units=units)
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
    else:
        logger.warning('extent has no units, assuming code units')
    if interp == 'cross_section':
        if z_slice is None:
            z_slice = 0 * plonk_units('cm')
        if isinstance(z_slice, Quantity):
            z_slice = (z_slice / snap.units['length']).to_base_units().magnitude
        else:
            logger.warning('z_slice has no units, assuming code units')
    if units is None:
        _units = {
            'quantity': 1 * snap[quantity].units,
            'extent': 1 * snap['position'].units,
            'projection': 1 * snap['position'].units,
        }
    else:
        qunit = 1 * plonk_units(units.get('quantity', str(snap[quantity].units)))
        eunit = 1 * plonk_units(units.get('extent', str(snap['position'].units)))
        punit = 1 * plonk_units(units.get('projection', str(snap['position'].units)))
        _units = {'quantity': qunit, 'extent': eunit, 'projection': punit}

    interpolation_kwargs = ('number_of_pixels', 'density_weighted')
    __kwargs = {key: val for key, val in _kwargs.items() if key in interpolation_kwargs}
    for key in __kwargs:
        _kwargs.pop(key)

    # Interpolate in code units
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
    # Convert back to physical units
    interpolated_data, extent = _convert_units_for_interpolation(
        snap=snap,
        quantity=quantity,
        interpolated_data=interpolated_data,
        extent=extent,
        units=_units,
        interp=interp,
    )

    if kind is None:
        if interpolated_data.ndim == 2:
            kind = 'image'
        elif interpolated_data.ndim == 3:
            kind = 'quiver'

    show_colorbar = _kwargs.pop('show_colorbar', kind == 'image')

    vmin, vmax = _kwargs.get('vmin', None), _kwargs.get('vmax', None)
    vmin = _convert_units_for_cmap(vmin, 'vmin', _units, interp)
    vmax = _convert_units_for_cmap(vmax, 'vmax', _units, interp)
    if vmin is not None:
        _kwargs['vmin'] = vmin
    if vmax is not None:
        _kwargs['vmax'] = vmax

    if kind in ('image', 'contour', 'quiver', 'streamplot'):
        plot_object = _kind_to_function[kind](
            interpolated_data=interpolated_data, extent=extent, ax=ax, **_kwargs,
        )

    else:
        raise ValueError('Cannot determine plot type')

    if show_colorbar:
        cbar = fig.colorbar(plot_object, ax=ax, **colorbar_kwargs)

        if interp == 'projection':
            qunit = _units['quantity'] * _units['projection']
        elif interp == 'cross_section':
            qunit = _units['quantity']
        if np.allclose(qunit.magnitude, 1.0):
            qunit = qunit.units
        cbar.set_label(f'{quantity} [{qunit:~P}]')

    ax.set_xlim(*extent[:2])
    ax.set_ylim(*extent[2:])

    eunit = _units['extent']
    if np.allclose(eunit.magnitude, 1.0):
        eunit = eunit.units
    ax.set_xlabel(f'{x} [{eunit:~P}]')
    ax.set_ylabel(f'{y} [{eunit:~P}]')

    ratio = (extent[1] - extent[0]) / (extent[3] - extent[2])
    if not max(ratio, 1 / ratio) > 10.0:
        ax.set_aspect('equal')

    return ax


def particle_plot(
    snap: SnapLike,
    *,
    x: str = 'x',
    y: str = 'y',
    c: str = None,
    s: str = None,
    units: Dict[str, str] = None,
    xscale: str = None,
    yscale: str = None,
    ax: Any = None,
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
    units
        The units of the plot as a dictionary with keys 'x', 'y', 'c',
        and 's'. The values are strings representing units, e.g.
        'g/cm^3'.
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

    Examples
    --------
    Show the particles in xy-plane.

    >>> plonk.particle_plot(snap=snap)

    Plot density against x.

    >>> plonk.particle_plot(snap=snap, x='x', y='density')

    Color particles by density in xy-plane.

    >>> plonk.particle_plot(
    ...     snap=snap, x='x', y='y', c='density'
    ... )

    Set units for the plot.

    >>> units = {
    ...     'x': 'au', 'y': 'au', 'cunit': 'g/cm^3',
    ... }

    >>> plonk.particle_plot(
    ...     snap=snap, x='x', y='y', c='density', units=units
    ... )
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
        'units': units,
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
    units=None,
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

    if units is None:
        _units = {
            'x': 1 * snap[x].units,
            'y': 1 * snap[y].units,
        }
        if c is not None:
            _units['c'] = 1 * snap[c].units
        if s is not None:
            _units['s'] = 1 * snap[s].units
    else:
        xunit = units.get('x', str(snap[x].units))
        yunit = units.get('y', str(snap[y].units))
        _units = {
            'x': 1 * plonk_units(xunit),
            'y': 1 * plonk_units(yunit),
        }
        if c is not None:
            cunit = units.get('c', str(snap[c].units))
            _units['c'] = 1 * plonk_units(cunit)
        if s is not None:
            sunit = units.get('s', str(snap[s].units))
            _units['s'] = 1 * plonk_units(sunit)

    _x = _x.to(_units['x']).magnitude
    _y = _y.to(_units['y']).magnitude
    if _c is not None:
        _c = _c.to(_units['c']).magnitude
    if _s is not None:
        _s = _s.to(_units['s']).magnitude

    h: ndarray = snap['smoothing_length'].m
    mask = h > 0

    _x = _x[mask]
    _y = _y[mask]
    if _c is not None:
        _c = _c[mask]
    if _s is not None:
        _s = _s[mask]
        _s = 100 * _s / _s.max()
        _s = _s.magnitude

    show_colorbar = _kwargs.pop('show_colorbar', _c is not None)

    if _s is None and _c is None:
        plots.plot(x=_x, y=_y, ax=ax, **_kwargs)

    else:
        plot_object = plots.scatter(x=_x, y=_y, c=_c, s=_s, ax=ax, **_kwargs)
        if show_colorbar:
            cbar = fig.colorbar(plot_object, ax=ax, **colorbar_kwargs)
            cunit = _units['c']
            if np.allclose(cunit.magnitude, 1.0):
                cunit = cunit.units
            cbar.set_label(f'{c} [{cunit:~P}]')

    if xscale is not None:
        ax.set_xscale(xscale)
    if yscale is not None:
        ax.set_yscale(yscale)

    ratio = (_x.max() - _x.min()) / (_y.max() - _y.min())
    if not max(ratio, 1 / ratio) > 10.0:
        ax.set_aspect('equal')

    xunit, yunit = _units['x'], _units['y']
    if np.allclose(xunit.magnitude, 1.0):
        xunit = xunit.units
    if np.allclose(yunit.magnitude, 1.0):
        yunit = yunit.units
    ax.set_xlabel(f'{x} [{xunit:~P}]')
    ax.set_ylabel(f'{y} [{yunit:~P}]')


def _convert_units_for_interpolation(
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
    quantity_unit = snap.get_array_code_unit(quantity)
    if interp == 'projection':
        proj_unit = units['quantity'] * units['projection']
        data = (interpolated_data * quantity_unit * snap.units['length']).to(
            proj_unit.units
        ).magnitude / proj_unit.magnitude
    elif interp == 'cross_section':
        data = (interpolated_data * quantity_unit).to(
            units['quantity'].units
        ).magnitude / units['quantity'].magnitude

    new_extent = tuple(
        (extent * snap.units['length']).to(units['extent'].units).magnitude
        / units['extent'].magnitude
    )

    return data, new_extent


def _convert_units_for_cmap(vm, name, units, interp):
    if interp == 'projection':
        quantity_unit = units['quantity'] * units['projection']
    elif interp == 'cross_section':
        quantity_unit = units['quantity']
    if vm is not None:
        if isinstance(vm, Quantity):
            if vm.dimensionality != quantity_unit.dimensionality:
                raise ValueError(f'{name} has incorrect units')
            vm = (vm / quantity_unit).to_base_units().magnitude
        else:
            logger.warning(f'{name} has no units, assuming same units as quantity')

    return vm
