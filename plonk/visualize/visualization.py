"""Visualize smoothed particle hydrodynamics data."""

from __future__ import annotations

import warnings
from copy import copy
from typing import TYPE_CHECKING, Any, Dict, Sequence, Tuple, Union

import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable
from numpy import ndarray

from .._logging import logger
from .._units import Quantity
from .._units import units as plonk_units
from . import plots
from .functions import get_extent_from_percentile
from .interpolation import interpolate

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
    slice_normal: Tuple[float, float, float] = None,
    z_slice: Union[Quantity, float] = None,
    extent: Quantity = None,
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
    slice_normal
        The normal vector to the plane in which to take the
        cross-section slice as a tuple (x, y, z). Default is
        (0, 0, 1).
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
    warnings.warn(
        'In Plonk v0.7.0, plonk.plot will be a particle plot. To make images or'
        'vector plots use plonk.image or plonk.vector.',
        DeprecationWarning,
    )
    logger.debug(f'Visualizing "{quantity}" on snap: {snap.file_path.name}')
    _kwargs = copy(kwargs)

    if kind is None:
        q: Quantity = snap[quantity]
        if q.ndim == 1:
            kind = 'image'
        elif q.ndim == 2:
            kind = 'quiver'
    if kind not in ('image', 'contour', 'quiver', 'streamplot'):
        raise ValueError('Cannot determine plot type')

    if ax is None:
        fig, ax = plt.subplots()
    else:
        fig = ax.figure

    # Interpolate data to plot
    interp_kwargs = {
        key: val
        for key, val in _kwargs.items()
        if key in ('number_of_pixels', 'density_weighted')
    }
    for key in interp_kwargs:
        _kwargs.pop(key)
    _data, _extent, _units = _plot_data(
        snap=snap,
        quantity=quantity,
        x=x,
        y=y,
        interp=interp,
        slice_normal=slice_normal,
        z_slice=z_slice,
        extent=extent,
        units=units,
        **interp_kwargs,
    )

    # Make the actual plot
    _plot_plot(
        interpolated_data=_data,
        extent=_extent,
        names={'quantity': quantity, 'x': x, 'y': y},
        kind=kind,
        interp=interp,
        units=_units,
        ax=ax,
        fig=fig,
        colorbar_kwargs=colorbar_kwargs,
        **_kwargs,
    )

    return ax


def _plot_data(
    snap, quantity, x, y, interp, slice_normal, z_slice, extent, units, **kwargs
):
    if extent is None:
        _extent = get_extent_from_percentile(snap=snap, x=x, y=y)
    else:
        _extent = extent
    if isinstance(_extent[0], Quantity):
        __extent = (
            (_extent[0] / snap.units['length']).to_base_units().magnitude,
            (_extent[1] / snap.units['length']).to_base_units().magnitude,
            (_extent[2] / snap.units['length']).to_base_units().magnitude,
            (_extent[3] / snap.units['length']).to_base_units().magnitude,
        )
    else:
        __extent = _extent
        logger.warning('extent has no units, assuming code units')
    if interp == 'cross_section':
        if slice_normal is None:
            slice_normal = (0, 0, 1)
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

    # Interpolate in code units
    interpolated_data = interpolate(
        snap=snap,
        quantity=quantity,
        x=x,
        y=y,
        interp=interp,
        slice_normal=slice_normal,
        z_slice=z_slice,
        extent=__extent,
        **kwargs,
    )
    # Convert back to physical units
    interpolated_data, ___extent = _convert_units_for_interpolation(
        snap=snap,
        quantity=quantity,
        interpolated_data=interpolated_data,
        extent=__extent,
        units=_units,
        interp=interp,
    )

    return interpolated_data, ___extent, _units


def _plot_plot(
    interpolated_data,
    extent,
    names,
    kind,
    interp,
    units,
    ax,
    fig,
    colorbar_kwargs,
    **kwargs,
):

    show_colorbar = kwargs.pop('show_colorbar', kind == 'image')

    vmin, vmax = kwargs.get('vmin', None), kwargs.get('vmax', None)
    vmin = _convert_units_for_cmap(vmin, 'vmin', units, interp)
    vmax = _convert_units_for_cmap(vmax, 'vmax', units, interp)
    if vmin is not None:
        kwargs['vmin'] = vmin
    if vmax is not None:
        kwargs['vmax'] = vmax

    plot_object = _kind_to_function[kind](
        interpolated_data=interpolated_data, extent=extent, ax=ax, **kwargs,
    )

    ax.set_xlim(*extent[:2])
    ax.set_ylim(*extent[2:])

    eunit = units['extent']
    if np.allclose(eunit.magnitude, 1.0):
        eunit = eunit.units
    ax.set_xlabel(f'{names["x"]} [{eunit:~P}]')
    ax.set_ylabel(f'{names["y"]} [{eunit:~P}]')

    ratio = (extent[1] - extent[0]) / (extent[3] - extent[2])
    if not max(ratio, 1 / ratio) > 10.0:
        ax.set_aspect('equal')

    if show_colorbar:
        divider = make_axes_locatable(ax)
        _kwargs = copy(colorbar_kwargs)
        position = _kwargs.pop('position', 'right')
        size = _kwargs.pop('size', '5%')
        pad = _kwargs.pop('pad', '2%')
        if position in ('top', 'bottom'):
            _kwargs.update({'orientation': 'horizontal'})
        cax = divider.append_axes(position=position, size=size, pad=pad)
        cbar = fig.colorbar(plot_object, cax, **_kwargs)

        if interp == 'projection':
            qunit = units['quantity'] * units['projection']
        elif interp == 'cross_section':
            qunit = units['quantity']
        if np.allclose(qunit.magnitude, 1.0):
            qunit = qunit.units
        cbar.set_label(f'{names["quantity"]} [{qunit:~P}]')


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

    >>> units = {'x': 'au', 'y': 'au', 'c': 'g/cm^3'}

    >>> plonk.particle_plot(
    ...     snap=snap, x='x', y='y', c='density', units=units
    ... )
    """
    warnings.warn(
        'In Plonk v0.7.0, plonk.particle_plot will removed in favor of plonk.plot.',
        DeprecationWarning,
    )
    logger.debug(f'Plotting particles "{x}" vs "{y}"" on snap: {snap.file_path.name}')
    _kwargs = copy(kwargs)

    if ax is None:
        fig, ax = plt.subplots()
    else:
        fig = ax.figure

    if c is None and s is None:
        # If color (c) and size (s) are not required we color each
        # particle type differently
        try:
            subsnaps: Sequence[SnapLike] = snap.subsnaps_as_list()
        except AttributeError:
            # Sinks do not have subsnaps
            subsnaps = [snap]
    else:
        # The subsnaps list is just a list with the original snap
        subsnaps = [snap]

    for subsnap in subsnaps:
        _x, _y, _c, _s, _units = _particle_plot_data(
            snap=subsnap, x=x, y=y, c=c, s=s, units=units
        )
        _particle_plot_plot(
            x=_x,
            y=_y,
            c=_c,
            s=_s,
            units=_units,
            names={'x': x, 'y': y, 'c': c, 's': s},
            xscale=xscale,
            yscale=yscale,
            fig=fig,
            ax=ax,
            colorbar_kwargs=colorbar_kwargs,
            **_kwargs,
        )

    return ax


def _particle_plot_data(snap, x, y, c, s, units):
    _x: Quantity = snap[x]
    _y: Quantity = snap[y]
    _c: Quantity = snap[c] if c is not None else None
    _s: Quantity = snap[s] if s is not None else None

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

    try:
        # ignore accreted particles
        h: ndarray = snap['smoothing_length'].m
        mask = h > 0
        _x = _x[mask]
        _y = _y[mask]
        if _c is not None:
            _c = _c[mask]
        if _s is not None:
            _s = _s[mask]
    except ValueError:
        # sink particles do not have smoothing length
        pass

    if _s is not None:
        _s = 100 * _s / _s.max()

    return _x, _y, _c, _s, _units


def _particle_plot_plot(
    x, y, c, s, units, names, xscale, yscale, fig, ax, colorbar_kwargs, **kwargs
):
    show_colorbar = kwargs.pop('show_colorbar', c is not None)

    if s is None and c is None:
        plots.plot(x=x, y=y, ax=ax, **kwargs)
    else:
        plot_object = plots.scatter(x=x, y=y, c=c, s=s, ax=ax, **kwargs)

    if xscale is not None:
        ax.set_xscale(xscale)
    if yscale is not None:
        ax.set_yscale(yscale)

    ratio = (x.max() - x.min()) / (y.max() - y.min())
    if not max(ratio, 1 / ratio) > 10.0:
        ax.set_aspect('equal')

    xunit, yunit = units['x'], units['y']
    if np.allclose(xunit.magnitude, 1.0):
        xunit = xunit.units
    if np.allclose(yunit.magnitude, 1.0):
        yunit = yunit.units
    ax.set_xlabel(f'{names["x"]} [{xunit:~P}]')
    ax.set_ylabel(f'{names["y"]} [{yunit:~P}]')

    if show_colorbar:
        divider = make_axes_locatable(ax)
        _kwargs = copy(colorbar_kwargs)
        position = _kwargs.pop('position', 'right')
        size = _kwargs.pop('size', '5%')
        pad = _kwargs.pop('pad', '2%')
        if position in ('top', 'bottom'):
            _kwargs.update({'orientation': 'horizontal'})
        cax = divider.append_axes(position=position, size=size, pad=pad)
        cbar = fig.colorbar(plot_object, cax, **_kwargs)

        cunit = units['c']
        if np.allclose(cunit.magnitude, 1.0):
            cunit = cunit.units
        cbar.set_label(f'{names["c"]} [{cunit:~P}]')


def _convert_units_for_interpolation(
    *,
    snap: SnapLike,
    quantity: str,
    interpolated_data: ndarray,
    extent: ndarray,
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
