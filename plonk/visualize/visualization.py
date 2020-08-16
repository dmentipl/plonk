"""Visualize smoothed particle hydrodynamics data."""

from __future__ import annotations

from copy import copy
from typing import TYPE_CHECKING, Any, Dict, Sequence, Tuple, Union

import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable
from numpy import ndarray

from .._logging import logger
from .._units import Quantity
from .._units import units as plonk_units
from ..utils.utils import pretty_array_name
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


def image(
    snap: SnapLike,
    quantity: str,
    *,
    x: str = 'x',
    y: str = 'y',
    interp: str = 'projection',
    slice_normal: Tuple[float, float, float] = None,
    slice_offset: Union[Quantity, float] = None,
    extent: Quantity = None,
    units: Dict[str, str] = None,
    ax: Any = None,
    ax_kwargs={},
    colorbar_kwargs={},
    **kwargs,
) -> Any:
    """Visualize scalar SPH data as an image.

    Visualize scalar smoothed particle hydrodynamics data by
    interpolation to a pixel grid.

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
    interp
        The interpolation type. Default is 'projection'.

        - 'projection' : 2d interpolation via projection to xy-plane
        - 'slice' : 3d interpolation via cross-section slice.
    slice_normal
        The normal vector to the plane in which to take the
        cross-section slice as a tuple (x, y, z). Default is
        (0, 0, 1).
    slice_offset
        The offset of the cross-section slice. Default is 0.0.
    extent
        The range in the x and y-coord as (xmin, xmax, ymin, ymax)
        where xmin, etc. can be floats or quantities with units of
        length. The default is to set the extent to a box of size such
        that 99% of particles are contained within.
    units
        The units of the plot as a dictionary. The keys correspond to
        quantities such as 'position', 'density', 'velocity', and so on.
        The values are strings representing units, e.g. 'g/cm^3' for
        density. There is a special key 'projection' that corresponds
        to the length unit in the direction of projection for projected
        interpolation plots.
    ax
        A matplotlib Axes handle.
    ax_kwargs
        Keyword arguments to pass to matplotlib Axes.
    colorbar_kwargs
        Keyword arguments to pass to matplotlib Colorbar.
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
    colorbars, etc. All other keyword arguments are passed to the
    appropriate matplotlib function.

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

    Examples
    --------
    Show an image of the surface density in xy-plane.

    >>> plonk.image(snap=snap, quantity='density')

    Alternatively, access the function as a method on the Snap object.

    >>> snap.image(quantity='density')

    Set units for the plot.

    >>> units = {'position': 'au', 'density': 'g/cm^3', 'projection': 'cm'}
    >>> snap.image(quantity='density', units=units)

    Show a slice image of the density in xy-plane at z=0.

    >>> snap.image(quantity='density', interp='slice')
    """
    return _interpolation_plot(
        snap=snap,
        quantity=quantity,
        x=x,
        y=y,
        kind='image',
        interp=interp,
        slice_normal=slice_normal,
        slice_offset=slice_offset,
        extent=extent,
        units=units,
        ax=ax,
        ax_kwargs=ax_kwargs,
        colorbar_kwargs=colorbar_kwargs,
        **kwargs,
    )


def vector(
    snap: SnapLike,
    quantity: str,
    *,
    x: str = 'x',
    y: str = 'y',
    interp: str = 'projection',
    slice_normal: Tuple[float, float, float] = None,
    slice_offset: Union[Quantity, float] = None,
    extent: Quantity = None,
    units: Dict[str, str] = None,
    ax: Any = None,
    ax_kwargs={},
    **kwargs,
) -> Any:
    """Visualize vector SPH data as a vector plot.

    Visualize scalar smoothed particle hydrodynamics data by
    interpolation to a pixel grid of arrows.

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
    interp
        The interpolation type. Default is 'projection'.

        - 'projection' : 2d interpolation via projection to xy-plane
        - 'slice' : 3d interpolation via cross-section slice.
    slice_normal
        The normal vector to the plane in which to take the
        cross-section slice as a tuple (x, y, z). Default is
        (0, 0, 1).
    slice_offset
        The offset of the cross-section slice. Default is 0.0.
    extent
        The range in the x and y-coord as (xmin, xmax, ymin, ymax)
        where xmin, etc. can be floats or quantities with units of
        length. The default is to set the extent to a box of size such
        that 99% of particles are contained within.
    units
        The units of the plot as a dictionary. The keys correspond to
        quantities such as 'position', 'density', 'velocity', and so on.
        The values are strings representing units, e.g. 'g/cm^3' for
        density. There is a special key 'projection' that corresponds
        to the length unit in the direction of projection for projected
        interpolation plots.
    ax
        A matplotlib Axes handle.
    ax_kwargs
        Keyword arguments to pass to matplotlib Axes.
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
    passed to lower level functions as required.

    See below for additional parameters for interpolation, vector
    properties, etc. All other keyword arguments are passed to the
    appropriate matplotlib function.

    Other Parameters
    ----------------
    number_of_pixels : tuple
        The number of pixels to interpolate particle quantities
        to as a tuple (nx, ny). Default is (512, 512).
    density_weighted : bool
        Whether to density weight the interpolation or not.
        Default is False.
    number_of_arrows : tuple
        The number of arrows to display by sub-sampling the
        interpolated data. Default is (25, 25).
    normalize_vectors : bool
        Whether to normalize the arrows to all have the same
        length. Default is False.

    Examples
    --------
    Show a vector plot of velocity in xy-plane.

    >>> plonk.vector(snap=snap, quantity='velocity')

    Alternatively, access the function as a method on the Snap object.

    >>> snap.vector(quantity='velocity')

    Set units for the plot.

    >>> units = {'position': 'au', 'velocity': 'km/s', 'projection': 'km'}
    >>> snap.vector(quantity='velocity', units=units)

    Show a slice plot of the velocity in xy-plane at z=0.

    >>> snap.vector(quantity='density', interp='slice')
    """
    return _interpolation_plot(
        snap=snap,
        quantity=quantity,
        x=x,
        y=y,
        kind='quiver',
        interp=interp,
        slice_normal=slice_normal,
        slice_offset=slice_offset,
        extent=extent,
        units=units,
        ax=ax,
        ax_kwargs=ax_kwargs,
        **kwargs,
    )


def _interpolation_plot(
    snap: SnapLike,
    quantity: str,
    *,
    x: str = 'x',
    y: str = 'y',
    kind: str = None,
    interp: str = 'projection',
    slice_normal: Tuple[float, float, float] = None,
    slice_offset: Union[Quantity, float] = None,
    extent: Quantity = None,
    units: Dict[str, str] = None,
    ax: Any = None,
    ax_kwargs={},
    colorbar_kwargs={},
    **kwargs,
) -> Any:
    logger.debug(f'Visualizing "{quantity}" on snap: {snap.file_path.name}')
    _kwargs = copy(kwargs)

    if interp not in ('projection', 'slice'):
        raise ValueError('interp must be "projection" or "slice"')

    q: Quantity = snap[quantity]
    if kind is None:
        if q.ndim == 1:
            kind = 'image'
        elif q.ndim == 2:
            kind = 'quiver'
    if kind not in ('image', 'contour', 'quiver', 'streamplot'):
        raise ValueError('Cannot determine plot type')
    if kind in ('image', 'contour') and not q.ndim == 1:
        raise ValueError('image and contour plots must be of 1-dimensional quantities')
    if kind in ('quiver', 'streamplot') and quantity not in snap._vector_arrays:
        raise ValueError('quiver and stream plots must be of vector quantities')

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
    _data, _extent, _units = _interpolated_data(
        snap=snap,
        quantity=quantity,
        x=x,
        y=y,
        interp=interp,
        slice_normal=slice_normal,
        slice_offset=slice_offset,
        extent=extent,
        units=units,
        **interp_kwargs,
    )

    # Make the actual plot
    _interpolated_plot(
        interpolated_data=_data,
        extent=_extent,
        names={'quantity': quantity, 'x': x, 'y': y},
        kind=kind,
        interp=interp,
        units=_units,
        ax=ax,
        fig=fig,
        ax_kwargs=ax_kwargs,
        colorbar_kwargs=colorbar_kwargs,
        **_kwargs,
    )

    return ax


def _interpolated_data(
    snap, quantity, x, y, interp, slice_normal, slice_offset, extent, units, **kwargs,
):
    if units is None:
        _units = {
            'quantity': 1 * snap[quantity].units,
            'extent': 1 * snap['position'].units,
            'projection': 1 * snap['position'].units,
        }
    else:
        qunit = 1 * plonk_units(
            units.get(snap.base_array_name(quantity), str(snap[quantity].units))
        )
        eunit = 1 * plonk_units(units.get('position', str(snap['position'].units)))
        punit = 1 * plonk_units(units.get('projection', str(snap['position'].units)))
        _units = {'quantity': qunit, 'extent': eunit, 'projection': punit}

    if extent is None:
        _extent = get_extent_from_percentile(snap=snap, x=x, y=y)
        _extent = np.array([e.magnitude for e in _extent]) * _extent[0].units
    else:
        _extent = extent
    if not isinstance(_extent[0], Quantity):
        _extent = (
            _extent[0] * _units['extent'],
            _extent[1] * _units['extent'],
            _extent[2] * _units['extent'],
            _extent[3] * _units['extent'],
        )
    __extent = (
        (_extent[0] / snap.code_units['length']).to_base_units().magnitude,
        (_extent[1] / snap.code_units['length']).to_base_units().magnitude,
        (_extent[2] / snap.code_units['length']).to_base_units().magnitude,
        (_extent[3] / snap.code_units['length']).to_base_units().magnitude,
    )

    if interp == 'slice':
        if slice_normal is None:
            slice_normal = (0, 0, 1)
        if slice_offset is None:
            slice_offset = 0 * plonk_units('cm')
        if isinstance(slice_offset, Quantity):
            slice_offset = (
                (slice_offset / snap.code_units['length']).to_base_units().magnitude
            )
        else:
            logger.warning('slice_offset has no units, assuming code units')

    # Interpolate in code units
    interpolated_data = interpolate(
        snap=snap,
        quantity=quantity,
        x=x,
        y=y,
        interp=interp,
        slice_normal=slice_normal,
        slice_offset=slice_offset,
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


def _interpolated_plot(
    interpolated_data,
    extent,
    names,
    kind,
    interp,
    units,
    ax,
    fig,
    ax_kwargs,
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
    xname, yname = pretty_array_name(names["x"]), pretty_array_name(names["y"])
    ax.set_xlabel(f'{xname} [{eunit:~P}]')
    ax.set_ylabel(f'{yname} [{eunit:~P}]')

    ratio = (extent[1] - extent[0]) / (extent[3] - extent[2])
    if not max(ratio, 1 / ratio) > 10.0:
        ax.set_aspect('equal')

    ax.set(**ax_kwargs)

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
        elif interp == 'slice':
            qunit = units['quantity']
        if np.allclose(qunit.magnitude, 1.0):
            qunit = qunit.units
        qname = pretty_array_name(names["quantity"])
        cbar.set_label(f'{qname} [{qunit:~P}]')


def plot(
    snap: SnapLike,
    *,
    x: str = 'x',
    y: str = 'y',
    c: str = None,
    s: str = None,
    units: Dict[str, str] = None,
    xlim: Quantity = None,
    ylim: Quantity = None,
    ax: Any = None,
    ax_kwargs={},
    colorbar_kwargs={},
    **kwargs,
) -> Any:
    """Visualize SPH data as a particle plot.

    Visualize SPH data by plotting the particles, or a subset of
    the particles, possibly with marker colors and different sizes.

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
        The units of the plot as a dictionary. The keys correspond to
        quantities such as 'position', 'density', 'velocity', and so on.
        The values are strings representing units, e.g. 'g/cm^3' for
        density.
    xlim
        The range in the x-coord as (xmin, xmax) where xmin/xmax can be
        floats or quantities with units of length.
    ylim
        The range in the y-coord as (ymin, ymax) where ymin/ymax can be
        floats or quantities with units of length.
    ax
        A matplotlib Axes handle.
    ax_kwargs
        Keyword arguments to pass to matplotlib Axes.
    colorbar_kwargs
        Keyword arguments to pass to matplotlib Colorbar.
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

    >>> plonk.plot(snap=snap)

    Alternatively, access the function as a method on the Snap object.

    >>> snap.plot()

    Plot density against x.

    >>> snap.plot(x='x', y='density')

    Color particles by density.

    >>> snap.plot(x='x', y='y', c='density')

    Set units for the plot.

    >>> units = {'position': 'au', 'density': 'g/cm^3'}
    >>> snap.plot(x='x', y='y', c='density', units=units)
    """
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
        _x, _y, _c, _s, _units = _plot_data(
            snap=subsnap, x=x, y=y, c=c, s=s, units=units
        )
        _plot_plot(
            x=_x,
            y=_y,
            c=_c,
            s=_s,
            units=_units,
            xlim=xlim,
            ylim=ylim,
            names={'x': x, 'y': y, 'c': c, 's': s},
            fig=fig,
            ax=ax,
            ax_kwargs=ax_kwargs,
            colorbar_kwargs=colorbar_kwargs,
            **_kwargs,
        )

    return ax


def _plot_data(snap, x, y, c, s, units):
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
        xunit = units.get(snap.base_array_name(x), str(snap[x].units))
        yunit = units.get(snap.base_array_name(y), str(snap[y].units))
        _units = {
            'x': 1 * plonk_units(xunit),
            'y': 1 * plonk_units(yunit),
        }
        if c is not None:
            cunit = units.get(snap.base_array_name(c), str(snap[c].units))
            _units['c'] = 1 * plonk_units(cunit)
        if s is not None:
            sunit = units.get(snap.base_array_name(s), str(snap[s].units))
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
    except KeyError:
        # sink particles do not have smoothing length
        pass

    if _s is not None:
        _s = 100 * _s / _s.max()

    return _x, _y, _c, _s, _units


def _plot_plot(
    x, y, c, s, units, xlim, ylim, names, fig, ax, ax_kwargs, colorbar_kwargs, **kwargs,
):
    show_colorbar = kwargs.pop('show_colorbar', c is not None)

    if s is None and c is None:
        plots.plot(x=x, y=y, ax=ax, **kwargs)
    else:
        plot_object = plots.scatter(x=x, y=y, c=c, s=s, ax=ax, **kwargs)

    ratio = (x.max() - x.min()) / (y.max() - y.min())
    if not max(ratio, 1 / ratio) > 10.0:
        ax.set_aspect('equal')

    xunit, yunit = units['x'], units['y']
    if np.allclose(xunit.magnitude, 1.0):
        xunit = xunit.units
    if np.allclose(yunit.magnitude, 1.0):
        yunit = yunit.units
    xname, yname = pretty_array_name(names["x"]), pretty_array_name(names["y"])
    ax.set_xlabel(f'{xname} [{xunit:~P}]')
    ax.set_ylabel(f'{yname} [{yunit:~P}]')

    ax.set(**ax_kwargs)

    if xlim is not None:
        if not isinstance(xlim, Quantity):
            _xlim = xlim * xunit
        else:
            _xlim = xlim.to(xunit)
        ax.set_xlim(_xlim.magnitude)
    if ylim is not None:
        if not isinstance(ylim, Quantity):
            _ylim = ylim * yunit
        else:
            _ylim = ylim.to(yunit)
        ax.set_ylim(_ylim.magnitude)

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
        cname = pretty_array_name(names["c"])
        cbar.set_label(f'{cname} [{cunit:~P}]')


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
    quantity_unit = snap.array_code_unit(quantity)
    if interp == 'projection':
        proj_unit = units['quantity'] * units['projection']
        data = (interpolated_data * quantity_unit * snap.code_units['length']).to(
            proj_unit.units
        ).magnitude / proj_unit.magnitude
    elif interp == 'slice':
        data = (interpolated_data * quantity_unit).to(
            units['quantity'].units
        ).magnitude / units['quantity'].magnitude

    new_extent = tuple(
        (extent * snap.code_units['length']).to(units['extent'].units).magnitude
        / units['extent'].magnitude
    )

    return data, new_extent


def _convert_units_for_cmap(vm, name, units, interp):
    if interp == 'projection':
        quantity_unit = units['quantity'] * units['projection']
    elif interp == 'slice':
        quantity_unit = units['quantity']
    if vm is not None:
        if isinstance(vm, Quantity):
            if vm.dimensionality != quantity_unit.dimensionality:
                raise ValueError(f'{name} has incorrect units')
            vm = (vm / quantity_unit).to_base_units().magnitude
        else:
            logger.warning(f'{name} has no units, assuming same units as quantity')

    return vm
