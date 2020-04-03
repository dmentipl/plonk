"""The Visualization class.

This class contains methods for visualizing smoothed particle
hydrodynamics simulation data.
"""

from __future__ import annotations

from copy import copy
from typing import Any, Dict, Optional

import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from numpy import ndarray

from .. import Quantity
from ..snap import SnapLike
from ..utils import get_extent_from_percentile
from . import plots
from .interpolation import Extent, interpolate

_kind_to_function = {
    'image': plots.imshow,
    'contour': plots.contour,
    'quiver': plots.quiver,
    'streamplot': plots.streamplot,
}


class Visualization:
    """Visualize scalar and vector smoothed particle hydrodynamics data.

    Visualize SPH data as a scatter plot, an image, a contour plot, a
    vector plot, or a stream plot.

    Parameters
    ----------
    snap
        The associated Snap (or SubSnap) object to visualize.

    Attributes
    ----------
    snap
        The associated Snap (or SubSnap) object to visualize.
    ax
        The matplotlib Axes object of the plot.
    """

    def __init__(self, snap: SnapLike):
        self.snap = snap
        self.ax: Any = None

    def plot(
        self,
        *,
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
        quantity
            The quantity to visualize. Must be a string to pass to Snap,
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
        _kwargs = copy(kwargs)

        if self.ax is None:
            if ax is None:
                fig, self.ax = plt.subplots()
            else:
                fig = ax.figure
                self.ax = ax
        else:
            if ax is not None:
                raise ValueError('Trying to change existing Axes attribute')

        if extent == (-1, -1, -1, -1):
            extent = get_extent_from_percentile(snap=self.snap, x=x, y=y)
        if isinstance(extent[0], Quantity):
            extent = (
                (extent[0] / self.snap.units['length']).to_base_units().magnitude,
                (extent[1] / self.snap.units['length']).to_base_units().magnitude,
                (extent[2] / self.snap.units['length']).to_base_units().magnitude,
                (extent[3] / self.snap.units['length']).to_base_units().magnitude,
            )
        if isinstance(z_slice, Quantity):
            z_slice = (z_slice / self.snap.units['length']).to_base_units().magnitude

        interpolation_kwargs = ('number_of_pixels', 'density_weighted')
        __kwargs = {
            key: val for key, val in _kwargs.items() if key in interpolation_kwargs
        }
        for key in __kwargs:
            _kwargs.pop(key)
        interpolated_data = interpolate(
            snap=self.snap,
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
                snap=self.snap,
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

        show_colorbar = _kwargs.pop('show_colorbar', True if kind == 'image' else False)

        if kind in ('image', 'contour', 'quiver', 'streamplot'):
            plot_object = _kind_to_function[kind](
                interpolated_data=interpolated_data,
                extent=extent,
                ax=self.ax,
                **_kwargs,
            )

        else:
            raise ValueError('Cannot determine plot type')

        if show_colorbar:
            divider = make_axes_locatable(self.ax)
            _kwargs = copy(colorbar_kwargs)
            position = _kwargs.pop('position', 'right')
            size = _kwargs.pop('size', '5%')
            pad = _kwargs.pop('pad', '2%')
            if position in ('top', 'bottom'):
                _kwargs.update({'orientation': 'horizontal'})
            cax = divider.append_axes(position=position, size=size, pad=pad)
            fig.colorbar(plot_object, cax, **_kwargs)

        self.ax.set_xlim(*extent[:2])
        self.ax.set_ylim(*extent[2:])

        ratio = (extent[1] - extent[0]) / (extent[3] - extent[2])
        if not max(ratio, 1 / ratio) > 10.0:
            self.ax.set_aspect('equal')

        return self.ax

    def particle_plot(
        self,
        *,
        x: str = 'x',
        y: str = 'y',
        color: Optional[str] = None,
        size: Optional[str] = None,
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
        x
            The x-coordinate for the visualization. Must be a string to
            pass to Snap. Default is 'x'.
        y
            The y-coordinate for the visualization. Must be a string to
            pass to Snap. Default is 'y'.
        color
            The quantity to color the particles. Must be a string to
            pass to Snap.
        size
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
        _kwargs = copy(kwargs)

        if self.ax is None:
            if ax is None:
                fig, self.ax = plt.subplots()
            else:
                fig = ax.figure
                self.ax = ax
        else:
            if ax is not None:
                raise ValueError('Trying to change existing Axes attribute')

        _x: ndarray = self.snap[x]
        _y: ndarray = self.snap[y]
        _color: ndarray = self.snap[color] if color is not None else None
        _size: ndarray = self.snap[size] if size is not None else None

        if self.snap._physical_units:
            if xunit is not None:
                _x = _x.to(xunit).magnitude
            else:
                _x = _x.magnitude
            if yunit is not None:
                _y = _y.to(yunit).magnitude
            else:
                _y = _y.magnitude
            if cunit is not None:
                _color = _color.to(cunit).magnitude
            else:
                _color = _color.magnitude

        h: ndarray = self.snap['smoothing_length']
        mask = h > 0

        _x = _x[mask]
        _y = _y[mask]
        if _color is not None:
            _color = _color[mask]
        if _size is not None:
            _size = _size[mask]
            _size = 100 * _size / _size.max()
            if self.snap._physical_units:
                _size = _size.magnitude

        show_colorbar = _kwargs.pop(
            'show_colorbar', True if _color is not None else False
        )

        if _size is None and _color is None:
            plot_object = plots.plot(x=_x, y=_y, ax=self.ax, **_kwargs)

        else:
            plot_object = plots.scatter(
                x=_x, y=_y, color=_color, size=_size, ax=self.ax, **_kwargs,
            )
            if show_colorbar:
                divider = make_axes_locatable(self.ax)
                _kwargs = copy(colorbar_kwargs)
                position = _kwargs.pop('position', 'right')
                size = _kwargs.pop('size', '5%')
                pad = _kwargs.pop('pad', '2%')
                if position in ('top', 'bottom'):
                    _kwargs.update({'orientation': 'horizontal'})
                cax = divider.append_axes(position=position, size=size, pad=pad)
                fig.colorbar(plot_object, cax, **_kwargs)

        if xscale is not None:
            self.ax.set_xscale(xscale)
        if yscale is not None:
            self.ax.set_yscale(yscale)

        ratio = (_x.max() - _x.min()) / (_y.max() - _y.min())
        if not max(ratio, 1 / ratio) > 10.0:
            self.ax.set_aspect('equal')

        return self.ax

    def __repr__(self):
        """Dunder repr method."""
        return self.__str__()

    def __str__(self):
        """Dunder str method."""
        return f'<plonk.Visualization "{self.snap.file_path.name}">'


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
