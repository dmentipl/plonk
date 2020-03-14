"""The Visualization class.

This class contains methods for visualizing smoothed particle
hydrodynamics simulation data.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional, Union

import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable
from numpy import ndarray

from .. import Quantity
from .. import units as plonk_units
from ..snap import SnapLike
from ..snap.snap import get_array_from_input
from ..utils import get_extent_from_percentile, is_documented_by
from . import plots
from .interpolation import Extent, interpolate

_kind_to_object = {
    'render': 'image',
    'contour': 'contour',
    'quiver': 'quiver',
    'stream': 'streamplot',
}

_kind_to_function = {
    'render': plots.render_plot,
    'contour': plots.contour_plot,
    'quiver': plots.quiver_plot,
    'stream': plots.stream_plot,
}


class Visualization:
    """Visualize scalar and vector smoothed particle hydrodynamics data.

    Visualize SPH data as a particle plot, a rendered image, a contour
    plot, a vector plot, or a stream plot.

    Parameters
    ----------
    snap
        The associated Snap (or SubSnap) object to visualize.

    Attributes
    ----------
    snap
        The associated Snap (or SubSnap) object to visualize.
    fig
        The matplotlib Figure object of the plot.
    ax
        The matplotlib Axes object of the plot.
    units
        The units of the plot: 'quantity', 'extent', 'projection'. The
        values are pint Unit objects.
    extent
        A tuple (xmin, xmax, ymin, ymax) of the extent of the plot.
    objects
        A dictionary containing the matplotlib plot objects:

        - 'lines' : list of matplotlib Line2D objects for particle plots
        - 'image' : matplotlib AxesImage object for rendered plots
        - 'colorbar' : matplotlib Colorbar object for rendered plots
        - 'contour' : matplotlib QuadContourSet object for contour plots
        - 'quiver' : matplotlib Quiver object for quiver (arrow) plots
        - 'streamplot' : matplotlib StreamplotSet for stream plots
    data
        A dictionary containing the data interpolated to a pixel grid:
        'render', 'contour', 'quiver', 'stream'.
    """

    def __init__(self, snap: SnapLike):
        self.snap = snap
        self.fig: Any = None
        self.ax: Any = None
        self.units: Dict[str, Any] = {
            'quantity': None,
            'extent': None,
            'projection': None,
        }
        self.extent: Extent = (-1, -1, -1, -1)
        self.objects: Dict[str, Any] = {
            'lines': None,
            'image': None,
            'colorbar': None,
            'contour': None,
            'quiver': None,
            'streamplot': None,
        }
        self.data: Dict[str, ndarray] = {
            'render': None,
            'contour': None,
            'quiver': None,
            'stream': None,
        }

    def plot(
        self,
        *,
        quantity: Optional[Union[str, ndarray]] = None,
        x: Union[str, ndarray] = 'x',
        y: Union[str, ndarray] = 'y',
        z: Union[str, ndarray] = 'z',
        kind: Optional[str] = None,
        interp: str = 'projection',
        z_slice: float = 0.0,
        extent: Extent = (-1, -1, -1, -1),
        units=None,
        ax: Optional[Any] = None,
        **kwargs,
    ) -> Visualization:
        """Visualize smoothed particle hydrodynamics data.

        Visualize SPH data by showing the particle positions, a rendered
        image or a contour plot for scalar data, a quiver (arrow) plot
        or stream plot for vector data.

        Parameters
        ----------
        quantity
            The quantity to visualize. Can be a string to pass to Snap,
            or a 1d array (N,) of scalar data, or a 2d array (N, 3) of
            vector data. If quantity is 2d, only the first two
            components are visualized, i.e. quantity[:, 0] and
            quantity[:, 1]. Default is None.
        x
            The x-coordinate for the visualization. Can be a string to
            pass to Snap, or a 1d array (N,). Default is 'x'.
        y
            The y-coordinate for the visualization. Can be a string to
            pass to Snap, or a 1d array (N,). Default is 'y'.
        z
            The z-coordinate for the visualization. Can be a string to
            pass to Snap, or a 1d array (N,). This is only required for
            cross-section plots. Default is 'z'.
        kind
            The type of plot.

            - 'particle' : particle plot (default if quantity is None)
            - 'render' : rendered image (default for scalar quantities)
            - 'contour' : contour plot (scalar quantity)
            - 'quiver' : quiver (arrow) plot (default for vector
              quantities)
            - 'stream' : stream plot (vector quantity)
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
            The default is to set the extent to a box of size set by an
            inscribed sphere containing 99% of particles.
        units
            The units of the plot as a dictionary with keys 'quantity',
            'extent', 'projection'. The values are Pint Unit objects.
        ax
            A matplotlib Axes handle.
        **kwargs
            Additional keyword arguments to pass to interpolation and
            matplotlib functions.

        Notes
        -----
        Additional parameters passed as keyword arguments will be
        passed to lower level functions as required. E.g. Plonk uses
        matplotlib's imshow for a render plot, so additional arguments
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
        fmt : str
            This is the matplotlib ax.plot method positional
            argument format string. Default is 'k.'.
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
        Plot the particles.

        >>> viz = plonk.visualize.plot(snap=snap)

        Render the surface density in xy-plane.

        >>> viz = plonk.visualize.plot(
        ...     snap=snap,
        ...     quantity='density',
        ... )

        Get the interpolation to grid directly (without plotting).

        >>> grid_data = plonk.visualize.interpolate(
        ...     snap=snap,
        ...     quantity='density',
        ...     interp='projection',
        ...     extent=(-100, 100, -100, 100),
        ... )

        Make an animation of multiple snaps.

        >>> plonk.visualize.animation(
        ...     snaps=snaps,
        ...     quantity='density',
        ...     extent=(-100, 100, -100, 100),
        ...     vmin=0.0,
        ...     vmax=1.0,
        ...     filename='animation.mp4',
        ... )

        Set units for the plot.

        >>> units = {
        ...     'quantity': plonk.units('g / cm ** 3'),
        ...     'extent': plonk.units('au'),
        ...     'projection': plonk.units('cm'),
        ... }

        >>> viz = plonk.visualize.plot(
        ...     snap=snap,
        ...     quantity='density',
        ...     units=units,
        ... )

        You can also use the utility function to generate the units
        dictionary.

        >>> units = plonk.visualize.str_to_units('g/cm^3', 'au', 'cm')
        """
        if self.ax is None:
            if ax is None:
                self.fig, self.ax = plt.subplots()
            else:
                self.fig = ax.get_figure()
                self.ax = ax
        else:
            if ax is not None:
                raise ValueError('Trying to change existing Axes attribute')

        if extent == (-1, -1, -1, -1):
            extent = get_extent_from_percentile(self.snap)

        quantity_str: Optional[str] = quantity if isinstance(quantity, str) else None
        quantity, x, y, z, kind = _check_input(
            snap=self.snap, quantity=quantity, x=x, y=y, z=z, kind=kind
        )

        if quantity is not None:
            interpolation_kwargs = ('number_of_pixels', 'density_weighted')
            _kwargs = {
                key: val for key, val in kwargs.items() if key in interpolation_kwargs
            }
            for key in _kwargs:
                kwargs.pop(key)
            interpolated_data = interpolate(
                snap=self.snap,
                quantity=quantity,
                x=x,
                y=y,
                z=z,
                interp=interp,
                z_slice=z_slice,
                extent=extent,
                **_kwargs,
            )
            if units is not None:
                if quantity_str is None:
                    raise ValueError(
                        'Cannot set units when passing in arrays. '
                        'Instead, use strings to access\n'
                        'quantities on Snap. E.g. plot(..., quantity="density", ...).'
                    )
                interpolated_data, extent = _convert_units(
                    snap=self.snap,
                    quantity_str=quantity_str,
                    interpolated_data=interpolated_data,
                    extent=extent,
                    units=units,
                    interp=interp,
                )
                self.units.update(units)

        self.extent = extent

        show_colorbar = kwargs.pop('show_colorbar', True if kind == 'render' else False)

        if kind == 'particle':
            self.objects['lines'] = plots.particle_plot(
                snap=self.snap, x=x, y=y, extent=extent, ax=self.ax, **kwargs,
            )

        elif kind in ('render', 'contour', 'quiver', 'stream'):
            self.data[kind] = interpolated_data
            self.objects[_kind_to_object[kind]] = _kind_to_function[kind](
                interpolated_data=interpolated_data,
                extent=extent,
                ax=self.ax,
                **kwargs,
            )

        else:
            raise ValueError('Cannot determine plot type')

        if show_colorbar:
            divider = make_axes_locatable(self.ax)
            cax = divider.append_axes("right", size="5%", pad=0.05)
            self.objects['colorbar'] = self.fig.colorbar(self.objects['image'], cax)

        self.ax.set_xlim(*extent[:2])
        self.ax.set_ylim(*extent[2:])
        self.ax.set_aspect('equal')

        return self

    def __repr__(self):
        """Dunder repr method."""
        return '<plonk.Visualization>'


class MultiVisualization:
    """Visualize multiple snaps.

    Attributes
    ----------
    snaps
        A list of snaps.
    options
        A dictionary of arguments passed to Visualization plot method.
    visualization
        The Visualization object.
    ax
        The matplotlib Axes object.

    Parameters
    ----------
    snaps
        A list of Snap objects.
    **kwargs
        Keyword arguments to pass to Visualization plot method.
    """

    def __init__(self, snaps, **kwargs):
        self.snaps = snaps
        self.options = kwargs
        self.visualization = Visualization(snap=snaps[0]).plot(**kwargs)
        self.ax = self.visualization.ax

        self._len = -1
        self._where = 0

    def _fn(self, idx):
        self.ax.clear()
        cbar = self.visualization.objects['colorbar']
        if cbar is not None:
            cbar.remove()
        viz = Visualization(snap=self.snaps[idx]).plot(ax=self.ax, **self.options)
        self.visualization = viz
        return viz

    def next(self, number: int = 1):
        """Visualize next snap."""
        idx = self._where + number
        if idx < len(self):
            self.visualization = self._fn(idx)
            self._where += number
        else:
            print('Too far forward. Going to last snap.')
            self.visualization = self._fn(len(self) - 1)
            self._where = len(self) - 1

    def prev(self, number: int = 1):
        """Visualize previous snap."""
        idx = self._where - number
        if idx > 0:
            self.visualization = self._fn(idx)
            self._where -= number
        else:
            print('Too far back. Going to first snap.')
            self.visualization = self._fn(0)
            self._where = 0

    def goto(self, idx):
        """Visualize particular snap by index."""
        if -len(self) < idx < len(self) - 1:
            self.visualization = self._fn(idx)
            self._where = np.mod(idx, len(self))
        else:
            raise ValueError('out of range')

    @property
    def index(self):
        """Current snap index."""
        return self._where

    def __len__(self):
        """Length as number of snaps."""
        if self._len == -1:
            self._len = len(self.snaps)
        return self._len


def _check_input(*, snap, quantity, x, y, z, kind):

    if quantity is not None:
        quantity = get_array_from_input(snap, quantity)

    x = get_array_from_input(snap, x, 'x')
    y = get_array_from_input(snap, y, 'y')
    z = get_array_from_input(snap, z, 'z')

    if quantity is not None:
        if quantity.ndim > 2:
            raise ValueError('Cannot interpret quantity')
        if kind in ('render', 'contour') and quantity.ndim != 1:
            raise ValueError('quantity is wrong shape for render or contour')
        if kind in ('quiver', 'stream') and quantity.ndim != 2:
            raise ValueError('quantity is wrong shape for quiver or streamplot')
        if kind is None:
            if quantity.ndim == 1:
                kind = 'render'
            elif quantity.ndim == 2:
                kind = 'quiver'
    else:
        if kind is None:
            kind = 'particle'
        else:
            raise ValueError(f'No quantity: can only do particle plot')

    return quantity, x, y, z, kind


def _convert_units(
    *,
    snap: SnapLike,
    quantity_str: str,
    interpolated_data: ndarray,
    extent: Extent,
    units: Dict[str, Any],
    interp: str,
):
    required_keys = {'extent', 'projection', 'quantity'}
    if not set(units) == required_keys:
        raise ValueError(f'units dictionary requires: {required_keys}')
    _quantity_unit = snap.get_array_unit(quantity_str)
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


@is_documented_by(Visualization.plot)
def plot(
    *,
    snap: SnapLike,
    quantity: Optional[Union[str, ndarray]] = None,
    x: Union[str, ndarray] = 'x',
    y: Union[str, ndarray] = 'y',
    z: Union[str, ndarray] = 'z',
    kind: Optional[str] = None,
    interp: str = 'projection',
    z_slice: float = 0.0,
    extent: Extent = (-1, -1, -1, -1),
    ax: Optional[Any] = None,
    **kwargs,
) -> Visualization:
    viz = Visualization(snap)
    viz.plot(
        quantity=quantity,
        x=x,
        y=y,
        z=z,
        kind=kind,
        interp=interp,
        z_slice=z_slice,
        extent=extent,
        ax=ax,
        **kwargs,
    )
    return viz


def plot_snaps(snaps: List[SnapLike], **kwargs) -> MultiVisualization:
    """Visualize multiple snaps.

    Parameters
    ----------
    snaps
        A list of Snap objects.
    **kwargs
        Keyword arguments to pass to Visualization plot method.

    Returns
    -------
    MultiVisualization

    Examples
    --------
    Initialize object passing in plotting parameters.

    >>> vi = plot_snaps(
    ...     snaps=sim.snaps,
    ...     quantity='density',
    ... )

    Go forwards and backwards through snaps.

    >>> vi.next()
    >>> vi.prev()

    Go to a particular snap, or skip ahead.

    >>> vi.goto(10)
    >>> vi.next(5)
    """
    return MultiVisualization(snaps, **kwargs)


def str_to_units(quantity, extent, projection):
    """Convert string to plonk units.

    Parameters
    ----------
    quantity
        The units string for the quantity.
    extent
        The units string for the plot extent.
    projection
        The units string for projection interpolation.
    """
    if isinstance(quantity, str):
        quantity = plonk_units(quantity)
    elif isinstance(quantity, Quantity):
        pass
    else:
        raise ValueError(f'Cannot determine quantity unit')
    if isinstance(extent, str):
        extent = plonk_units(extent)
    elif isinstance(extent, Quantity):
        pass
    else:
        raise ValueError(f'Cannot determine extent unit')
    if isinstance(projection, str):
        projection = plonk_units(projection)
    elif isinstance(projection, Quantity):
        pass
    else:
        raise ValueError(f'Cannot determine projection unit')

    if extent.dimensionality != plonk_units('cm').dimensionality:
        raise ValueError('extent has incorrect dimensions')
    if projection.dimensionality != plonk_units('cm').dimensionality:
        raise ValueError('projection has incorrect dimensions')

    return {'quantity': quantity, 'extent': extent, 'projection': projection}
