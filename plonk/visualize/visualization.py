"""The Visualization class.

This class contains methods for visualizing smoothed particle
hydrodynamics simulation data.
"""

from __future__ import annotations

from typing import Any, Dict, Optional, Tuple, Union

import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from numpy import ndarray

from ..snap.snap import SnapLike, get_array_from_input
from ..utils import is_documented_by
from . import _plots
from .interpolation import interpolate

_kind_to_object = {
    'render': 'image',
    'contour': 'contour',
    'quiver': 'quiver',
    'stream': 'streamplot',
}

_kind_to_function = {
    'render': _plots._render_plot,
    'contour': _plots._contour_plot,
    'quiver': _plots._quiver_plot,
    'stream': _plots._stream_plot,
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
    axis
        The matplotlib Axis object of the plot.
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
        self.axis: Any = None
        self.extent: Tuple[float, float, float, float] = (-1, -1, -1, -1)
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
        extent: Tuple[float, float, float, float],
        axis: Optional[Any] = None,
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
            The interpolation type.
            - 'projection' : 2d interpolation via projection to xy-plane
            - 'cross_section' : 3d interpolation via cross-section in
              z-direction
            Default is 'projection'.
        z_slice
            The z-coordinate value of the cross-section slice. Default
            is 0.0.
        extent
            The range in the x and y-coord as (xmin, xmax, ymin, ymax).
        axis
            A matplotlib axis handle.
        **kwargs
            Additional key word arguments to pass to interpolation and
            matplotlib functions.

        Notes
        -----
        Additional parameters passed as key word arguments will be
        passed to lower level functions as required. E.g. Plonk uses
        matplotlib's imshow for a render plot, so additional arguments
        to imshow can be passed this way.

        Other Parameters
        ----------------
        This is a list of parameters that are passed in as key word
        arguments to other functions inside this method.

        Parameters for interpolation:

            number_of_pixels : tuple
                The number of pixels to interpolate particle quantities
                to as a tuple (nx, ny). Default is (512, 512).
            density_weighted : bool
                Whether to density weight the interpolation or not.
                Default is False.

        Parameters for particle plots are passed to axis.plot except:

            fmt : str
                This is the matplotlib axis.plot method positional
                argument format string. Default is 'k.'.

        Parameters for kind='render' plots are passed to axis.imshow
        except:

            show_colorbar : bool
                Whether or not to display a colorbar. Default is True.

        Parameters for kind='contour' plots are passed to axis.contour.

        Parameters for kind='quiver' plots are passed to axis.quiver
        except:

            number_of_arrows : tuple
                The number of arrows to display by sub-sampling the
                interpolated data. Default is (25, 25).
            normalize_vectors : bool
                Whether to normalize the arrows to all have the same
                length. Default is False.

        Parameters for kind='stream' plots are passed to
        axis.streamplot.
        """
        if self.axis is None:
            if axis is None:
                self.fig, self.axis = plt.subplots()
            else:
                self.fig = axis.get_figure()
                self.axis = axis
        else:
            if axis is not None:
                raise ValueError('Trying to change existing axis attribute')

        self.extent = extent

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

        show_colorbar = kwargs.pop('show_colorbar', True if kind == 'render' else False)

        if kind == 'particle':
            self.objects['lines'] = _plots._particle_plot(
                snap=self.snap, x=x, y=y, extent=extent, axis=self.axis, **kwargs,
            )

        elif kind in ('render', 'contour', 'quiver', 'stream'):
            self.data[kind] = interpolated_data
            self.objects[_kind_to_object[kind]] = _kind_to_function[kind](
                interpolated_data=interpolated_data,
                extent=extent,
                axis=self.axis,
                **kwargs,
            )

        else:
            raise ValueError('Cannot determine plot type')

        if show_colorbar:
            divider = make_axes_locatable(self.axis)
            cax = divider.append_axes("right", size="5%", pad=0.05)
            self.objects['colorbar'] = self.fig.colorbar(self.objects['image'], cax)

        self.axis.set_xlim(*extent[:2])
        self.axis.set_ylim(*extent[2:])
        self.axis.set_aspect('equal')

        return self

    def __repr__(self):
        """Dunder repr method."""
        return '<plonk.Visualization>'


def _check_input(*, snap, quantity, x, y, z, kind):

    try:
        quantity = get_array_from_input(snap, quantity)
    except ValueError:
        quantity = None

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
    extent: Tuple[float, float, float, float],
    axis: Optional[Any] = None,
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
        axis=axis,
        **kwargs,
    )
    return viz
