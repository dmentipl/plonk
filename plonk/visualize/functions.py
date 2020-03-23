"""Functions for visualization."""

from typing import Any, Dict, List, Optional, Union

from numpy import ndarray

from .. import Quantity
from .. import units as plonk_units
from ..snap import SnapLike
from ..utils import is_documented_by
from .interpolation import Extent
from .multi import MultiVisualization
from .visualization import Visualization


@is_documented_by(Visualization.plot)
def plot(
    *,
    snap: SnapLike,
    quantity: Union[str, ndarray],
    x: Union[str, ndarray] = 'x',
    y: Union[str, ndarray] = 'y',
    z: Union[str, ndarray] = 'z',
    kind: Optional[str] = None,
    interp: str = 'projection',
    z_slice: float = 0.0,
    extent: Extent = (-1, -1, -1, -1),
    units: Dict[str, Any] = None,
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
        units=units,
        ax=ax,
        **kwargs,
    )
    return viz


@is_documented_by(Visualization.particle_plot)
def particle_plot(
    *,
    snap: SnapLike,
    x: Union[str, ndarray] = 'x',
    y: Union[str, ndarray] = 'y',
    color: Optional[Union[str, ndarray]] = None,
    size: Optional[Union[str, ndarray]] = None,
    units: Dict[str, Any] = None,
    ax: Optional[Any] = None,
    **kwargs,
) -> Visualization:
    viz = Visualization(snap)
    viz.particle_plot(
        x=x, y=y, color=color, size=size, units=units, ax=ax, **kwargs,
    )
    return viz


def plot_snaps(
    snaps: List[SnapLike], quantity: Union[str, ndarray], **kwargs
) -> MultiVisualization:
    """Visualize multiple snaps.

    Parameters
    ----------
    snaps
        A list of Snap objects.
    quantity
        The quantity to visualize.
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
    return MultiVisualization(snaps=snaps, quantity=quantity, **kwargs)


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
