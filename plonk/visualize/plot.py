"""Plot helper functions."""

from typing import Any, Optional, Tuple, Union

from numpy import ndarray

from ..snap.snap import SnapLike
from ..utils import is_documented_by
from .interpolation import interpolate as _interpolate
from .visualization import Visualization


@is_documented_by(Visualization.plot)
def plot(
    *,
    snap: SnapLike,
    data: Optional[Union[str, ndarray]] = None,
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
        data=data,
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
