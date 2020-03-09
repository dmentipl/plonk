"""Animations of visualizations."""

import pathlib
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

import matplotlib.pyplot as plt
from matplotlib import animation as _animation
from numpy import ndarray

from ..snap.snap import SnapLike
from .visualization import plot

# TODO:
# - fix color scale
# - remove axis labels
# - black background


def animation(
    *,
    filename: Union[str, Path],
    snaps: List[SnapLike],
    data: Optional[Union[str, ndarray]] = None,
    extent: Tuple[float, float, float, float],
    animation_kwargs: Dict[str, Any] = {},
    **kwargs,
):
    """Generate an animation.

    Parameters
    ----------
    filename
        The file name to save the animation to.
    snaps
        A list of Snap objects to animate.
    data
        The data to visualize. Can be a string to pass to Snap, or
        a 1d array (N,) of scalar data, or a 2d array (N, 3) of
        vector data. Default is None.
    extent
        The range in the x and y-coord as (xmin, xmax, ymin, ymax).
    animation_kwargs : optional
        Key word arguments to pass to matplotlib Animation.save.
    **kwargs
        Arguments to pass to visualize.plot.
    """
    filepath = pathlib.Path(filename)
    if filepath.suffix != '.mp4':
        raise ValueError('filename should end in ".mp4"')

    fig, ax = plt.subplots()
    viz = plot(snap=snaps[0], data=data, extent=extent, **kwargs)
    interp_data = _get_data(viz.data)
    im = ax.imshow(interp_data)

    # If plotting colorbar, must fix the render range
    vmin, vmax = None, None
    show_colorbar = kwargs.pop('show_colorbar', True)
    if show_colorbar is None or show_colorbar is True:
        vmin, vmax = (interp_data.min(), interp_data.max())

    def animate(idx):
        print(f'Rendering image: {idx}')
        viz = plot(
            snap=snaps[idx], data=data, extent=extent, vmin=vmin, vmax=vmax, **kwargs
        )
        interp_data = _get_data(viz.data)
        im.set_data(interp_data)
        return [im]

    anim = _animation.FuncAnimation(fig, animate, frames=len(snaps))
    anim.save(
        filepath, extra_args=['-vcodec', 'libx264'], **animation_kwargs,
    )
    plt.close()


def _get_data(data):
    for k, v in data.items():
        if v is not None:
            return v
