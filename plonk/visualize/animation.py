"""Animations of visualizations."""

import pathlib
from copy import copy
from pathlib import Path
from typing import Any, Dict, List, Tuple, Union

import matplotlib.pyplot as plt
from matplotlib import animation as _animation
from numpy import ndarray

from ..snap.snap import Snap, SubSnap
from .plot import render

SnapLike = Union[Snap, SubSnap]

# TODO:
# - fix color scale
# - remove axis labels
# - black background


def animation(
    snaps: List[SnapLike],
    quantity: Union[str, ndarray],
    extent: Tuple[float, float, float, float],
    filename: Union[str, Path],
    scalar_options: Dict[Any, Any] = None,
    interpolation_options: Dict[Any, Any] = None,
    animation_options: Dict[Any, Any] = None,
):
    """Generate an animation.

    Parameters
    ----------
    snaps
        A list of Snap objects to animate.
    quantity
        The quantity on the Snap objects to render.
    extent
        The extent of the image like (xmin, xmax, ymin, ymax).
    filename
        The file name to save the animation to.
    scalar_options : optional
        Scalar render options to pass to the render function.
    interpolation_options : optional
        Interpolation options to pass to the render function.
    animation_options : optional
        Options to pass to matplotlib Animation.save.
    """
    filepath = pathlib.Path(filename)
    if filepath.suffix != '.mp4':
        raise ValueError('filename should end in ".mp4"')

    if scalar_options is None:
        scalar_options = {}
    _scalar_options = copy(scalar_options)

    if animation_options is None:
        animation_options = {}

    fig, ax = plt.subplots()
    viz = render(
        snaps[0],
        quantity,
        extent=extent,
        axis=ax,
        scalar_options=_scalar_options,
        interpolation_options=interpolation_options,
    )
    data = viz.data['scalar']
    im = ax.imshow(data)

    # If plotting colorbar, must fix the render range
    plot_colorbar = _scalar_options.get('plot_colorbar')
    if plot_colorbar is None or plot_colorbar is True:
        _scalar_options['render_range'] = (data.min(), data.max())

    def animate(idx):
        print(f'Rendering image: {idx}')
        viz = render(
            snaps[idx],
            quantity,
            extent=extent,
            axis=ax,
            scalar_options=_scalar_options,
            interpolation_options=interpolation_options,
        )
        data = viz.data['scalar']
        im.set_data(data)
        return [im]

    anim = _animation.FuncAnimation(fig, animate, frames=len(snaps))
    anim.save(
        filepath, extra_args=['-vcodec', 'libx264'], **animation_options,
    )
    plt.close()
