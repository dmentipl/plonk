"""Animations of visualizations."""

import pathlib
from pathlib import Path
from typing import List, Optional, Tuple, Union

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
    x: Union[str, ndarray] = 'x',
    y: Union[str, ndarray] = 'y',
    z: Union[str, ndarray] = 'z',
    kind: Optional[str] = None,
    interp: str = 'projection',
    z_slice: float = 0.0,
    extent: Tuple[float, float, float, float],
    animation_kwargs,
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
        - 'particle' : particle plot (default if data is None)
        - 'render' : rendered image (default for scalar data)
        - 'contour' : contour plot (scalar data)
        - 'quiver' : quiver (arrow) plot (default for vector data)
        - 'stream' : stream plot (vector data)
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
        The extent of the image like (xmin, xmax, ymin, ymax).
    animation_kwargs : optional
        Key word arguments to pass to matplotlib Animation.save.
    """
    filepath = pathlib.Path(filename)
    if filepath.suffix != '.mp4':
        raise ValueError('filename should end in ".mp4"')

    fig, ax = plt.subplots()
    viz = plot(
        snap=snaps[0],
        data=data,
        x=x,
        y=y,
        z=z,
        kind=kind,
        interp=interp,
        z_slice=z_slice,
        extent=extent,
        **kwargs,
    )
    scalar_data = viz.data['scalar']
    im = ax.imshow(scalar_data)

    # If plotting colorbar, must fix the render range
    vmin, vmax = None, None
    show_colorbar = kwargs.pop('show_colorbar', True)
    if show_colorbar is None or show_colorbar is True:
        vmin, vmax = (scalar_data.min(), scalar_data.max())

    def animate(idx):
        print(f'Rendering image: {idx}')
        viz = plot(
            snap=snaps[idx],
            data=data,
            x=x,
            y=y,
            z=z,
            kind=kind,
            interp=interp,
            z_slice=z_slice,
            extent=extent,
            vmin=vmin,
            vmax=vmax,
            **kwargs,
        )
        scalar_data = viz.data['scalar']
        im.set_data(scalar_data)
        return [im]

    anim = _animation.FuncAnimation(fig, animate, frames=len(snaps))
    anim.save(
        filepath, extra_args=['-vcodec', 'libx264'], **animation_kwargs,
    )
    plt.close()
