"""Animations of visualizations."""

import pathlib
from pathlib import Path
from typing import Any, Dict, List, Optional, Union

import matplotlib.pyplot as plt
from matplotlib import animation as _animation

from .. import logger
from ..analysis import Profile
from ..snap import SnapLike
from .functions import get_extent_from_percentile
from .interpolation import interpolate
from .visualization import plot

_interp_kwargs = ('number_of_pixels', 'density_weighted')


def animation(
    *,
    filename: Union[str, Path],
    snaps: List[SnapLike],
    quantity: str,
    image: Optional[Any] = None,
    text: List[str] = None,
    text_kwargs: Dict[str, Any] = {},
    func_animation_kwargs: Dict[str, Any] = {},
    save_kwargs: Dict[str, Any] = {},
    **kwargs,
):
    """Generate an animation of images.

    Parameters
    ----------
    filename
        The file name to save the animation to.
    snaps
        A list of Snap objects to animate.
    quantity
        The quantity to visualize. Must be a string to pass to Snap.
    image
        The matplotlib AxesImage object that represents the quantity to
        animate. If None, generate a new AxesImage object.
    text
        List of strings to display per snap.
    text_kwargs : optional
        Keyword arguments to pass to matplotlib text.
    func_animation_kwargs : optional
        Keyword arguments to pass to matplotlib FuncAnimation.
    save_kwargs : optional
        Keyword arguments to pass to matplotlib Animation.save.
    **kwargs
        Arguments to pass to visualize.plot.

    Returns
    -------
    anim
        The matplotlib FuncAnimation object.

    Examples
    --------
    Make an animation of multiple snaps.

    >>> plonk.visualize.animation(
    ...     snaps=snaps,
    ...     quantity='density',
    ...     extent=(-100, 100, -100, 100),
    ...     vmin=0.0,
    ...     vmax=1.0,
    ...     filename='animation.mp4',
    ... )
    """
    filepath = pathlib.Path(filename)
    if filepath.suffix != '.mp4':
        raise ValueError('filename should end in ".mp4"')

    interp = kwargs.get('interp', 'projection')
    x = kwargs.get('x', 'x')
    y = kwargs.get('x', 'y')

    if image is None:
        fig, ax = plt.subplots()
        ax = plot(snap=snaps[0], quantity=quantity, ax=ax, **kwargs)
        im = ax.images[0]
    else:
        im = image
        ax = image.axes
        fig = ax.figure

    if im is None:
        raise NotImplementedError(
            'Can only currently produce animations of images. (Or profiles\n'
            'with animation_profile.)'
        )
    if text is not None:
        _text = ax.text(
            0.9, 0.9, text[0], ha='right', transform=ax.transAxes, **text_kwargs
        )

    def animate(idx):
        logger.info(f'Visualizing snap: {idx}')
        _kwargs = {k: v for k, v in kwargs.items() if k in _interp_kwargs}
        extent = kwargs.get('extent', get_extent_from_percentile(snaps[idx], x, y))
        interp_data = interpolate(
            snap=snaps[idx], interp=interp, quantity=quantity, extent=extent, **_kwargs,
        )
        im.set_data(interp_data)
        im.set_extent(extent)
        im.axes.set_xlim(extent[:2])
        im.axes.set_ylim(extent[2:])
        vmin = kwargs.get('vmin', interp_data.min())
        vmax = kwargs.get('vmax', interp_data.max())
        im.set_clim(vmin, vmax)
        if text is not None:
            _text.set_text(text[idx])
        return [im]

    anim = _animation.FuncAnimation(
        fig, animate, frames=len(snaps), **func_animation_kwargs
    )
    anim.save(filepath, extra_args=['-vcodec', 'libx264'], **save_kwargs)
    plt.close()

    return anim


def animation_profiles(
    *,
    filename: Union[str, Path],
    profiles: List[Profile],
    quantity: str,
    line: Optional[Any] = None,
    text: List[str] = None,
    text_kwargs: Dict[str, Any] = {},
    func_animation_kwargs: Dict[str, Any] = {},
    save_kwargs: Dict[str, Any] = {},
    **kwargs,
):
    """Generate an animation of a radial profile.

    Parameters
    ----------
    filename
        The file name to save the animation to.
    profiles
        A list of Profile objects to animate.
    quantity
        The quantity to profile. Must be a string to pass to Profile.
    line
        The matplotlib Line2D object that represents the quantity to
        animate. If None, generate a new Line2D object.
    text
        List of strings to display per profile plot.
    text_kwargs : optional
        Keyword arguments to pass to matplotlib text.
    func_animation_kwargs : optional
        Keyword arguments to pass to matplotlib FuncAnimation.
    save_kwargs : optional
        Keyword arguments to pass to matplotlib Animation.save.
    **kwargs
        Arguments to pass to plot function.

    Returns
    -------
    anim
        The matplotlib FuncAnimation object.
    """
    filepath = pathlib.Path(filename)
    if filepath.suffix != '.mp4':
        raise ValueError('filename should end in ".mp4"')

    if line is None:
        fig, ax = plt.subplots()
        [line] = ax.plot(profiles[0]['radius'], profiles[0][quantity], **kwargs)
    else:
        ax = line.axes
        fig = ax.figure

    if text is not None:
        _text = ax.text(
            0.9, 0.9, text[0], ha='right', transform=ax.transAxes, **text_kwargs
        )

    def animate(idx):
        logger.info(f'Plotting profile: {idx}')
        line.set_data(profiles[idx]['radius'], profiles[idx][quantity])
        if text is not None:
            _text.set_text(text[idx])
        return [line]

    anim = _animation.FuncAnimation(
        fig, animate, frames=len(profiles), **func_animation_kwargs
    )
    anim.save(filepath, extra_args=['-vcodec', 'libx264'], **save_kwargs)
    plt.close()

    return anim
