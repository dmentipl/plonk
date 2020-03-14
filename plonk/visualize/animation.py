"""Animations of visualizations."""

import pathlib
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

import matplotlib.pyplot as plt
from matplotlib import animation as _animation
from numpy import ndarray

from ..analysis import Profile
from ..snap import SnapLike
from ..utils import get_extent_from_percentile
from .visualization import interpolate, plot

_interp_kwargs = ('number_of_pixels', 'density_weighted')


def animation(
    *,
    filename: Union[str, Path],
    snaps: List[SnapLike],
    quantity: Optional[Union[str, ndarray]] = None,
    text: List[str] = None,
    text_kwargs: Dict[str, Any] = {},
    func_animation_kwargs: Dict[str, Any] = {},
    save_kwargs: Dict[str, Any] = {},
    **kwargs,
):
    """Generate an animation.

    Parameters
    ----------
    filename
        The file name to save the animation to.
    snaps
        A list of Snap objects to animate.
    quantity
        The quantity to visualize. Can be a string to pass to Snap, or
        a 1d array (N,) of scalar data, or a 2d array (N, 3) of
        vector data. If quantity is 2d, only the first two components
        are interpolated, i.e. quantity[:, 0] and quantity[:, 1].
        Default is None.
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
    """
    filepath = pathlib.Path(filename)
    if filepath.suffix != '.mp4':
        raise ValueError('filename should end in ".mp4"')

    interp = kwargs.get('interp', 'projection')

    fig, ax = plt.subplots()
    viz = plot(snap=snaps[0], quantity=quantity, ax=ax, **kwargs)
    im = viz.objects['image']
    if im is None:
        raise NotImplementedError(
            'Can only currently produce animations of rendered images. (Or profiles\n'
            'with animation_profile.)'
        )
    if text is not None:
        _text = ax.text(
            0.9, 0.9, text[0], ha='right', transform=ax.transAxes, **text_kwargs
        )

    def animate(idx):
        print(f'Visualizing snap: {idx}')
        _kwargs = {k: v for k, v in kwargs.items() if k in _interp_kwargs}
        extent = kwargs.get('extent', get_extent_from_percentile(snaps[idx]))
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
    text: List[str] = None,
    text_kwargs: Dict[str, Any] = {},
    func_animation_kwargs: Dict[str, Any] = {},
    save_kwargs: Dict[str, Any] = {},
    **kwargs,
):
    """Generate an animation.

    Parameters
    ----------
    filename
        The file name to save the animation to.
    profiles
        A list of Profile objects to animate.
    quantity
        The quantity to profile. Can be a string to pass to Profile.
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
    """
    filepath = pathlib.Path(filename)
    if filepath.suffix != '.mp4':
        raise ValueError('filename should end in ".mp4"')

    fig, ax = plt.subplots()
    [line] = ax.plot(profiles[0]['radius'], profiles[0][quantity], **kwargs)
    if text is not None:
        _text = ax.text(
            0.9, 0.9, text[0], ha='right', transform=ax.transAxes, **text_kwargs
        )

    def animate(idx):
        print(f'Plotting profile: {idx}')
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
