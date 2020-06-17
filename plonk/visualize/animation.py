"""Animations of visualizations."""

import pathlib
from pathlib import Path
from typing import Any, Dict, List, Tuple, Union

import matplotlib.pyplot as plt
from matplotlib import animation as _animation
from tqdm import tqdm

from ..analysis import Profile
from ..snap import SnapLike
from .functions import get_extent_from_percentile
from .interpolation import interpolate
from .visualization import particle_plot, plot

_interp_kwargs = ('number_of_pixels', 'density_weighted')


def animation(
    *,
    filename: Union[str, Path],
    snaps: List[SnapLike],
    quantity: str,
    image: Any = None,
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
    image : optional
        The matplotlib AxesImage object that represents the quantity to
        animate. If None, generate a new AxesImage object.
    text : optional
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
            'with animation_profile or particles with animation_particles.)'
        )
    if text is not None:
        _text = ax.text(
            0.9, 0.9, text[0], ha='right', transform=ax.transAxes, **text_kwargs
        )

    pbar = tqdm(total=len(snaps))

    def animate(idx):
        pbar.update(n=1)
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
    pbar.close()

    return anim


def animation_profiles(
    *,
    filename: Union[str, Path],
    profiles: List[Profile],
    quantity: Union[str, List[str]],
    fig: Any = None,
    adaptive_limits: bool = True,
    xlim: Tuple[float, float] = None,
    ylim: Tuple[float, float] = None,
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
        The quantity, or quantities, to profile. Must be a string to
        pass to Profile, or a list of such strings.
    fig : optional
        A matplotlib Figure object to animate. If None, generate a new
        Figure.
    adaptive_limits : optional
        If True, adapt plot limits during animation. If False, the plot
        limits are fixed by the initial plot.
    xlim : optional
        The x-axis range as a tuple (xmin, xmax). Only used if fig not
        passed in.
    ylim : optional
        The y-axis range as a tuple (ymin, ymax). Only used if fig not
        passed in.
    text : optional
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

    if isinstance(quantity, list):
        quantities = quantity
    elif isinstance(quantity, str):
        quantities = [quantity]
    else:
        raise ValueError('Cannot determine quantity')

    if fig is None:
        fig, ax = plt.subplots()
        for quantity in quantities:
            ax.plot(profiles[0]['radius'], profiles[0][quantity], **kwargs)
            ax.set(xlim=xlim, ylim=ylim)
        lines = ax.lines
        if text is not None:
            texts = [
                ax.text(
                    0.9, 0.9, text[0], ha='right', transform=ax.transAxes, **text_kwargs
                )
            ]
    else:
        lines = [line for ax in fig.axes for line in ax.lines]
        if text is not None:
            texts = [text for ax in fig.axes for text in ax.texts]

    pbar = tqdm(total=len(profiles))

    def animate(idx):
        pbar.update(n=1)
        for line, quantity in zip(lines, quantities):
            x, y = profiles[idx]['radius'], profiles[idx][quantity]
            line.set_data(x, y)
            if adaptive_limits:
                ax = line.axes
                ylim = _get_range(y, (0, 0))
                ax.set(ylim=ylim)
        if text is not None:
            for _text in texts:
                _text.set_text(text[idx])
        return [line]

    anim = _animation.FuncAnimation(
        fig, animate, frames=len(profiles), **func_animation_kwargs
    )
    anim.save(filepath, extra_args=['-vcodec', 'libx264'], **save_kwargs)
    plt.close()
    pbar.close()

    return anim


def animation_particles(
    *,
    filename: Union[str, Path],
    snaps: List[SnapLike],
    x: str,
    y: Union[str, List[str]],
    fig: Any = None,
    adaptive_limits: bool = True,
    xlim: Tuple[float, float] = None,
    ylim: Tuple[float, float] = None,
    text: List[str] = None,
    text_kwargs: Dict[str, Any] = {},
    func_animation_kwargs: Dict[str, Any] = {},
    save_kwargs: Dict[str, Any] = {},
    **kwargs,
):
    """Generate an animation of particle plots.

    Parameters
    ----------
    filename
        The file name to save the animation to.
    snaps
        A list of Snap objects to animate.
    x
        The quantity for the x-axis. Must be a string to pass to Snap.
    y
        The quantity for the y-axis, or list of quantities. Must be a
        string (or list of strings) to pass to Snap.
    fig : optional
        A matplotlib Figure object to animate. If None, generate a new
        Figure.
    adaptive_limits : optional
        If True, adapt plot limits during animation. If False, the plot
        limits are fixed by the initial plot.
    xlim : optional
        The x-axis range as a tuple (xmin, xmax). Only used if fig not
        passed in.
    ylim : optional
        The y-axis range as a tuple (ymin, ymax). Only used if fig not
        passed in.
    text : optional
        List of strings to display per snap.
    text_kwargs : optional
        Keyword arguments to pass to matplotlib text.
    func_animation_kwargs : optional
        Keyword arguments to pass to matplotlib FuncAnimation.
    save_kwargs : optional
        Keyword arguments to pass to matplotlib Animation.save.
    **kwargs
        Arguments to pass to visualize.particle_plot.

    Returns
    -------
    anim
        The matplotlib FuncAnimation object.

    Examples
    --------
    Make an animation of x vs density on the particles.

    >>> plonk.visualize.animation(
    ...     snaps=snaps,
    ...     x='x',
    ...     y='density',
    ...     filename='animation.mp4',
    ... )

    Make an animation of x vs density x vs velocity_x on the
    particles on the same axes, specifying the x- and y-plot limits.

    >>> plonk.visualize.animation(
    ...     snaps=snaps,
    ...     x='x',
    ...     y=['density', 'velocity_x'],
    ...     xlim=(-100, 100),
    ...     ylim=(0, 1),
    ...     filename='animation.mp4',
    ... )
    """
    filepath = pathlib.Path(filename)
    if filepath.suffix != '.mp4':
        raise ValueError('filename should end in ".mp4"')

    if isinstance(y, list):
        ys = y
    elif isinstance(y, str):
        ys = [y]
    else:
        raise ValueError('Cannot determine y')

    if fig is None:
        fig, ax = plt.subplots()
        for y in ys:
            particle_plot(snap=snaps[0], x=x, y=y, ax=ax, **kwargs)
            ax.set(xlim=xlim, ylim=ylim)
        lines = ax.lines
        if text is not None:
            texts = [
                ax.text(
                    0.9, 0.9, text[0], ha='right', transform=ax.transAxes, **text_kwargs
                )
            ]
    else:
        lines = [line for ax in fig.axes for line in ax.lines]
        if text is not None:
            texts = [text for ax in fig.axes for text in ax.texts]

    pbar = tqdm(total=len(snaps))

    def animate(idx):
        pbar.update(n=1)
        subsnaps = snaps[idx].subsnaps_as_list()
        num_subsnaps = len(subsnaps)
        for idxi, y in enumerate(ys):
            _xlim, _ylim = (0.0, 0.0), (0.0, 0.0)
            for idxj, subsnap in enumerate(subsnaps):
                line = lines[idxi * num_subsnaps + idxj]
                _x = subsnap[x]
                _y = subsnap[y]
                line.set_data(_x, _y)
                if adaptive_limits:
                    ax = line.axes
                    _xlim, _ylim = _get_range(_x, _xlim), _get_range(_y, _ylim)
                    ax.set(xlim=_xlim, ylim=_ylim)
        if text is not None:
            for _text in texts:
                _text.set_text(text[idx])
        return lines

    anim = _animation.FuncAnimation(
        fig, animate, frames=len(snaps), **func_animation_kwargs
    )
    anim.save(filepath, extra_args=['-vcodec', 'libx264'], **save_kwargs)
    plt.close()
    pbar.close()

    return anim


def _get_range(vals, vlim):
    vmin, vmax = vals.min(), vals.max()
    dv = 0.05 * (vmax - vmin)
    return (
        vmin - dv if vmin - dv < vlim[0] else vlim[0],
        vmax + dv if vmax + dv > vlim[1] else vlim[1],
    )
