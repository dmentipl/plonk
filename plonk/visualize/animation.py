"""Animations of visualizations."""

from __future__ import annotations

import pathlib
from pathlib import Path
from typing import TYPE_CHECKING, Any, Dict, List, Tuple, Union

import matplotlib.pyplot as plt
from matplotlib import animation as _animation

try:
    from tqdm import tqdm
except ImportError:
    tqdm = None

from .._logging import logger
from .functions import get_extent_from_percentile
from .interpolation import interpolate
from .visualization import particle_plot, plot

if TYPE_CHECKING:
    from ..analysis.profile import Profile
    from ..snap.snap import SnapLike

_interp_kwargs = ('number_of_pixels', 'density_weighted')

Extent = Tuple[float, float, float, float]


def animation(
    *,
    filename: Union[str, Path],
    snaps: List[SnapLike],
    quantity: Union[str, List[str]],
    x: Union[str, List[str]] = 'x',
    y: Union[str, List[str]] = 'y',
    extent: Union[Extent, List[Extent]] = (-1, -1, -1, -1),
    fig: Any = None,
    adaptive_colorbar: bool = False,
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
        The quantity, or quantities, to visualize. Must be a string, or
        list of strings, to pass to Snap.
    x
        The x-axis as a string to pass to Snap. Can be a list of strings
        if multiple quantities are animated.
    y
        The y-axis as a string to pass to Snap. Can be a list of strings
        if multiple quantities are animated.
    extent
        The extent as a tuple (xmin, xmax, ymin, ymax). Can be a list
        of tuples if multiple quantities are animated.
    fig : optional
        A matplotlib Figure object to animate. If None, generate a new
        Figure.
    adaptive_colorbar : optional
        If True, adapt colorbar range during animation. If False, the
        colorbar range is fixed by the initial plot. Default is False.
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

    >>> plonk.animation(
    ...     snaps=snaps,
    ...     quantity='density',
    ...     x='x',
    ...     y='y',
    ...     extent=(-100, 100, -100, 100),
    ...     vmin=0.0,
    ...     vmax=1.0,
    ...     filename='animation.mp4',
    ... )
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
    if isinstance(x, list):
        xs = x
    elif isinstance(x, str):
        xs = [x]
    else:
        raise ValueError('Cannot determine x')
    if isinstance(y, list):
        ys = y
    elif isinstance(y, str):
        ys = [y]
    else:
        raise ValueError('Cannot determine y')
    if isinstance(extent, list):
        extents = extent
    elif isinstance(extent, tuple):
        extents = [extent]
    else:
        raise ValueError('Cannot determine extent')
    if len(quantities) > 1:
        if len(xs) == 1:
            xs = [xs[0] for _ in range(len(quantities))]
        if len(ys) == 1:
            ys = [ys[0] for _ in range(len(quantities))]
        if len(extents) == 1:
            extents = [extents[0] for _ in range(len(quantities))]
    if len(xs) != len(quantities):
        raise ValueError('x must have same length as quantity or be a single value')
    if len(ys) != len(quantities):
        raise ValueError('y must have same length as quantity or be a single value')
    if len(extents) != len(quantities):
        raise ValueError('extent have same length as quantity or be a single value')

    interp = kwargs.get('interp', 'projection')

    if fig is None:
        fig, axs = plt.subplots(ncols=len(quantities), squeeze=False)
        axs = axs.flatten()
        images = list()
        for quantity, x, y, extent, ax in zip(quantities, xs, ys, extents, axs):
            plot(
                snap=snaps[0],
                quantity=quantity,
                x=x,
                y=y,
                extent=extent,
                ax=ax,
                **kwargs,
            )
            images += ax.images
    else:
        images = [image for ax in fig.axes for image in ax.images]

    if text is not None:
        _text = ax.text(
            0.9, 0.9, text[0], ha='right', transform=ax.transAxes, **text_kwargs
        )

    if tqdm is not None:
        pbar = tqdm(total=len(snaps))
    else:
        logger.info(
            'progress bar not available\ntry pip install tqdm --or-- conda install tqdm'
        )

    def animate(idx):
        if tqdm is not None:
            pbar.update(n=1)
        _kwargs = {k: v for k, v in kwargs.items() if k in _interp_kwargs}
        for quantity, x, y, extent, image in zip(quantities, xs, ys, extents, images):
            if extent == (-1, -1, -1, -1):
                _extent = get_extent_from_percentile(snaps[idx], x, y)
            else:
                _extent = extent
            interp_data = interpolate(
                snap=snaps[idx],
                interp=interp,
                quantity=quantity,
                x=x,
                y=y,
                extent=_extent,
                **_kwargs,
            )
            image.set_data(interp_data)
            image.set_extent(_extent)
            image.axes.set_xlim(_extent[:2])
            image.axes.set_ylim(_extent[2:])
            if adaptive_colorbar:
                vmin = kwargs.get('vmin', interp_data.min())
                vmax = kwargs.get('vmax', interp_data.max())
                image.set_clim(vmin, vmax)
        if text is not None:
            _text.set_text(text[idx])
        return images

    anim = _animation.FuncAnimation(
        fig, animate, frames=len(snaps), **func_animation_kwargs
    )
    anim.save(filepath, extra_args=['-vcodec', 'libx264'], **save_kwargs)
    plt.close()
    if tqdm is not None:
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
    """Generate an animation of a profile.

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
        limits are fixed by the initial plot. Default is True.
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

    if tqdm is not None:
        pbar = tqdm(total=len(profiles))
    else:
        logger.info(
            'progress bar not available\ntry pip install tqdm --or-- conda install tqdm'
        )

    def animate(idx):
        if tqdm is not None:
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
    if tqdm is not None:
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
        limits are fixed by the initial plot. Default is True.
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

    >>> plonk.animation(
    ...     snaps=snaps,
    ...     x='x',
    ...     y='density',
    ...     filename='animation.mp4',
    ... )

    Make an animation of x vs density x vs velocity_x on the
    particles on the same axes, specifying the x- and y-plot limits.

    >>> plonk.animation(
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

    if tqdm is not None:
        pbar = tqdm(total=len(snaps))
    else:
        logger.info(
            'progress bar not available\ntry pip install tqdm --or-- conda install tqdm'
        )

    def animate(idx):
        if tqdm is not None:
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
    if tqdm is not None:
        pbar.close()

    return anim


def _get_range(vals, vlim):
    vmin, vmax = vals.min(), vals.max()
    dv = 0.05 * (vmax - vmin)
    return (
        vmin - dv if vmin - dv < vlim[0] else vlim[0],
        vmax + dv if vmax + dv > vlim[1] else vlim[1],
    )
