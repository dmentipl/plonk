"""Animations of visualizations."""

import pathlib
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

import matplotlib.pyplot as plt
from matplotlib import animation as _animation

from .. import logger
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
    line : optional
        The matplotlib Line2D object that represents the quantity to
        animate. If None, generate a new Line2D object.
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


def animation_particles(
    *,
    filename: Union[str, Path],
    snaps: List[SnapLike],
    x: str,
    y: Union[str, List[str]],
    xlim: Tuple[float, float] = None,
    ylim: Tuple[float, float] = None,
    lines: List[Any] = None,
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
    xlim : optional
        The x-axis range as (xmin, xmax).
    ylim : optional
        The y-axis range as (ymin, ymax).
    lines : optional
        A list of matplotlib Line2D object that represents the quantity
        to animate. If None, generate a new Line2D object.
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
    ...     xlim=(-100, 100),
    ...     ylim=(0, 1),
    ...     filename='animation.mp4',
    ... )

    Make an animation of x vs density and x vs velocity_x on the
    particles on the same plot.

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

    if isinstance(y, (list, tuple)):
        ys = y
    elif isinstance(y, str):
        ys = [y]
    else:
        raise ValueError('Cannot determine y')

    if lines is None:
        fig, ax = plt.subplots()
        for y in ys:
            particle_plot(snap=snaps[0], x=x, y=y, ax=ax, **kwargs)
        ax.set(xlim=xlim, ylim=ylim)
        lines = ax.lines
    else:
        ax = lines[0].axes
        fig = ax.figure

    if text is not None:
        _text = ax.text(
            0.9, 0.9, text[0], ha='right', transform=ax.transAxes, **text_kwargs
        )

    def animate(idx):
        logger.info(f'Visualizing snap: {idx}')
        subsnaps = snaps[idx].subsnaps_as_list()
        num_subsnaps = len(subsnaps)
        _xlim, _ylim = [0, 0], [0, 0]
        for idxi, y in enumerate(ys):
            for idxj, subsnap in enumerate(snaps[idx].subsnaps_as_list()):
                _x = subsnap[x]
                _y = subsnap[y]
                ax.lines[idxi * num_subsnaps + idxj].set_data(_x, _y)
                if xlim is None:
                    _xlim = _get_range(_x, _xlim)
                    ax.set_xlim(_xlim)
                if ylim is None:
                    _ylim = _get_range(_y, _ylim)
                    ax.set_ylim(_ylim)
        if text is not None:
            _text.set_text(text[idx])
        return lines

    anim = _animation.FuncAnimation(
        fig, animate, frames=len(snaps), **func_animation_kwargs
    )
    anim.save(filepath, extra_args=['-vcodec', 'libx264'], **save_kwargs)
    plt.close()

    return anim


def _get_range(vals, vlim):
    vmin, vmax = vals.min(), vals.max()
    dv = 0.05 * (vmax - vmin)
    return (
        vmin - dv if vmin - dv < vlim[0] else vlim[0],
        vmax + dv if vmax + dv > vlim[1] else vlim[1],
    )
