"""Animations of visualizations."""

from __future__ import annotations

import pathlib
from pathlib import Path
from typing import TYPE_CHECKING, Any, Dict, List, Union

import matplotlib.pyplot as plt
from matplotlib import animation as _animation

try:
    from tqdm import tqdm
except ImportError:
    tqdm = None

from .._logging import logger
from .._units import Quantity
from .._units import units as plonk_units
from .visualization import _particle_plot_data, _plot_data, particle_plot, plot

if TYPE_CHECKING:
    from ..analysis.profile import Profile
    from ..snap.snap import SnapLike

_interp_kwargs = ('number_of_pixels', 'density_weighted')


def animation(
    *,
    filename: Union[str, Path],
    snaps: List[SnapLike],
    quantity: Union[str, List[str]],
    x: Union[str, List[str]] = 'x',
    y: Union[str, List[str]] = 'y',
    units: Dict[str, str] = None,
    extent: Union[Quantity, List[Quantity]] = None,
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
    units
        The units of the plot as a dictionary with keys 'quantity',
        'extent', 'projection'. The values are strings representing
        units, e.g. 'g/cm^3'.
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
    Make an animation of projected density.

    >>> plonk.animation(
    ...     filename='animation.mp4',
    ...     snaps=snaps,
    ...     quantity='density',
    ...     units={'extent': 'au', 'quantity': 'g/cm^3'},
    ...     adaptive_colorbar=False,
    ...     save_kwargs={'fps': 10, 'dpi': 300},
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
    elif extent is None:
        extents = [None for _ in quantities]
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
                units=units,
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
        snap = snaps[idx]
        interp_kwargs = {k: v for k, v in kwargs.items() if k in _interp_kwargs}
        interp = kwargs.get('interp', 'projection')
        slice_normal = kwargs.get('slice_normal')
        z_slice = kwargs.get('z_slice')
        for quantity, x, y, extent, image in zip(quantities, xs, ys, extents, images):
            interp_data, _extent, _units = _plot_data(
                snap=snap,
                quantity=quantity,
                x=x,
                y=y,
                interp=interp,
                slice_normal=slice_normal,
                z_slice=z_slice,
                extent=extent,
                units=units,
                **interp_kwargs,
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
    x: str,
    y: Union[str, List[str]],
    units: Dict[str, Any] = None,
    fig: Any = None,
    adaptive_limits: bool = True,
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
    x
        The quantity for the x-axis. Must be a string to pass to Snap.
    y
        The quantity for the y-axis, or list of quantities. Must be a
        string (or list of strings) to pass to Snap.
    units
        The units of the plot as a dictionary with keys 'x', 'y'.
        The values are strings representing units, e.g. 'g/cm^3'.
    fig : optional
        A matplotlib Figure object to animate. If None, generate a new
        Figure.
    adaptive_limits : optional
        If True, adapt plot limits during animation. If False, the plot
        limits are fixed by the initial plot. Default is True.
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

    Examples
    --------
    Make an animation of radius vs surface density.

    >>> plonk.animation_profiles(
    ...     filename='animation.mp4',
    ...     profiles=profiles,
    ...     x='radius',
    ...     y='surface_density',
    ...     units={'x': 'au', 'y': 'g/cm^2'},
    ...     adaptive_limits=False,
    ...     save_kwargs={'fps': 10, 'dpi': 300},
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
        raise ValueError('Cannot determine quantity')

    if units is None:
        _units = {'x': profiles[0][x].units, 'y': [profiles[0][y].units for y in ys]}
    else:
        _units = units
    if isinstance(_units['y'], str):
        yu = _units['y']
        _units['y'] = [yu for _ in range(len(ys))]
    units_str = _units
    units_val = {
        'x': plonk_units(_units['x']),
        'y': [plonk_units(yu) for yu in _units['y']],
    }

    if fig is None:
        fig, ax = plt.subplots()
        for yidx, y in enumerate(ys):
            __units = {'x': units_str['x'], 'y': units_str['y'][yidx]}
            profiles[0].plot(x=x, y=y, units=__units, ax=ax, **kwargs)
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
        xlim, ylim = (0, 0), (0, 0)
        _x = (profiles[idx][x] / units_val['x']).to_base_units().magnitude
        for yidx, (line, y) in enumerate(zip(lines, ys)):
            _y = (profiles[idx][y] / units_val['y'][yidx]).to_base_units().magnitude
            line.set_data(_x, _y)
            if adaptive_limits:
                ax = line.axes
                xlim, ylim = _get_range(_x, xlim), _get_range(_y, ylim)
                ax.set(xlim=xlim, ylim=ylim)
        if text is not None:
            for _text in texts:
                _text.set_text(text[idx])
        return lines

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
    units: Dict[str, Any] = None,
    fig: Any = None,
    adaptive_limits: bool = True,
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
    units
        The units of the plot as a dictionary with keys 'x', 'y'.
        The values are strings representing units, e.g. 'g/cm^3'.
    fig : optional
        A matplotlib Figure object to animate. If None, generate a new
        Figure.
    adaptive_limits : optional
        If True, adapt plot limits during animation. If False, the plot
        limits are fixed by the initial plot. Default is True.
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
    ...     filename='animation.mp4',
    ...     snaps=snaps,
    ...     x='x',
    ...     y='density',
    ...     adaptive_limits=False,
    ...     save_kwargs={'fps': 10, 'dpi': 300},
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

    if units is None:
        _units = {'x': snaps[0][x].units, 'y': [snaps[0][y].units for y in ys]}
    else:
        _units = units
    if isinstance(_units['y'], str):
        yu = _units['y']
        _units['y'] = [yu for _ in range(len(ys))]
    units_str = _units

    if fig is None:
        fig, ax = plt.subplots()
        for yidx, y in enumerate(ys):
            __units = {'x': units_str['x'], 'y': units_str['y'][yidx]}
            particle_plot(snap=snaps[0], x=x, y=y, ax=ax, units=__units, **kwargs)
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
        snap = snaps[idx]
        subsnaps = snap.subsnaps_as_list()
        if tqdm is not None:
            pbar.update(n=1)
        xlim, ylim = (0, 0), (0, 0)
        for yidx, y in enumerate(ys):
            __units = {'x': units_str['x'], 'y': units_str['y'][yidx]}
            for line, subsnap in zip(
                lines[yidx * len(subsnaps) : (yidx + 1) * len(subsnaps)], subsnaps
            ):
                _x, _y, _, _, _ = _particle_plot_data(
                    snap=subsnap, x=x, y=y, c=None, s=None, units=__units
                )
                line.set_data(_x, _y)
                if adaptive_limits:
                    xlim, ylim = _get_range(_x, xlim), _get_range(_y, ylim)
                    line.axes.set(xlim=xlim, ylim=ylim)
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
