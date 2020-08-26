"""Animations of visualizations."""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING, Any, Dict, List, Union

import matplotlib.pyplot as plt
from matplotlib import animation as _animation

try:
    from tqdm import tqdm
except ImportError:
    tqdm = None

from .._logging import logger
from .._units import units as plonk_units
from . import visualization as viz

if TYPE_CHECKING:
    from ..analysis.profile import Profile
    from ..snap.snap import SnapLike


def animate(
    filename: Union[str, Path],
    *,
    snaps: List[SnapLike] = None,
    profiles: List[Profile] = None,
    quantity: str = None,
    x: str = None,
    y: str = None,
    fig: Any = None,
    adaptive_colorbar: bool = True,
    adaptive_limits: bool = True,
    text: List[str] = None,
    text_kwargs: Dict[str, Any] = {},
    func_animation_kwargs: Dict[str, Any] = {},
    save_kwargs: Dict[str, Any] = {},
    **kwargs,
):
    """Animate a Snap or Profile visualization.

    Pass in a list of Snap objects or Profile objects. If snaps are
    passed in and the quantity is None, then the animation will be of
    particle plots, otherwise it will be of images. If profiles are
    passed in the animation will be of profiles.

    Parameters
    ----------
    filename
        The file name to save the animation to.
    snaps
        A list of Snap objects to animate.
    profiles
        A list of Profile objects to animate.
    quantity : optional
        The quantity to visualize. Must be a string to pass to Snap.
    x : optional
        The quantity for the x-axis. Must be a string to pass to Snap.
        For interpolated plots must be 'x', 'y', or 'z'.
    y : optional
        The quantity for the y-axis. Must be a string to pass to Snap.
        For interpolated plots must be 'x', 'y', or 'z'.
    fig : optional
        A matplotlib Figure object to animate. If None, generate a new
        Figure.
    adaptive_colorbar : optional
        If True, adapt colorbar range during animation. If False, the
        colorbar range is fixed by the initial image. Default is True.
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
        Arguments to pass to visualize.image.

    Returns
    -------
    anim
        The matplotlib FuncAnimation object.

    Examples
    --------
    Make an image animation of projected density.

    >>> units = {
    ...     'position': 'au',
    ...     'density': 'g/cm^3',
    ...     'projection': 'cm'
    ... }
    >>> plonk.animate(
    ...     filename='animation.mp4',
    ...     snaps=snaps,
    ...     quantity='density',
    ...     units=units,
    ...     save_kwargs={'fps': 10, 'dpi': 300},
    ... )

    Make a particle animation of x vs density.

    >>> units = {'position': 'au', 'density': 'g/cm^3'}
    >>> plonk.animate(
    ...     filename='animation.mp4',
    ...     snaps=snaps,
    ...     x='x',
    ...     y='density',
    ...     units=units,
    ...     adaptive_limits=False,
    ...     save_kwargs={'fps': 10, 'dpi': 300},
    ... )

    Make a profile animation of radius vs surface density.

    >>> units={'position': 'au', 'surface_density': 'g/cm^2'}
    >>> plonk.animate(
    ...     filename='animation.mp4',
    ...     profiles=profiles,
    ...     x='radius',
    ...     y='surface_density',
    ...     units=units,
    ...     adaptive_limits=False,
    ...     save_kwargs={'fps': 10, 'dpi': 300},
    ... )
    """
    if snaps is not None:
        if profiles is not None:
            raise ValueError('Cannot pass in snaps and profiles')
        if quantity is None:
            anim = animation_particles(
                filename=filename,
                snaps=snaps,
                x=x,
                y=y,
                fig=fig,
                adaptive_limits=adaptive_limits,
                text=text,
                text_kwargs=text_kwargs,
                func_animation_kwargs=func_animation_kwargs,
                save_kwargs=save_kwargs,
                **kwargs,
            )
        else:
            anim = animation_images(
                filename=filename,
                snaps=snaps,
                quantity=quantity,
                fig=fig,
                adaptive_colorbar=adaptive_colorbar,
                text=text,
                text_kwargs=text_kwargs,
                func_animation_kwargs=func_animation_kwargs,
                save_kwargs=save_kwargs,
                **kwargs,
            )
    elif profiles is not None:
        if x is None:
            raise ValueError('x must be passed in for profile animations')
        if y is None:
            raise ValueError('y must be passed in for profile animations')
        anim = animation_profiles(
            filename=filename,
            profiles=profiles,
            x=x,
            y=y,
            fig=fig,
            adaptive_limits=adaptive_limits,
            text=text,
            text_kwargs=text_kwargs,
            func_animation_kwargs=func_animation_kwargs,
            save_kwargs=save_kwargs,
            **kwargs,
        )
    else:
        raise ValueError('Must pass in snaps or profiles')

    return anim


def animation_images(
    *,
    filename: Union[str, Path],
    snaps: List[SnapLike],
    quantity: str,
    fig: Any = None,
    adaptive_colorbar: bool = True,
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
    fig : optional
        A matplotlib Figure object to animate. If None, generate a new
        Figure.
    adaptive_colorbar : optional
        If True, adapt colorbar range during animation. If False, the
        colorbar range is fixed by the initial image. Default is True.
    text : optional
        List of strings to display per snap.
    text_kwargs : optional
        Keyword arguments to pass to matplotlib text.
    func_animation_kwargs : optional
        Keyword arguments to pass to matplotlib FuncAnimation.
    save_kwargs : optional
        Keyword arguments to pass to matplotlib Animation.save.
    **kwargs
        Arguments to pass to visualize.image.

    Returns
    -------
    anim
        The matplotlib FuncAnimation object.

    Examples
    --------
    Make an animation of projected density.

    >>> animation_images(
    ...     filename='animation.mp4',
    ...     snaps=snaps,
    ...     quantity='density',
    ...     units={'position': 'au', 'density': 'g/cm^3'},
    ...     save_kwargs={'fps': 10, 'dpi': 300},
    ... )
    """
    filepath = Path(filename)
    if filepath.suffix != '.mp4':
        raise ValueError('filename should end in ".mp4"')

    if fig is None:
        ax = snaps[0].image(quantity=quantity, **kwargs)
        fig = ax.figure
        image = ax.images[0]
    else:
        ax = fig.axes[0]
        image = ax.images[0]

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

    x = kwargs.get('x', 'x')
    y = kwargs.get('y', 'y')
    interp = kwargs.get('interp', 'projection')
    slice_normal = kwargs.get('slice_normal')
    slice_offset = kwargs.get('slice_offset')
    extent = kwargs.get('extent')
    units = kwargs.get('units')
    weighted = kwargs.get('weighted')
    num_pixels = kwargs.get('num_pixels')

    def animate(idx):
        if tqdm is not None:
            pbar.update(n=1)
        snap = snaps[idx]
        interp_data, _extent, _units = viz._interpolated_data(
            snap=snap,
            quantity=quantity,
            x=x,
            y=y,
            interp=interp,
            slice_normal=slice_normal,
            slice_offset=slice_offset,
            extent=extent,
            units=units,
            weighted=weighted,
            num_pixels=num_pixels,
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
        return [image]

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
        Arguments to pass to Profile.plot function.

    Returns
    -------
    anim
        The matplotlib FuncAnimation object.

    Examples
    --------
    Make an animation of radius vs surface density.

    >>> animation_profiles(
    ...     filename='animation.mp4',
    ...     profiles=profiles,
    ...     x='radius',
    ...     y='surface_density',
    ...     units={'position': 'au', 'surface_density': 'g/cm^2'},
    ...     adaptive_limits=False,
    ...     save_kwargs={'fps': 10, 'dpi': 300},
    ... )
    """
    filepath = Path(filename)
    if filepath.suffix != '.mp4':
        raise ValueError('filename should end in ".mp4"')

    if isinstance(y, str):
        ys = [y]
    elif isinstance(y, (tuple, list)):
        ys = y
    else:
        raise ValueError('y must be a string or list of strings')

    prof = profiles[0]
    if fig is None:
        ax = prof.plot(x=x, y=ys, **kwargs)
        fig = ax.figure
        lines = ax.lines
    else:
        ax = fig.axes[0]
        lines = [line for line in ax.lines]

    if text is not None:
        _text = ax.text(
            0.9, 0.9, text[0], ha='right', transform=ax.transAxes, **text_kwargs
        )

    if tqdm is not None:
        pbar = tqdm(total=len(profiles))
    else:
        logger.info(
            'progress bar not available\ntry pip install tqdm --or-- conda install tqdm'
        )

    units = kwargs.get('units')
    if units is not None:
        xunit = plonk_units(units.get(prof.base_profile_name(x), str(prof[x].units)))
        yunits = [
            plonk_units(units.get(prof.base_profile_name(y), str(prof[y].units)))
            for y in ys
        ]
    else:
        xunit = prof[x].units
        yunits = [prof[y].units for y in ys]

    def animate(idx):
        if tqdm is not None:
            pbar.update(n=1)
        profile = profiles[idx]
        xlim, ylim = (0, 0), (0, 0)
        _x = (profile[x] / xunit).to_base_units().magnitude
        for line, y, yunit in zip(lines, ys, yunits):
            _y = (profile[y] / yunit).to_base_units().magnitude
            line.set_data(_x, _y)
            if adaptive_limits:
                ax = line.axes
                xlim, ylim = _get_range(_x, xlim), _get_range(_y, ylim)
                ax.set(xlim=xlim, ylim=ylim)
        if text is not None:
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
        Arguments to pass to visualize.plot.

    Returns
    -------
    anim
        The matplotlib FuncAnimation object.

    Examples
    --------
    Make an animation of x vs density on the particles.

    >>> animation_particles(
    ...     filename='animation.mp4',
    ...     snaps=snaps,
    ...     x='x',
    ...     y='density',
    ...     adaptive_limits=False,
    ...     save_kwargs={'fps': 10, 'dpi': 300},
    ... )
    """
    filepath = Path(filename)
    if filepath.suffix != '.mp4':
        raise ValueError('filename should end in ".mp4"')

    if fig is None:
        ax = snaps[0].plot(**kwargs)
        fig = ax.figure
        lines = ax.lines
    else:
        ax = fig.axes[0]
        lines = [line for line in ax.lines]

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

    x = kwargs.get('x', 'x')
    y = kwargs.get('y', 'y')
    c = kwargs.get('c')
    s = kwargs.get('s')
    units = kwargs.get('units')

    def animate(idx):
        snap = snaps[idx]
        if tqdm is not None:
            pbar.update(n=1)
        xlim, ylim = (0, 0), (0, 0)

        if c is None and s is None:
            try:
                subsnaps = snap.subsnaps_as_list(squeeze=False)
            except AttributeError:
                subsnaps = [snap]
        else:
            subsnaps = [snap]

        for line, subsnap in zip(lines, subsnaps):
            _x, _y, _c, _s, _units = viz._plot_data(
                snap=subsnap, x=x, y=y, c=c, s=s, units=units
            )
            line.set_data(_x, _y)
            if adaptive_limits:
                xlim, ylim = _get_range(_x, xlim), _get_range(_y, ylim)
                line.axes.set(xlim=xlim, ylim=ylim)
        if text is not None:
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
