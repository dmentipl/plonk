"""
image.py

Daniel Mentiplay, 2019.
"""

import matplotlib as mpl
import matplotlib.colors as colors
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable

from ..particles import I_DUST, I_GAS
from ..utils import normalize_vector, rotate_vector_arbitrary_axis
from .interpolation import scalar_interpolation, vector_interpolation

# --- Set default plot option dictionary.
plot_options = {
    'accelerate': False,
    'colorbar': True,
    'colormap': 'gist_heat',
    'colorscale': 'linear',
    'density_weighted': False,
    'dscreen': None,
    'fontfamily': 'sans-serif',
    'fontsize': 12,
    'normalize': False,
    'zobserver': None,
}


# --- Main plotting function.
def plot(
    dump,
    render=None,
    vector=None,
    particle_types=None,
    itype=None,
    horizontal_range=None,
    vertical_range=None,
    rotation_axis=None,
    rotation_angle=None,
    position_angle=None,
    inclination=None,
    stream=None,
    cross_section=None,
    zslice=None,
    opacity=None,
    number_pixels=None,
    image_range=-1,
    render_scale=None,
    render_min=None,
    render_max=None,
    render_fraction_max=None,
    title=None,
    ax=None,
    newfig=None,
    colorbar=None,
    colormap=None,
):
    """
    Visualize a dump as a particle plot, a rendered plot, or a vector
    plot.

    Parameters
    ----------
    dump : Dump object
        The Plonk ``Dump`` object to visualize.
    render : str, default ``None``
        Scalar quantity to render. Associated options are:
        ``cross_section``, ``zslice``, ``opacity``, ``render_scale``,
        ``render_min``, ``render_max``, ``render_fraction_max``.
    vector : str, default ``None``
        Vector quantity to be represented as arrows or stream function.
        See also: ``stream``.
    particle_types : list of str, default ``None``
        Particle types to plot. See also: ``itype``.
    itype : int, default ``None``
        Particle type to plot, represented as an integer type.
    horizontal_range : list of float (len=2), default ``None``
        The range of values for the horizontal (x) axis.
    vertical_range : list of float (len=2), default ``None``
        The range of values for the vertical (y) axis.
    rotation_axis : list of float (len=3), default ``None``
        A 3-dimensional vector specifying an axis around which to
        rotate the reference frame. Must also specify ``rotation_angle``.
    rotation_angle : float, default ``None``
        An angle (radians) to rotate the reference frame around a
        ``rotation_axis``.
    position_angle : float, default ``None``
        An angle (radians) East of North specifying an axis around
        which to incline the reference frame as specified by
        ``inclination``.
    inclination : float, default ``None``
        An angle (radians) of inclination specified relative to a
        ``position_angle``.
    stream : bool, default ``None``
        If true, plot vector plots as stream functions.
    cross_section : bool, default ``None``
        If true, plot a cross section rather than column density with
        slice thickness specified by ``zslice``.
    zslice : float, default ``None``
        Thickness of a cross sectional slice.
    opacity : float, default ``None``
        Opacity for opacity plots.
    number_pixels : list of float (len=2), default ``None``
        The number of pixels in the horizontal and vertical directions,
        like [npixx, npixy], for interpolation. This determines the
        resolution of the image.
    image_range : float, default ``None``
        Specify the horizontal and vertical image range to both equal
        ``image_range``.
    render_scale : str, default ``None``
        Render scale options include: 'log', 'linear'.
    render_min : float, default ``None``
        Minimum value of the rendered quantity.
    render_max : float, default ``None``
        Maximum value of the rendered quantity.
    render_fraction_max : float, default ``None``
        Maximum value of the rendered quantity specified as a fraction
        of the maximum in the data.
    title : str, default ``None``
        Plot title.
    ax : Axes, default ``None``
        Matplotlib Axes object to plot to.
    newfig : bool, default ``None``
        If true, create new figure using Matplotlib.
    colorbar : bool, default ``None``
        If true, plot colorbar.
    colormap : str, default ``None``
        Specify the colormap.

    Returns
    -------

    Examples
    --------
    Rendering density.
    >>> plonk.plot(dump, render='rho')

    Plotting velocity vectors.
    >>> plonk.plot(dump, vector='velocity')
    """

    # TODO: add options
    # TODO: add docs
    # TODO: choose fluid type: gas, dust1, dust2, ...
    # TODO: add more checks on input
    # TODO: physical units
    # TODO: calculated quantities

    # --- Render and vector field

    particles = False
    if render is None and vector is None:
        particles = True

    if render is not None and render not in dump.particles:
        raise ValueError(f'{render} not available for rendering')

    if vector is not None and vector not in ['v', 'velocity']:
        raise ValueError(f'{vector} not available for vector field overlay')

    if stream is None:
        stream = False

    # --- Rotate frame

    rotate = False

    if (rotation_axis is not None or rotation_angle is not None) and (
        position_angle is not None and inclination is not None
    ):
        raise ValueError(
            'Cannot set rotation_axis/rotation_angle and '
            + ' position_angle/inclination at the same time'
        )

    if rotation_axis is not None:
        if rotation_angle is not None:
            rotate = True
            rotation_axis = normalize_vector(rotation_axis)
        else:
            raise ValueError('Must specify rotation_angle')
        if isinstance(rotation_axis, list):
            rotation_axis = np.array(rotation_axis)
        rotation_axis = normalize_vector(rotation_axis)

    if rotation_angle is not None and rotation_axis is None:
        raise ValueError('Must specify rotation_axis')

    if position_angle is not None:
        if inclination is not None:
            rotate = True
            rotation_angle = inclination
            rotation_axis = np.array(
                [np.cos(position_angle), np.sin(position_angle), 0]
            )
        else:
            raise ValueError('Must specify inclination')

    if inclination is not None and position_angle is None:
        raise ValueError('Must specify position_angle')

    # --- Render type

    if cross_section is None:
        cross_section = False
    if cross_section:
        if zslice is None:
            zslice = 0.0

    if opacity is None:
        opacity = False

    # --- Particle type

    itypes = list()

    if itype is not None and particle_types is not None:
        raise ValueError('Cannot set itype and particle_types together')

    if particle_types is None and itype is None:
        itypes = [I_GAS]

    if itype is not None:
        itypes = [itype]

    if particle_types is not None:

        if 'gas' in particle_types:
            itypes.append(I_GAS)

        # TODO: can only plot one type at the moment
        if 'dust' in particle_types:
            for i in range(dump.parameters['ndustlarge']):
                itypes.append(I_DUST + i)

    if len(itypes) > 1:
        raise ValueError('plotting multiple types at once is not working')

    # --- Number of pixels

    if number_pixels is None:
        npix = [512, 512]
    else:
        npix = number_pixels

    # --- Get options

    normalize = plot_options['normalize']
    zobserver = plot_options['zobserver']
    dscreen = plot_options['dscreen']
    accelerate = plot_options['accelerate']

    # --- Dataframe subsets

    pd = dump.particles.loc[dump.particles['itype'].isin(itypes)].copy()

    positions = np.array(pd[['x', 'y', 'z']])
    smoothing_length = np.array(pd['h'])
    particle_mass = np.array(pd['m'])

    if 'vx' in pd:
        velocities = np.array(pd[['vx', 'vy', 'vz']])

    if render:
        render_data = np.array(pd[render])

    if vector:
        vector_data = velocities

    # --- Interpolation weights

    weights = _interpolation_weights(
        plot_options['density_weighted'], pd, dump.parameters['hfact']
    )

    # --- Rotate frame

    if rotate:
        print(
            f'Rotating {rotation_angle*180/np.pi:.0f} deg around '
            f'[{rotation_axis[0]:.2f},'
            f' {rotation_axis[1]:.2f},'
            f' {rotation_axis[2]:.2f}]'
        )
        positions, velocities = _rotate_frame(
            positions, velocities, rotation_axis, rotation_angle
        )

    # --- Image window range

    if image_range > 0:
        if horizontal_range is not None or vertical_range is not None:
            raise ValueError(
                'Cannot set image_range and horizontal_range '
                + '(or vertical_range) at the same time'
            )
        horizontal_range = [-image_range, image_range]
        vertical_range = [-image_range, image_range]

    if horizontal_range is None and vertical_range is None:
        range_min = min(positions[:, 0].min(), positions[:, 1].min())
        range_max = max(positions[:, 0].max(), positions[:, 1].max())
        horizontal_range = [range_min, range_max]
        vertical_range = [range_min, range_max]

    if horizontal_range is None:
        horizontal_range = [positions[:, 0].min(), positions[:, 0].max()]

    if vertical_range is None:
        vertical_range = [positions[:, 1].min(), positions[:, 1].max()]

    extent = horizontal_range + vertical_range

    # --- Interpolate scalar data

    if render:

        print(f'Rendering {render}')

        image_data = scalar_interpolation(
            positions,
            smoothing_length,
            weights,
            render_data,
            particle_mass,
            horizontal_range,
            vertical_range,
            npix,
            cross_section,
            zslice,
            opacity,
            normalize,
            zobserver,
            dscreen,
            accelerate,
        )

    # --- Interpolate vector data

    if vector:

        print(f'Vector field {vector}')

        vector_data = vector_interpolation(
            positions,
            smoothing_length,
            weights,
            vector_data,
            horizontal_range,
            vertical_range,
            npix,
            cross_section,
            zslice,
            normalize,
            zobserver,
            dscreen,
        )

        xvector_data = vector_data[0]
        yvector_data = vector_data[1]

    # --- Physical units

    ############################################################################
    # TODO: temporary; testing phase
    physical_units = False

    if physical_units:
        extent = [val * dump.units['distance'] for val in extent]
        horizontal_range = [
            val * dump.units['distance'] for val in horizontal_range
        ]
        vertical_range = [
            val * dump.units['distance'] for val in vertical_range
        ]
        image_data *= dump.units['surface_density']
    ############################################################################

    # --- Render settings

    if render:

        cmap = plot_options['colormap']
        if colormap is not None:
            cmap = colormap

        if render_max is None:
            vmax = image_data.max()

        if render_fraction_max is not None:
            vmax = image_data.max() * render_fraction_max

        if render_min is None:
            vmin = image_data.min()

        if render_scale is None:
            render_scale = plot_options['colorscale']

        if render_scale == 'log':
            norm = colors.Sym_log_norm(1e-1 * vmax, clip=True)
        elif render_scale == 'linear':
            norm = colors.Normalize(vmin=vmin, vmax=vmax, clip=True)
        else:
            raise ValueError("Unknown color render_scale: " + render_scale)

    # --- Font settings

    mpl.rcParams['font.family'] = plot_options['fontfamily']
    mpl.rcParams['font.size'] = plot_options['fontsize']

    # --- Figure and axis handles

    if ax is None:
        if newfig:
            plt.figure()
        plt.clf()
        ax = plt.gca()

    # --- Rendered image

    if render:

        img = ax.imshow(
            image_data, norm=norm, origin='lower', extent=extent, cmap=cmap
        )

    # --- Plot particles

    if particles:

        print('Plotting particles')
        marker_size = 0.01
        ax.scatter(positions[:, 0], positions[:, 1], s=marker_size, c='k')
        ax.set_aspect('equal', 'box')

    # --- Vector field

    if vector:

        X, Y = np.meshgrid(
            np.linspace(*horizontal_range, len(xvector_data)),
            np.linspace(*vertical_range, len(yvector_data)),
        )

        ########################################################################
        # TODO: temporary; testing phase
        # set stride such that number of vector arrows is ~ 15-25 x 15-25
        stride = 25

        vector_color = 'black'
        if render:
            vector_color = 'white'
        ########################################################################

        if stream:
            ax.streamplot(
                X[::stride, ::stride],
                Y[::stride, ::stride],
                xvector_data[::stride, ::stride],
                yvector_data[::stride, ::stride],
                color=vector_color,
            )

        else:
            ax.quiver(
                X[::stride, ::stride],
                Y[::stride, ::stride],
                xvector_data[::stride, ::stride],
                yvector_data[::stride, ::stride],
                color=vector_color,
            )

        ax.set_aspect('equal', 'box')

    # --- Axis labels/limits, title

    if not rotate:
        ax.set_xlabel('x')
        ax.set_ylabel('y')

    ax.set_xlim(horizontal_range[0], horizontal_range[1])
    ax.set_ylim(vertical_range[0], vertical_range[1])

    if title is not None:
        ax.set_title(title)

    # TODO: make render_label respond to settings/options
    render_label = r'$\int$ ' + f'{render}' + ' dz'

    # --- Colorbar

    if render:

        cb = None
        if colorbar is not None:
            colorbar_ = colorbar
        else:
            colorbar_ = plot_options['colorbar']
        if colorbar_:
            divider = make_axes_locatable(ax)
            cax = divider.append_axes("right", size="5%", pad=0.05)
            cb = plt.colorbar(img, cax=cax)
            if render_label:
                cb.set_label(render_label)


def _convert_units():
    """Convert units."""

    # TODO: write this


def _calculate_quantity():
    """Calculate an extra quantity."""

    # TODO: write this


def set_colorbar(cb, vmin=None, vmax=None):
    """Set colorbar limits."""

    if vmin is None:
        vmin = cb._colorbar.vmin

    if vmax is None:
        vmax = cb._colorbar.vmax

    cb._image.set_clim([vmin, vmax])


def _interpolation_weights(density_weighted, particles, hfact):
    """Calculate interpolation weights."""

    if density_weighted:
        return np.array(particles['m'] / particles['h'] ** 2)

    return np.full_like(particles['h'], 1 / hfact)


def _rotate_frame(positions, velocities, axis, theta):
    """Rotate around axis."""

    return (
        rotate_vector_arbitrary_axis(positions, axis, theta),
        rotate_vector_arbitrary_axis(velocities, axis, theta),
    )
