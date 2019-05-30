"""
Plot options.
"""

from collections import namedtuple

FigureOptions = namedtuple(
    'FigureOptions',
    [
        'axis',
        'colorbar',
        'colormap',
        'figure',
        'font_family',
        'font_size',
        'title',
    ],
)

ImageRangeOptions = namedtuple(
    'ImageRangeOptions', ['xrange', 'yrange', 'extent']
)

InterpolationOptions = namedtuple(
    'InterpolationOptions',
    [
        'accelerate',
        'cross_section',
        'density_weighted',
        'normalize',
        'number_pixels',
        'observer_distance',
        'opacity',
        'perspective',
        'slice_position',
    ],
)

RenderOptions = namedtuple(
    'RenderOptions',
    ['render_scale', 'render_min', 'render_max', 'render_fraction_max'],
)

RotationOptions = namedtuple(
    'RotationOptions',
    ['rotation_axis', 'rotation_angle', 'position_angle', 'inclination'],
)

VectorOptions = namedtuple(
    'VectorOptions', ['stream', 'stride', 'vector_color']
)

PlotOptions = namedtuple(
    'PlotOptions',
    [
        'FigureOptions',
        'ImageRangeOptions',
        'InterpolationOptions',
        'RenderOptions',
        'RotationOptions',
        'VectorOptions',
    ],
)

DEFAULT_OPTIONS = PlotOptions(
    FigureOptions(
        axis=None,
        colorbar=True,
        colormap='gist_heat',
        figure=None,
        font_family='sans-serif',
        font_size=12,
        title=None,
    ),
    ImageRangeOptions(xrange=None, yrange=None, extent=None),
    InterpolationOptions(
        accelerate=False,
        cross_section=False,
        density_weighted=False,
        normalize=False,
        number_pixels=(512, 512),
        observer_distance=0.0,
        opacity=False,
        perspective=False,
        slice_position=0.0,
    ),
    RenderOptions(
        render_scale='linear',
        render_min=None,
        render_max=None,
        render_fraction_max=None,
    ),
    RotationOptions(
        rotation_axis=None,
        rotation_angle=None,
        position_angle=None,
        inclination=None,
    ),
    VectorOptions(stream=False, stride=25, vector_color='black'),
)


def plot_options(**kwargs):
    """
    Create a dictionary with plot options.

    Parameters
    ----------
    **kwargs
        Any values in PlotOptions namedtuple.

    Returns
    -------
    dict
        A dictionary of visualization options.
    """

    options = {
        key: opts._asdict() for key, opts in DEFAULT_OPTIONS._asdict().items()
    }
    for key, item in kwargs.items():
        if key in options:
            options[key] = item

    return options
