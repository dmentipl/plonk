"""
Plot options.
"""

from collections import namedtuple

FigureOptions = namedtuple(
    'FigureOptions',
    [
        'axis',
        'cbar_axis',
        'colorbar',
        'colormap',
        'figure',
        'font_family',
        'font_size',
        'render_label',
        'title',
    ],
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
    'RenderOptions', ['render_scale', 'render_min', 'render_max']
)

RotationOptions = namedtuple(
    'RotationOptions', ['rotation_axis', 'rotation_angle']
)

UnitsOptions = namedtuple('UnitsOptions', ['units', 'integrated_z'])

VectorOptions = namedtuple(
    'VectorOptions', ['stream', 'stride', 'vector_color']
)

PlotOptions = namedtuple(
    'PlotOptions',
    [
        'FigureOptions',
        'InterpolationOptions',
        'RenderOptions',
        'RotationOptions',
        'UnitsOptions',
        'VectorOptions',
    ],
)

DEFAULT_OPTIONS = PlotOptions(
    FigureOptions(
        axis=None,
        cbar_axis=None,
        colorbar=True,
        colormap='gist_heat',
        figure=None,
        font_family='sans-serif',
        font_size=12,
        render_label=None,
        title=None,
    ),
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
    RenderOptions(render_scale='linear', render_min=None, render_max=None),
    RotationOptions(
        rotation_axis=None,
        rotation_angle=None,
    ),
    UnitsOptions(units=None, integrated_z=None),
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
