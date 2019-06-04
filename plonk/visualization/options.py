"""
Plot options.
"""

FIGURE_OPTIONS = {
    'axis': None,
    'cbar_axis': None,
    'colorbar': True,
    'colormap': 'gist_heat',
    'figure': None,
    'font_family': 'sans-serif',
    'font_size': 12,
    'render_label': None,
    'title': None,
}

INTERPOLATION_OPTIONS = {
    'accelerate': False,
    'cross_section': False,
    'density_weighted': False,
    'normalize': False,
    'number_pixels': (512, 512),
    'observer_distance': 0.0,
    'opacity': False,
    'perspective': False,
    'slice_position': 0.0,
}

RENDER_OPTIONS = {
    'render_scale': 'linear',
    'render_min': None,
    'render_max': None,
}

ROTATION_OPTIONS = {'rotation_axis': None, 'rotation_angle': None}

UNITS_OPTIONS = {'units': None, 'integrated_z': None}

VECTOR_OPTIONS = {'stream': False, 'stride': 25, 'vector_color': 'black'}

DEFAULT_OPTIONS = {
    'figure': FIGURE_OPTIONS,
    'interpolation': INTERPOLATION_OPTIONS,
    'render': RENDER_OPTIONS,
    'rotation': ROTATION_OPTIONS,
    'units': UNITS_OPTIONS,
    'vector': VECTOR_OPTIONS,
}
