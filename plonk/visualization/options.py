"""
Plot options.
"""

from dataclasses import dataclass, field
from typing import Any, List


@dataclass
class FigureOptions:
    axis: Any = None
    cbar_axis: Any = None
    colorbar: bool = True
    colormap: str = 'gist_heat'
    figure: Any = None
    font_family: str = 'sans-serif'
    font_size: int = 12
    render_label: str = None
    title: str = None


@dataclass
class InterpolationOptions:
    accelerate: bool = False
    cross_section: bool = False
    density_weighted: bool = False
    normalize: bool = False
    number_pixels: tuple = (512, 512)
    observer_distance: float = 0.0
    opacity: bool = False
    perspective: bool = False
    slice_position: float = 0.0


@dataclass
class RenderOptions:
    render_scale: str = 'linear'
    render_min: float = None
    render_max: float = None


@dataclass
class RotationOptions:
    rotation_axis: List[float] = None
    rotation_angle: float = None


@dataclass
class UnitsOptions:
    units: Any = None
    integrated_z: Any = None


@dataclass
class VectorOptions:
    stream: bool = False
    stride: int = 25
    vector_color: str = 'black'


def make_figure_options():
    return FigureOptions()


def make_interpolation_options():
    return InterpolationOptions()


def make_render_options():
    return RenderOptions()


def make_rotation_options():
    return RotationOptions()


def make_units_options():
    return UnitsOptions()


def make_vector_options():
    return VectorOptions()


@dataclass
class PlotOptions:

    figure: Any = field(default_factory=make_figure_options)
    interpolation: Any = field(default_factory=make_interpolation_options)
    render: Any = field(default_factory=make_render_options)
    rotation: Any = field(default_factory=make_rotation_options)
    units: Any = field(default_factory=make_units_options)
    vector: Any = field(default_factory=make_vector_options)
