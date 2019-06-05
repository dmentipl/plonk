"""
Plot options.
"""

from dataclasses import dataclass
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


@dataclass
class PlotOptions:
    figure: Any
    interpolation: Any
    render: Any
    rotation: Any
    units: Any
    vector: Any


DEFAULT_OPTIONS = PlotOptions(
    FigureOptions(),
    InterpolationOptions(),
    RenderOptions(),
    RotationOptions(),
    UnitsOptions(),
    VectorOptions(),
)
