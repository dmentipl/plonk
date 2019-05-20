"""
This is the visualization module.

It contains the Plonk implementation for rendering smoothed particle
hydrodynamics simulations using Splash. This is accessed by the plot
function.

Examples
--------
Rendering density.

>>> plonk.visualization.plot(dump, render='density')
"""

from .image import Visualization

visualization = Visualization()
plot = visualization.plot

__all__ = ['plot']
