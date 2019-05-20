"""
This is the visualization module.

It contains the Plonk implementation for visualizing smoothed particle
hydrodynamics simulations using Splash for interpolation. This is
accessed by the plot function.

Examples
--------
Rendering density.

>>> plonk.visualization.plot(dump, render='density')

Notes
-----
For the Splash source code see https://github.com/danieljprice/splash.
The user guide is at http://users.monash.edu.au/~dprice/splash. When
using any software derived from Splash for academic purposes you should
cite: Price, 2007, Publ. Astron. Soc. Aust., 24, 159-173.
"""

from .image import Visualization

visualization = Visualization()
plot = visualization.plot

__all__ = ['plot']
