"""The visualization sub-package.

The Plonk implementation for visualizing smoothed particle hydrodynamics
simulations using kernel density estimation based interpolation.

Examples
--------
Plot the particles.

>>> viz = plonk.visualize.plot(
...     snap=snap,
...     x='x',
...     y='y',
...     extent=(-100, 100, -100, 100),
... )

Render the surface density in xz-plane.

>>> viz = plonk.visualize.plot(
...     snap=snap,
...     quantity='density',
...     x='x',
...     y='z',
...     extent=(-100, 100, -25, 25),
... )

Get the interpolation to grid directly (without plotting).

>>> grid_data = plonk.visualize.interpolate(
...     snap=snap,
...     quantity='density',
...     interp='projection',
...     extent=(-100, 100, -100, 100),
... )

Make an animation of multiple snaps.

>>> plonk.visualize.animation(
...     snaps=snaps,
...     quantity='density',
...     extent=(-100, 100, -100, 100),
...     filename='animation.mp4',
... )

Set units for the plot.

>>> units = {
...     'quantity': plonk.units('g / cm ** 3'),
...     'extent': plonk.units('au'),
...     'projection': plonk.units('cm'),
... }

>>> viz = plonk.visualize.plot(
...     snap=snap,
...     quantity='density',
...     extent=(-100, 100, -100, 100),
...     units=units,
... )
"""

from .animation import animation
from .interpolation import interpolate
from .visualization import Visualization, plot

__all__ = ['Visualization', 'animation', 'interpolate', 'plot']
