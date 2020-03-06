"""The visualization sub-package.

The Plonk implementation for visualizing smoothed particle hydrodynamics
simulations using kernel density estimation based interpolation.

Examples
--------
Rendering the surface density.

>>> viz = plonk.visualize.plot(
...     scalar_data=density,
...     x_coordinate=x_position,
...     y_coordinate=y_position,
...     extent=extent,
...     particle_mass=particle_mass,
...     smoothing_length=smoothing_length,
... )

Or via the helper 'render' function.

>>> viz = plonk.visualize.render(snap=snap, quantity='density')
>>> viz = plonk.visualize.render(snap=snap, quantity=snap['density'])

Get the interpolation to grid directly (without plotting).

>>> grid_data = plonk.visualize.interpolate(
...     snap=snap,
...     quantity='density',
...     extent=(-100, 100, -100, 100),
... )

Make an animation of multiple snaps.

>>> plonk.visualize.animation(
...     snaps=snaps,
...     quantity='density',
...     extent=(-100, 100, -100, 100),
...     filename='animation.mp4',
... )
"""

from .animation import animation
from .plot import interpolate, plot, render
from .visualization import Visualization

__all__ = ['Visualization', 'animation', 'interpolate', 'plot', 'render']
