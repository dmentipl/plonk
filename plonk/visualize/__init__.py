"""The visualization sub-package.

The Plonk implementation for visualizing smoothed particle hydrodynamics
simulations using kernel density estimation based interpolation.

Examples
--------
Plot the particles.

>>> viz = plonk.visualize.plot(
...     snap=snap,
...     x_coordinate=snap['x'],
...     y_coordinate=snap['y'],
...     extent=(-100, 100, -100, 100),
... )

Render the surface density in xz-plane.

>>> viz = plonk.visualize.plot(
...     snap=snap,
...     scalar_data=snap['density'],
...     x_coordinate='x',
...     y_coordinate='z',
...     extent=(-100, 100, -25, 25),
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
