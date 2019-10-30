"""The visualization sub-package.

The Plonk implementation for visualizing smoothed particle hydrodynamics
simulations using kernel density estimation based interpolation.

Examples
--------
Rendering the surface density.

>>> viz = plonk.Visualization(
...     scalar_data=density,
...     x_coordinate=x_position,
...     y_coordinate=y_position,
...     extent=extent,
...     particle_mass=particle_mass,
...     smoothing_length=smoothing_length,
... )
"""

from .visualization import Visualization

__all__ = ['Visualization']
