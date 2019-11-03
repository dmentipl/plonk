"""The analysis sub-package.

It contains Plonk implementations of typical smoothed particle
hydrodynamics post-simulation analysis tasks.

Examples
--------
Create a radial profile in the xy-plane.

>>> p = Profile(snap, radius_min=10, radius_max=200, n_bins=100)
>>> p.plot('radius', 'density')
"""

from .profile import Profile

__all__ = ['Profile']
