"""Analysis of SPH data.

The analysis sub-package contains Plonk implementations of typical
smoothed particle hydrodynamics post-simulation analysis tasks.

Examples
--------
Create a radial profile in the xy-plane.

>>> p = Profile(snap, cmin=10, cmax=200, n_bins=100)
>>> p.plot('radius', 'density')

Calculate the angular momentum on the particles.

>>> angmom = particles.angular_momentum(snap)

Calculate the total angular momentum over all particles.

>>> angmom_tot = total.angular_momentum(snap)

Calculate the Roche sphere radius given two sink particles.

>>> s1 = snap.sinks[0]
>>> s2 = snap.sinks[1]
>>> separation = plonk.utils.math.norm(s1['position'] - s2['position'])
>>> Roche = sinks.Roche_sphere(s1['mass'], s2['mass'], separation)
"""

from . import discs, filters, particles, sinks, sph, total
from .profile import Profile, load_profile

__all__ = [
    'Profile',
    'discs',
    'filters',
    'load_profile',
    'particles',
    'sinks',
    'sph',
    'total',
]
