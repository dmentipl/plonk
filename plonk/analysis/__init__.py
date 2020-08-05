"""The analysis sub-package.

It contains Plonk implementations of typical smoothed particle
hydrodynamics post-simulation analysis tasks.

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

>>> m1 = snap.sinks['mass'][0]
>>> m2 = snap.sinks['mass'][1]
>>> separation = np.linalg.norm(
        snap.sinks['position'][0] - snap.sinks['position'][1]
    )
>>> Roche = sinks.Roche_sphere(m1, m2, separation)
"""

from . import particles, sinks, total
from .profile import Profile, load_profile

__all__ = ['Profile', 'load_profile', 'particles', 'sinks', 'total']
