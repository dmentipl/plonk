"""Radial profiles.

Heavily inspired by pynbody (https://pynbody.github.io/).
"""

from bisect import bisect
from typing import Dict, List, Optional, Tuple

import numpy as np
from numpy import ndarray

from ..core.dump import Dump


class Profile:
    """Radial profiles.

    Parameters
    ----------
    dump
        The Dump object.
    ndim, optional
        The dimension of the profile. For ndim == 2, the radial binning
        is cylindrical in the xy-plane. For ndim == 3, the radial
        binning is spherical. Default is 2.
    radius_min, optional
        The minimum radius for binning. Defaults to minimum on the
        particles.
    radius_max, optional
        The maximum radius for binning. Defaults to the 99 percentile
        distance.
    n_bins, optional
        The number of radial bins. Default is 100.
    ignore_accreted, optional
        Ignore particles accreted onto sinks. Default is True.
    mask, optional
        Select a subset of all particles via a NumPy mask array.
    """

    def __init__(
        self,
        dump: Dump,
        ndim: Optional[int] = 2,
        radius_min: Optional[float] = None,
        radius_max: Optional[float] = None,
        n_bins: int = 100,
        ignore_accreted: bool = True,
        mask: Optional[ndarray] = None,
    ):

        self.dump = dump
        self.ndim = ndim

        self._mask = self._setup_particle_mask(ignore_accreted, mask)
        self._x = self._calculate_x()
        self.range = self._set_range(radius_min, radius_max)
        self.n_bins = n_bins

        self.bin_edges, self.bin_sizes = self._setup_bins()
        self.bin_centers = 0.5 * (self.bin_edges[:-1] + self.bin_edges[1:])
        self._particle_bin = np.digitize(self._x, self.bin_edges)
        self.bin_indicies = self._set_particle_bin_indicies()

        self._profiles: Dict[str, ndarray] = {}
        self._profiles['radius'] = self.bin_centers
        n_per_bin, _ = np.histogram(self._x, self.bin_edges)
        self._profiles['number'] = n_per_bin

        self._profiles['mass'] = mass(self)
        self._profiles['density'] = density(self)
        self._profiles['angular momentum'] = angular_momentum(self)

    def _setup_particle_mask(
        self, ignore_accreted: bool, mask: Optional[ndarray]
    ) -> ndarray:
        if ignore_accreted is False:
            if mask is None:
                return np.ones(self.dump.particles.number, dtype=bool)
            return mask
        if mask is None:
            return self.dump.particles.arrays['h'][:] > 0.0
        return mask & self.dump.particles.arrays['h'][:] > 0.0

    def _calculate_x(self) -> ndarray:
        pos = self.dump.particles.arrays['xyz'][:]
        pos = pos[self._mask]
        if self.ndim == 2:
            return np.hypot(pos[:, 0], pos[:, 1])
        elif self.ndim == 3:
            return np.hypot(np.hypot(pos[:, 0], pos[:, 1]), pos[:, 2])

    def _set_range(
        self, radius_min: Optional[float], radius_max: Optional[float]
    ) -> Tuple[float, float]:
        rmin, rmax = (radius_min, radius_max)
        if radius_min is None:
            rmin = self._x.min()
        if radius_max is None:
            rmax = np.percentile(self._x, 99, axis=0)

        return rmin, rmax

    def _setup_bins(self) -> ndarray:
        bin_edges = np.linspace(self.range[0], self.range[1], self.n_bins + 1)
        if self.ndim == 2:
            bin_sizes = np.pi * (bin_edges[1:] ** 2 - bin_edges[:-1] ** 2)
        elif self.ndim == 3:
            bin_sizes = 4 / 3 * np.pi * (bin_edges[1:] ** 3 - bin_edges[:-1] ** 3)
        return bin_edges, bin_sizes

    def _set_particle_bin_indicies(self) -> List[ndarray]:
        sortind = self._particle_bin.argsort()
        sort_pind = self._particle_bin[sortind]
        binind = list()
        prev_index = bisect(sort_pind, 0)
        for i in range(self.n_bins):
            new_index = bisect(sort_pind, i + 1)
            binind.append(np.sort(sortind[prev_index:new_index]))
            prev_index = new_index
        return binind

    def __getitem__(self, name: str) -> ndarray:
        """Return the profile of a given kind."""
        if name in self._profiles:
            return self._profiles[name]
        else:
            raise ValueError('Profile not available')

    def __delitem__(self, name):
        """Delete a profile from memory."""
        del self._profiles[name]

    def __repr__(self):
        """Object repr method."""
        return (
            f'<plonk.Profile: {self.n_bins} bins; '
            f'{self.range[0]:.3g} to {self.range[1]:.3g}>'
        )

    def keys(self):
        """Return a listing of available profile types."""
        return self._profiles.keys()


def mass(profile: Profile) -> ndarray:
    """Mass profile."""
    mass = np.zeros(profile.n_bins)
    M = profile.dump.mass[profile._mask]
    for idx, bin_ind in enumerate(profile.bin_indicies):
        mass[idx] = M[bin_ind].sum()
    return mass


def density(profile: Profile) -> ndarray:
    """Density profile.

    Units are [mass / length ** ndim], which depends on ndim of profile.
    """
    return profile._profiles['mass'] / profile.bin_sizes


def angular_momentum(profile: Profile) -> ndarray:
    """Angular momentum profile.

    The angular momentum is with respect to the xy-plane.
    """
    mass = profile.dump.mass[:, np.newaxis]
    pos = profile.dump.particles.arrays['xyz'][:]
    vel = profile.dump.particles.arrays['vxyz'][:]
    L = mass * np.cross(pos, vel)
    L_z = L[:, 2]
    angular_momentum = np.zeros(profile.n_bins)
    for idx, bin_ind in enumerate(profile.bin_indicies):
        angular_momentum[idx] = L_z[bin_ind].mean()
    return angular_momentum
