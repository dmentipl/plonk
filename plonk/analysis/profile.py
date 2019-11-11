"""Radial profiles.

Heavily inspired by pynbody (https://pynbody.github.io/).
"""

from bisect import bisect
from typing import Any, Callable, Collection, Dict, List, Optional, Tuple, Union

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from numpy import ndarray
from pandas import DataFrame

from ..snap import Snap


class Profile:
    """Radial profiles.

    Parameters
    ----------
    snap
        The Snap object.
    ndim
        The dimension of the profile. For ndim == 2, the radial binning
        is cylindrical in the xy-plane. For ndim == 3, the radial
        binning is spherical. Default is 2.
    radius_min
        The minimum radius for binning. Defaults to minimum on the
        particles.
    radius_max
        The maximum radius for binning. Defaults to the 99 percentile
        distance.
    n_bins
        The number of radial bins. Default is 100.
    ignore_accreted
        Ignore particles accreted onto sinks. Default is True.
    mask
        Select a subset of all particles via a NumPy mask array.
    """

    _profile_functions: Dict[str, Callable] = {}

    def __init__(
        self,
        snap: Snap,
        ndim: Optional[int] = 2,
        radius_min: Optional[float] = None,
        radius_max: Optional[float] = None,
        n_bins: int = 100,
        ignore_accreted: bool = True,
        mask: Optional[ndarray] = None,
    ):

        self.snap = snap
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
        self._profiles['number'] = np.histogram(self._x, self.bin_edges)[0]

    def _setup_particle_mask(
        self, ignore_accreted: bool, mask: Optional[ndarray]
    ) -> ndarray:
        if ignore_accreted is False:
            if mask is None:
                return np.ones(len(self.snap), dtype=bool)
            return mask
        h: ndarray = self.snap['h']
        if mask is None:
            return h > 0
        return mask & h > 0

    def _calculate_x(self) -> ndarray:
        pos = self.snap['xyz']
        pos = pos[self._mask]
        if self.ndim == 2:
            return np.hypot(pos[:, 0], pos[:, 1])
        elif self.ndim == 3:
            return np.hypot(np.hypot(pos[:, 0], pos[:, 1]), pos[:, 2])

    def _set_range(
        self, radius_min: Optional[float], radius_max: Optional[float]
    ) -> Tuple[float, float]:
        if radius_min is None:
            rmin = self._x.min()
        else:
            rmin = radius_min
        if radius_max is None:
            rmax = np.percentile(self._x, 99, axis=0)
        else:
            rmax = radius_max

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

    def _get_profile(self, name: str, args: Optional[Tuple[Any, ...]] = None):
        """Get a profile by name."""
        if name in self._profiles:
            return self._profiles[name]

        elif name in Profile._profile_functions:
            if args is not None:
                self._profiles[name] = Profile._profile_functions[name](self, *args)
            self._profiles[name] = Profile._profile_functions[name](self)
            return self._profiles[name]

        else:
            raise ValueError('Profile not available')

    def __getitem__(self, name: str) -> ndarray:
        """Return the profile of a given kind."""
        if isinstance(name, tuple):
            name, *args = name
            return self._get_profile(name, args)
        return self._get_profile(name)

    def __setitem__(self, name: str, values: ndarray):
        """Set the profile directly."""
        if name in self._profiles:
            self._profiles[name] = values
        else:
            raise KeyError(f'{name} is not a valid profile')

    def __delitem__(self, name):
        """Delete a profile from memory."""
        del self._profiles[name]

    def __repr__(self):
        """Object repr method."""
        return f'<plonk.Profile: {self.n_bins} bins>'

    def loaded_keys(self):
        """Return a listing of loaded profiles."""
        return tuple(sorted(self._profiles.keys()))

    def available_keys(self):
        """Return a listing of available profiles."""
        loaded = list(self.loaded_keys())
        available = list(self._profile_functions.keys())
        return tuple(sorted(set(loaded + available)))

    @staticmethod
    def profile_property(fn):
        """Decorate profile functions."""
        Profile._profile_functions[fn.__name__] = fn
        return fn

    def plot(self, x: str, y: Union[str, Collection[str]]):
        """Plot profile.

        Parameters
        ----------
        x
            The x axis to plot as a string.
        y
            The y axis to plot. Can be multiple as a list or tuple.
        """
        if x.lower() not in self.available_keys():
            raise ValueError('Cannot determine x axis to plot')
        _x = self._get_profile(x)
        if isinstance(y, (list, tuple)):
            for yi in y:
                if yi.lower() not in self.available_keys():
                    raise ValueError('Cannot determine y axis to plot')
                _y = self._get_profile(yi)
                plt.plot(_x, _y)
        elif isinstance(y, str):
            if y.lower() not in self.available_keys():
                raise ValueError('Cannot determine y axis to plot')
            _y = self._get_profile(y)
            plt.plot(_x, _y)
        else:
            raise ValueError('Cannot determine y axis to plot')
        return plt.gcf(), plt.gca()

    def to_dataframe(self, all_available: bool = False) -> DataFrame:
        """Convert Profile to DataFrame.

        Parameters
        ----------
        all_available
            If True, this will calculate all available profiles before
            converting to a DataFrame.
        """
        if all_available:
            columns = self.available_keys()
        else:
            columns = self.loaded_keys()
        data = dict()
        for column in columns:
            data[column] = self[column]
        return pd.DataFrame(data)


@Profile.profile_property
def mass(self) -> ndarray:
    """Mass profile."""
    mass = np.zeros(self.n_bins)
    M = self.snap['mass'][self._mask]
    for idx, bin_ind in enumerate(self.bin_indicies):
        mass[idx] = M[bin_ind].sum()
    return mass


@Profile.profile_property
def density(self) -> ndarray:
    """Density profile.

    Units are [mass / length ** ndim], which depends on ndim of profile.
    """
    return self._get_profile('mass') / self.bin_sizes


@Profile.profile_property
def smooth(self) -> ndarray:
    """Smoothing length profile."""
    smooth = np.zeros(self.n_bins)
    _h = self.snap['h']
    h = _h[self._mask]
    for idx, bin_ind in enumerate(self.bin_indicies):
        smooth[idx] = h[bin_ind].mean()
    return smooth


@Profile.profile_property
def scale_height(self) -> ndarray:
    """Scale height profile."""
    scale_height = np.zeros(self.n_bins)
    pos = self.snap['xyz']
    z = pos[self._mask, 2]
    for idx, bin_ind in enumerate(self.bin_indicies):
        scale_height[idx] = z[bin_ind].std()
    return scale_height


@Profile.profile_property
def angmom_mag(self) -> ndarray:
    """Magnitude of specific angular momentum profile."""
    angmom_mag = np.zeros(self.n_bins)

    mass = self.snap['mass'][:, np.newaxis]
    pos = self.snap['xyz']
    vel = self.snap['vxyz']

    J = mass * np.cross(pos, vel)
    J_mag = np.hypot(np.hypot(J[:, 0], J[:, 1]), J[:, 2])

    for idx, bin_ind in enumerate(self.bin_indicies):
        angmom_mag[idx] = J_mag[bin_ind].mean() / self._get_profile('mass')[idx]
    return angmom_mag


@Profile.profile_property
def angmom_theta(self) -> ndarray:
    """Angle between angular momentum and xy-plane."""
    angmom_theta = np.zeros(self.n_bins)

    mass = self.snap['mass'][:, np.newaxis]
    pos = self.snap['xyz']
    vel = self.snap['vxyz']

    J = mass * np.cross(pos, vel)
    J_z = J[:, 2]
    J_mag = np.hypot(np.hypot(J[:, 0], J[:, 1]), J[:, 2])

    for idx, bin_ind in enumerate(self.bin_indicies):
        j_z = J_z[bin_ind].mean() / self._get_profile('mass')[idx]
        j_mag = J_mag[bin_ind].mean() / self._get_profile('mass')[idx]
        angmom_theta[idx] = np.arccos(j_z / j_mag)
    return angmom_theta


@Profile.profile_property
def angmom_phi(self) -> ndarray:
    """Angle between angular momentum and x-axis in xy-plane."""
    angmom_phi = np.zeros(self.n_bins)

    mass = self.snap['mass'][:, np.newaxis]
    pos = self.snap['xyz']
    vel = self.snap['vxyz']

    J = mass * np.cross(pos, vel)
    J_x, J_y = (J[:, 0], J[:, 1])

    for idx, bin_ind in enumerate(self.bin_indicies):
        j_x = J_x[bin_ind].mean() / self._get_profile('mass')[idx]
        j_y = J_y[bin_ind].mean() / self._get_profile('mass')[idx]
        angmom_phi[idx] = np.arctan2(j_y, j_x)
    return angmom_phi


@Profile.profile_property
def eccentricity(self, gravitational_parameter: float):
    """Orbital eccentricity profile."""
    mu = gravitational_parameter

    pos = self.snap['xyz'][self._mask]
    vel = self.Snap['vxyz'][self._mask]

    r = np.hypot(np.hypot(pos[:, 0], pos[:, 1]), pos[:, 2])
    v = np.hypot(np.hypot(vel[:, 0], vel[:, 1]), vel[:, 2])

    h = np.cross(pos, vel)
    h_mag = np.hypot(np.hypot(h[:, 0], h[:, 1]), h[:, 2])

    ke = 0.5 * v ** 2
    pe = -mu / r
    e = ke + pe
    term = 2 * e * h_mag ** 2 / mu ** 2
    ecc = np.sqrt(1 + term)

    eccentricity = np.zeros(self.n_bins)
    for idx, bin_ind in enumerate(self.bin_indicies):
        eccentricity[idx] = ecc[bin_ind].mean()
    return eccentricity
