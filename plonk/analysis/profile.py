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
from .particles import eccentricity as _eccentricity
from .particles import specific_angular_momentum


class Profile:
    """Radial profiles.

    A radial profile is a binning of particles in cylindrical or
    spherical shells around the origin, i.e. (0, 0, 0). For cylindrical
    profiles the cylindrical cells are perpendicular to the xy-plane,
    i.e. the bins are averaged azimuthally and in the z-direction.

    Parameters
    ----------
    snap
        The Snap object.
    ndim : optional
        The dimension of the profile. For ndim == 2, the radial binning
        is cylindrical in the xy-plane. For ndim == 3, the radial
        binning is spherical. Default is 2.
    radius_min : optional
        The minimum radius for binning. Defaults to minimum on the
        particles.
    radius_max : optional
        The maximum radius for binning. Defaults to the 99 percentile
        distance.
    n_bins : optional
        The number of radial bins. Default is 100.
    spacing : optional
        The spacing of radial bins. Can be 'linear' or 'log'. Default is
        'linear'.
    ignore_accreted : optional
        Ignore particles accreted onto sinks. Default is True.

    Examples
    --------
    Generate profile from snapshot.

    >>> prof = plonk.Profile(snap=snap)
    >>> prof = plonk.Profile(snap=snap, n_bins=300)
    >>> prof = plonk.Profile(snap=snap, radius_min=1, radius_max=300)

    To access a profile.

    >>> prof['density']
    >>> prof['scale_height']

    To set a new profile.

    >>> prof['aspect_ratio'] = prof['scale_height'] / prof['radius']

    Alternatively use the profile_property decorator.

    >>> @Profile.profile_property
    ... def mass(prof):
    ...     M = prof.snap['mass']
    ...     return prof.particles_to_binned_quantity(np.sum, M)
    """

    _profile_functions: Dict[str, Callable] = {}

    @staticmethod
    def profile_property(fn: Callable) -> Callable:
        """Decorate function to add profile to Profile.

        Parameters
        ----------
        fn
            A function that returns the profile as an array. The name of
            the function is the string with which to reference the array.

        Returns
        -------
        Callable
            The function which returns the array.
        """
        Profile._profile_functions[fn.__name__] = fn
        return fn

    def __init__(
        self,
        snap: Snap,
        ndim: Optional[int] = 2,
        radius_min: Optional[float] = None,
        radius_max: Optional[float] = None,
        n_bins: int = 100,
        spacing: str = 'linear',
        ignore_accreted: bool = True,
    ):

        self.snap = snap
        self.ndim = ndim
        self.spacing = self._check_spacing(spacing)
        self.properties: Dict[str, Any] = {}

        self._profiles: Dict[str, ndarray] = {}

        self._mask = self._setup_particle_mask(ignore_accreted)
        self._x = self._calculate_x()
        self.range = self._set_range(radius_min, radius_max)
        self.n_bins = n_bins

        self.bin_edges, self['size'] = self._setup_bins()
        self.bin_centers = 0.5 * (self.bin_edges[:-1] + self.bin_edges[1:])
        self._particle_bin = np.digitize(self._x, self.bin_edges)
        self.bin_indicies = self._set_particle_bin_indicies()

        self._profiles['radius'] = self.bin_centers
        self._profiles['number'] = np.histogram(self._x, self.bin_edges)[0]

    def _check_spacing(self, spacing):
        if spacing.lower() in ('lin', 'linear'):
            spacing = 'linear'
        elif spacing.lower() in ('log', 'logarithm', 'logarithmic'):
            spacing = 'log'
        else:
            raise ValueError('Cannot determine spacing')
        return spacing

    def _setup_particle_mask(self, ignore_accreted: bool) -> ndarray:
        if ignore_accreted is False:
            return np.ones(len(self.snap), dtype=bool)
        h: ndarray = self.snap['h']
        return h > 0

    def _calculate_x(self) -> ndarray:
        pos: ndarray = self.snap['xyz']
        pos = pos[self._mask]
        if self.ndim == 2:
            return np.sqrt(pos[:, 0] ** 2 + pos[:, 1] ** 2)
        elif self.ndim == 3:
            return np.sqrt(pos[:, 0] ** 2 + pos[:, 1] ** 2 + pos[:, 2] ** 2)

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
        if self.spacing == 'linear':
            bin_edges = np.linspace(self.range[0], self.range[1], self.n_bins + 1)
        elif self.spacing == 'log':
            bin_edges = np.logspace(
                np.log10(self.range[0]), np.log10(self.range[1]), self.n_bins + 1
            )
        else:
            raise ValueError('Cannot determine spacing to setup bins')
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
            else:
                self._profiles[name] = Profile._profile_functions[name](self)
            return self._profiles[name]

        else:
            raise ValueError('Profile not available')

    def __getitem__(self, name: str) -> ndarray:
        """Return the profile of a given kind."""
        return self._get_profile(name)

    def __setitem__(self, name: str, item: ndarray):
        """Set the profile directly."""
        if not isinstance(item, ndarray):
            raise ValueError('"item" must be ndarray')
        if item.shape[0] != self.n_bins:
            raise ValueError('Length of array does not match number of bins')
        if name in self.loaded_keys():
            raise ValueError(
                'Attempting to overwrite existing profile. To do so, first delete the '
                'profile\nwith del prof["profile"], then try again.'
            )
        elif name in self.available_keys():
            raise ValueError(
                'Attempting to set profile already available. '
                'See prof.available_keys().'
            )
        self._profiles[name] = item

    def __delitem__(self, name):
        """Delete a profile from memory."""
        del self._profiles[name]

    def __len__(self):
        """Length as number of bins."""
        return self.n_bins

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

    def particles_to_binned_quantity(self, function, *args) -> ndarray:
        """Calculate binned quantities from particles.

        This takes care of the bin indices and ignoring accreted
        particles (if requested in instantiating the profile).

        Parameters
        ----------
        function
            The function that acts on particles in the radial bin. The
            function should produce a single number, i.e. it is a
            reduction. E.g. to calculate the mean of a particle quantity
            in a radial bin use np.mean.
        *args
            The arguments to the function. Typically this is one or more
            particle arrays or other quantities on a snapshot.
        """
        binned_quantity = np.zeros(self.n_bins)
        for idx, bin_ind in enumerate(self.bin_indicies):
            binned_quantity[idx] = function(
                tuple([arg[self._mask][bin_ind] for arg in args])
            )
        return binned_quantity


@Profile.profile_property
def mass(prof) -> ndarray:
    """Mass profile."""
    M = prof.snap['mass']
    return prof.particles_to_binned_quantity(np.sum, M)


@Profile.profile_property
def density(prof) -> ndarray:
    """Density profile.

    Units are [mass / length ** ndim], which depends on ndim of profile.
    """
    return prof['mass'] / prof['size']


@Profile.profile_property
def smooth(prof) -> ndarray:
    """Smoothing length profile."""
    h = prof.snap['smooth']
    return prof.particles_to_binned_quantity(np.mean, h)


@Profile.profile_property
def scale_height(prof) -> ndarray:
    """Scale height profile."""
    z = prof.snap['z']
    return prof.particles_to_binned_quantity(np.std, z)


@Profile.profile_property
def aspect_ratio(prof) -> ndarray:
    """Aspect ratio profile."""
    H = prof['scale_height']
    R = prof['radius']
    return H / R


@Profile.profile_property
def angmom_mag(prof) -> ndarray:
    """Magnitude of specific angular momentum profile."""
    angmom = specific_angular_momentum(snap=prof.snap)
    angmom_mag = np.linalg.norm(angmom, axis=1)

    return prof.particles_to_binned_quantity(np.mean, angmom_mag)


@Profile.profile_property
def angmom_x(prof) -> ndarray:
    """x-component of specific angular momentum profile."""
    angmom = specific_angular_momentum(snap=prof.snap)
    angmom_x = angmom[:, 0]

    return prof.particles_to_binned_quantity(np.mean, angmom_x)


@Profile.profile_property
def angmom_y(prof) -> ndarray:
    """y-component of specific angular momentum profile."""
    angmom = specific_angular_momentum(snap=prof.snap)
    angmom_y = angmom[:, 1]

    return prof.particles_to_binned_quantity(np.mean, angmom_y)


@Profile.profile_property
def angmom_z(prof) -> ndarray:
    """z-component of specific angular momentum profile."""
    angmom = specific_angular_momentum(snap=prof.snap)
    angmom_z = angmom[:, 2]

    return prof.particles_to_binned_quantity(np.mean, angmom_z)


@Profile.profile_property
def angmom_theta(prof) -> ndarray:
    """Angle between specific angular momentum and xy-plane."""
    angmom_z = prof['angmom_z']
    angmom_mag = prof['angmom_mag']

    return np.arccos(angmom_z / angmom_mag)


@Profile.profile_property
def angmom_phi(prof) -> ndarray:
    """Angle between specific angular momentum and x-axis in xy-plane."""
    angmom_x = prof['angmom_x']
    angmom_y = prof['angmom_y']
    return np.arctan2(angmom_y, angmom_x)


@Profile.profile_property
def eccentricity(prof) -> ndarray:
    """Orbital eccentricity profile.

    This profile assumes the central mass is at (0, 0, 0) and that the
    gravitational parameter (G*M) is set.
    """
    try:
        gravitational_parameter = prof.properties['gravitational_parameter']
    except KeyError:
        raise ValueError(
            'To generate eccentricity profile, first set the gravitational parameter\n'
            'via prof.properties["gravitational_parameter"].'
        )
    ecc = _eccentricity(
        snap=prof.snap, gravitational_parameter=gravitational_parameter,
    )
    return prof.particles_to_binned_quantity(np.mean, ecc)
