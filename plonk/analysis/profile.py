"""Radial profiles.

Heavily inspired by pynbody (https://pynbody.github.io/).
"""

from bisect import bisect
from typing import Any, Callable, Dict, List, Optional, Tuple, Union

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from numpy import ndarray
from pandas import DataFrame

from .. import Quantity
from .. import units as plonk_units
from ..snap import Snap, gravitational_constant_in_code_units
from ..utils.math import norm


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
    >>> prof = plonk.Profile(snap=snap, spacing='log')

    If snap has physical units and setting 'radius_min' or 'radius_max'
    must use physical quantities.

    >>> radius_min = plonk.Quantity('1 au')
    >>> radius_max = plonk.Quantity('300 au')
    >>> prof = plonk.Profile(
    ...     snap=snap,
    ...     radius_min=radius_min,
    ...     radius_max=radius_max
    ... )

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

    Plot one or many quantities on the profile.

    >>> prof.plot('radius', 'density')
    >>> prof.plot('radius', ['angular_momentum_x', 'angular_momentum_y'])

    Plot a quantity on the profile with units.

    >>> prof.plot('radius', 'density', x_unit='au', y_unit='g/cm**2')
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
        radius_min: Optional[Any] = None,
        radius_max: Optional[Any] = None,
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
        if self.snap._physical_units:
            self._particle_bin = np.digitize(
                self._x.magnitude, self.bin_edges.magnitude
            )
        else:
            self._particle_bin = np.digitize(self._x, self.bin_edges)
        self.bin_indicies = self._set_particle_bin_indicies()

        self._profiles['radius'] = self.bin_centers
        if self.snap._physical_units:
            self._profiles['number'] = np.histogram(
                self._x.magnitude, self.bin_edges.magnitude
            )[0]
        else:
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
        self, radius_min: Optional[Any], radius_max: Optional[Any]
    ) -> Tuple[float, float]:
        if self.snap._physical_units:
            return self._set_range_physical_units(radius_min, radius_max)
        else:
            return self._set_range_code_units(radius_min, radius_max)

    def _set_range_code_units(
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

    def _set_range_physical_units(
        self, radius_min: Optional[Any], radius_max: Optional[Any]
    ) -> Tuple[float, float]:
        if radius_min is None:
            rmin = self._x.min()
        else:
            if not isinstance(radius_min, Quantity):
                raise ValueError(
                    'Snap has physical units: must use dimensional radius_min'
                )
            rmin = radius_min.to_base_units()
        if radius_max is None:
            rmax = np.percentile(self._x.magnitude, 99, axis=0) * self._x.units
        else:
            if not isinstance(radius_max, Quantity):
                raise ValueError(
                    'Snap has physical units: must use dimensional radius_max'
                )
            rmax = radius_max.to_base_units()

        return rmin, rmax

    def _setup_bins(self) -> ndarray:
        if self.snap._physical_units:
            bin_edges = self._bin_edges_physical_units()
        else:
            bin_edges = self._bin_edges_code_units()
        if self.ndim == 2:
            bin_sizes = np.pi * (bin_edges[1:] ** 2 - bin_edges[:-1] ** 2)
        elif self.ndim == 3:
            bin_sizes = 4 / 3 * np.pi * (bin_edges[1:] ** 3 - bin_edges[:-1] ** 3)
        return bin_edges, bin_sizes

    def _bin_edges_physical_units(self):
        if self.spacing == 'linear':
            bin_edges = (
                np.linspace(
                    self.range[0].magnitude, self.range[1].magnitude, self.n_bins + 1
                )
                * self.range[0].units
            )
        elif self.spacing == 'log':
            bin_edges = (
                np.logspace(
                    np.log10(self.range[0].magnitude),
                    np.log10(self.range[1].magnitude),
                    self.n_bins + 1,
                )
                * self.range[0].units
            )
        else:
            raise ValueError('Cannot determine spacing to setup bins')
        return bin_edges

    def _bin_edges_code_units(self):
        if self.spacing == 'linear':
            bin_edges = np.linspace(self.range[0], self.range[1], self.n_bins + 1)
        elif self.spacing == 'log':
            bin_edges = np.logspace(
                np.log10(self.range[0]), np.log10(self.range[1]), self.n_bins + 1
            )
        else:
            raise ValueError('Cannot determine spacing to setup bins')
        return bin_edges

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

        elif name in Profile._profile_functions:
            self._profiles[name] = Profile._profile_functions[name](self)
            return self._profiles[name]

        elif name in self.snap.available_arrays():
            array: ndarray = self.snap[name]
            if array.ndim == 1:
                self._profiles[name] = self.particles_to_binned_quantity(np.mean, array)
                if name not in self.snap.loaded_arrays():
                    del self.snap[name]
                return self._profiles[name]
            else:
                if name.split('_')[0] == 'dust' or name in ('stopping_time',):
                    raise ValueError(
                        'To access dust profiles try, for example, '
                        'prof["stopping_time_001"] or\nprof["dust_mass_total"]'
                    )
                else:
                    raise ValueError('Cannot determine profile to return')

        elif '_'.join(name.split('_')[:-1]) in self.snap.available_arrays():
            name_root = '_'.join(name.split('_')[:-1])
            name_suffix = name.split('_')[-1]
            if name_suffix == 'x':
                array = self.snap[name_root][:, 0]
            elif name_suffix == 'y':
                array = self.snap[name_root][:, 1]
            elif name_suffix == 'z':
                array = self.snap[name_root][:, 2]
            elif name_suffix == 'magnitude':
                array = norm(self.snap[name_root], axis=1)
            elif name_suffix == 'total':
                array = np.sum(self.snap[name_root], axis=1)
            elif _str_is_int(name_suffix):
                array = self.snap[name_root][:, int(name_suffix) - 1]
                name = name_root + f'_{int(name_suffix):03}'
            else:
                raise ValueError('Cannot determine profile')
            if array.ndim == 1:
                self._profiles[name] = self.particles_to_binned_quantity(np.mean, array)
                if name_root not in self.snap.loaded_arrays():
                    del self.snap[name_root]
                return self._profiles[name]
            else:
                raise ValueError(
                    'Cannot take profile of vector quantity. Try, for example,\n'
                    'prof["velocity_x"] or prof["velocity_magnitude"], etc.'
                )

        else:
            raise ValueError('Profile not available')

    def __setitem__(self, name: str, item: ndarray):
        """Set the profile directly."""
        if not isinstance(item, (ndarray, Quantity)):
            raise ValueError('"item" must be ndarray or pint Quantity')
        if item.shape[0] != self.n_bins:
            raise ValueError('Length of array does not match number of bins')
        if name in self.loaded_profiles():
            raise ValueError(
                'Attempting to overwrite existing profile. To do so, first delete the '
                'profile\nwith del prof["profile"], then try again.'
            )
        elif name in self.available_profiles():
            raise ValueError(
                'Attempting to set profile already available. '
                'See prof.available_profiles().'
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

    def loaded_profiles(self):
        """Return a listing of loaded profiles."""
        return tuple(sorted(self._profiles.keys()))

    def available_profiles(self):
        """Return a listing of available profiles."""
        loaded = list(self.loaded_profiles())
        available = list(self._profile_functions.keys())
        snap_arrays = list(self.snap.available_arrays())
        return tuple(sorted(set(loaded + available + snap_arrays)))

    def plot(
        self,
        x: str,
        y: Union[str, List[str]],
        x_unit: Optional[str] = None,
        y_unit: Optional[Union[str, List[str]]] = None,
        ax: Any = None,
        **kwargs,
    ):
        """Plot profile.

        Parameters
        ----------
        x
            The x axis to plot as a string.
        y
            The y axis to plot. Can be string or multiple as a list of
            strings.
        x_unit : optional
            The x axis quantity unit as a string. Only works if using
            physical units.
        y_unit : optional
            The y axis quantity unit as a string or list of strings.
            Only works if using physical units.
        ax : optional
            A matplotlib Axes object to plot to.
        **kwargs
            Keyword arguments to pass to Axes plot method.
        """
        if isinstance(y, str):
            y = [y]

        if not self.snap._physical_units:
            if x_unit is not None or y_unit is not None:
                raise ValueError('Cannot set unit if snap is not in physical units')
        else:
            if y_unit is not None:
                if isinstance(y_unit, str):
                    y_unit = [y_unit]
                if len(y) != len(y_unit):
                    raise ValueError('Length of y does not match length of y_unit')

        _x = self[x]
        if x_unit is not None:
            _x = _x.to(x_unit)

        if ax is None:
            fig, ax = plt.subplots()
        else:
            fig = ax.get_figure()

        xlabel = x.capitalize().replace('_', ' ')
        if self.snap._physical_units:
            xlabel = ' '.join([xlabel, f'[{_x.units:~P}]'])
        ax.set_xlabel(xlabel)

        for idx, yi in enumerate(y):
            _y = self[yi]
            label = yi.capitalize().replace('_', ' ')
            if self.snap._physical_units:
                if y_unit is not None:
                    _y = _y.to(y_unit[idx])
                label = ' '.join([label, f'[{_y.units:~P}]'])
                _x = _x.magnitude
                _y = _y.magnitude
            ax.plot(_x, _y, label=label, **kwargs)

        ax.legend()

        return fig, ax

    def to_dataframe(
        self, columns: Union[Tuple[str, ...], List[str]] = None
    ) -> DataFrame:
        """Convert Profile to DataFrame.

        Parameters
        ----------
        columns : optional
            A list of columns to add to the data frame. Default is
            None.

        Returns
        -------
        DataFrame
        """
        data = dict()
        if columns is None:
            columns = self.loaded_profiles()
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
            val = function(*tuple([arg[self._mask][bin_ind] for arg in args]))
            if self.snap._physical_units:
                binned_quantity[idx] = val.magnitude
            else:
                binned_quantity[idx] = val
        if self.snap._physical_units:
            binned_quantity *= val.units
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
def angular_momentum_theta(prof) -> ndarray:
    """Angle between specific angular momentum and xy-plane."""
    angular_momentum_z = prof['angular_momentum_z']
    angular_momentum_magnitude = prof['angular_momentum_magnitude']

    return np.arccos(angular_momentum_z / angular_momentum_magnitude)


@Profile.profile_property
def angular_momentum_phi(prof) -> ndarray:
    """Angle between specific angular momentum and x-axis in xy-plane."""
    angular_momentum_x = prof['angular_momentum_x']
    angular_momentum_y = prof['angular_momentum_y']
    return np.arctan2(angular_momentum_y, angular_momentum_x)


@Profile.profile_property
def toomre_Q(prof) -> ndarray:
    """Toomre Q parameter."""
    if not prof.snap._physical_units:
        G = gravitational_constant_in_code_units(prof.snap)
    else:
        G = (1 * plonk_units.newtonian_constant_of_gravitation).to_base_units()
    return (
        prof['sound_speed']
        * prof['keplerian_frequency']
        / (np.pi * G * prof['density'])
    )


def _str_is_int(string):
    try:
        int(string)
        return True
    except ValueError:
        return False
