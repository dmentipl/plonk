"""Profiles.

Heavily inspired by pynbody (https://pynbody.github.io/).
"""

from __future__ import annotations

from bisect import bisect
from functools import partial
from typing import TYPE_CHECKING, Any, Callable, Dict, List, Tuple, Union

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from numpy import ndarray
from pandas import DataFrame

from .._logging import logger
from .._units import Quantity
from .._units import units as plonk_units
from ..utils.math import average
from ..utils.utils import is_documented_by

if TYPE_CHECKING:
    from ..snap.snap import SnapLike

_aggregations = ('average', 'mean', 'median', 'std', 'sum')


class Profile:
    """Profiles.

    A profile is a binning of particles in Cartesian slices, or
    cylindrical or spherical shells around the origin, i.e. (0, 0, 0).
    For cylindrical profiles the cylindrical cells are perpendicular to
    the xy-plane, i.e. the bins are averaged azimuthally and in the
    z-direction. Cartesian profiles can be in the x-, y-, and
    z-directions.

    Parameters
    ----------
    snap
        The Snap object.
    ndim : optional
        The dimension of the profile. For ndim == 2, the radial binning
        is cylindrical in the xy-plane. For ndim == 3, the radial
        binning is spherical. For ndim == 1, the radial binning is
        Cartesian along the x-axis. Default is 2.
    cmin : optional
        The minimum coordinate for binning. Can be a string, e.g.
        '10 au', or a quantity with units, e.g. plonk.units['10 au'].
        Defaults to minimum on the particles.
    cmax : optional
        The maximum coordinate for binning. Can be a string, e.g.
        '10 au', or a quantity with units, e.g. plonk.units['10 au'].
        Defaults to the 99 percentile distance.
    n_bins : optional
        The number of radial bins. Default is 100.
    aggregation : optional
        The method to aggregate particle quantities in bins by. Options
        are 'average', 'mean', or 'median'. Here 'average' is a
        mass-weighted average. Default is 'average'.
    spacing : optional
        The spacing of radial bins. Can be 'linear' or 'log'. Default is
        'linear'.
    coordinate : optional
        The coordinate ('x', 'y', or 'z') for Cartesian profiles only,
        i.e. when ndim==1. Default is 'x'. For cylindrical and spherical
        profiles the coordinate is 'radius'.
    ignore_accreted : optional
        Ignore particles accreted onto sinks. Default is True.

    Examples
    --------
    Generate profile from snapshot.

    >>> prof = plonk.load_profile(snap=snap)
    >>> prof = plonk.load_profile(snap=snap, n_bins=300)
    >>> prof = plonk.load_profile(
    ...     snap=snap, cmin='10 au', cmax='300 au'
    ... )
    >>> prof = plonk.load_profile(snap=snap, spacing='log')

    To access a profile.

    >>> prof['surface_density']
    >>> prof['scale_height']

    To set a new profile.

    >>> prof['aspect_ratio'] = prof['scale_height'] / prof['radius']

    Alternatively use the profile_property decorator.

    >>> @Profile.profile_property
    ... def mass(prof):
    ...     M = prof.snap['mass']
    ...     return prof.particles_to_binned_quantity('sum', M)

    Plot one or many quantities on the profile.

    >>> prof.plot('radius', 'density')
    >>> prof.plot(
    ...     'radius', ['angular_momentum_x', 'angular_momentum_y'],
    ... )

    Plot a quantity on the profile with units.

    >>> prof.plot(
    ...     'radius', 'surface_density', x_unit='au', y_unit='g/cm^2'
    ... )
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
        snap: SnapLike,
        ndim: int = 2,
        cmin: Any = None,
        cmax: Any = None,
        n_bins: int = 100,
        aggregation: str = 'average',
        spacing: str = 'linear',
        coordinate: str = 'x',
        ignore_accreted: bool = True,
    ):

        self.snap = snap
        self.ndim = ndim
        self.aggregation = _check_aggregation(aggregation)
        self.spacing = _check_spacing(spacing)
        self.properties: Dict[str, Any] = {}

        self._profiles: Dict[str, ndarray] = {}

        self._weights = self.snap['mass']
        self._mask = self._setup_particle_mask(ignore_accreted)
        self._x = self._calculate_x(coordinate)
        self.range = self._set_range(cmin, cmax)
        self.n_bins = n_bins

        self.bin_edges, self['size'] = self._setup_bins()
        self.bin_centers = 0.5 * (self.bin_edges[:-1] + self.bin_edges[1:])
        self._particle_bin = np.digitize(self._x.magnitude, self.bin_edges.magnitude)
        self.bin_indicies = self._set_particle_bin_indicies()

        if ndim == 1:
            self._profiles[coordinate] = self.bin_centers
        else:
            self._profiles['radius'] = self.bin_centers
        self._profiles['number'] = np.histogram(
            self._x.magnitude, self.bin_edges.magnitude
        )[0]

        n_dust = len(self.snap.properties.get('grain_size', []))
        _generate_profiles(n_dust)

    def _setup_particle_mask(self, ignore_accreted: bool) -> ndarray:
        if ignore_accreted is False:
            return np.ones(len(self.snap), dtype=bool)
        h: ndarray = self.snap['h']
        return h > 0

    def _calculate_x(self, coordinate) -> ndarray:
        pos: ndarray = self.snap['xyz']
        pos = pos[self._mask]
        if self.ndim == 1:
            if coordinate == 'x':
                return pos[:, 0]
            if coordinate == 'y':
                return pos[:, 1]
            if coordinate == 'z':
                return pos[:, 2]
            raise ValueError('coordinate must be "x", "y", or "z" for ndim==1')
        if self.ndim == 2:
            return np.sqrt(pos[:, 0] ** 2 + pos[:, 1] ** 2)
        if self.ndim == 3:
            return np.sqrt(pos[:, 0] ** 2 + pos[:, 1] ** 2 + pos[:, 2] ** 2)
        raise ValueError('Unknown ndim: cannot calculate x array')

    def _set_range_code_units(self, cmin: float, cmax: float) -> Tuple[float, float]:
        if cmin is None:
            rmin = self._x.min()
        else:
            rmin = cmin
        if cmax is None:
            rmax = np.percentile(self._x, 99, axis=0)
        else:
            rmax = cmax

        return rmin, rmax

    def _set_range(self, cmin: Any, cmax: Any) -> Tuple[float, float]:
        if cmin is None:
            rmin = self._x.min()
        else:
            rmin = Quantity(cmin)
            if not rmin.dimensionality == Quantity('cm').dimensionality:
                raise ValueError('must specify cmin units, e.g. cmin="10 au"')
            rmin = rmin.to_base_units()
        if cmax is None:
            rmax = np.percentile(self._x.magnitude, 99, axis=0) * self._x.units
        else:
            rmax = Quantity(cmax)
            if not rmax.dimensionality == Quantity('cm').dimensionality:
                raise ValueError('must specify cmin units, e.g. cmax="100 au"')
            rmax = rmax.to_base_units()

        return rmin, rmax

    def _setup_bins(self) -> ndarray:
        bin_edges = self._bin_edges()
        if self.ndim == 1:
            bin_sizes = bin_edges[1:] - bin_edges[:-1]
        elif self.ndim == 2:
            bin_sizes = np.pi * (bin_edges[1:] ** 2 - bin_edges[:-1] ** 2)
        elif self.ndim == 3:
            bin_sizes = 4 / 3 * np.pi * (bin_edges[1:] ** 3 - bin_edges[:-1] ** 3)
        return bin_edges, bin_sizes

    def _bin_edges(self):
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

    def _getitem(self, name: str) -> ndarray:
        """Return the profile of a given kind."""
        name_root = '_'.join(name.split('_')[:-1])
        name_suffix = name.split('_')[-1]

        if name in self._profiles:
            return self._profiles[name]
        if name in Profile._profile_functions:
            self._profiles[name] = Profile._profile_functions[name](self)
            return self._profiles[name]

        if name_suffix in _aggregations:
            aggregation = name_suffix
            array_name = name_root
        else:
            aggregation = self.aggregation
            array_name = name
        try:
            array: ndarray = self.snap[array_name]
            if array.ndim == 1:
                self._profiles[name] = self.particles_to_binned_quantity(
                    aggregation, array
                )
                return self._profiles[name]
            raise ValueError(
                'Requested profile has array dimension > 1.\nTo access x-, y-, or '
                'z-components, or magnitude of vector quantities,\ntry, for '
                'example, prof["velocity_x"] or prof["momentum_magnitude"].\nTo '
                'access dust profiles, try, for example, prof["stopping_time_001"] '
                'or\nprof["dust_mass_sum"].'
            )
        except ValueError:
            logger.error('Profile unavailable.')

    def __getitem__(self, name: str) -> ndarray:
        """Return the profile of a given kind."""
        return self._getitem(name)

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
        if name in self.available_profiles():
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
        """Dunder repr method."""
        return self.__str__()

    def __str__(self):
        """Dunder str method."""
        return f'<plonk.Profile "{self.snap.file_path.name}">'

    def _ipython_key_completions_(self):
        """Tab completion for IPython __getitem__ method."""
        return self.available_profiles()

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
        x_unit: str = None,
        y_unit: Union[str, List[str]] = None,
        std_dev_shading: bool = False,
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
            The x axis quantity unit as a string.
        y_unit : optional
            The y axis quantity unit as a string or list of strings.
        std_dev_shading : optional
            Add shading for standard deviation of profile.
        ax : optional
            A matplotlib Axes object to plot to.
        **kwargs
            Keyword arguments to pass to Axes plot method.

        Returns
        -------
        ax
            The matplotlib Axes object.
        """
        if isinstance(y, str):
            y = [y]

        if y_unit is not None:
            if isinstance(y_unit, str):
                y_unit = [y_unit]
            if len(y) != len(y_unit):
                raise ValueError('Length of y does not match length of y_unit')

        _x = self[x]
        if x_unit is not None:
            _x = _x.to(x_unit)

        if ax is None:
            _, ax = plt.subplots()

        xlabel = x.capitalize().replace('_', ' ')
        xlabel = ' '.join([xlabel, f'[{_x.units:~P}]'])
        _x = _x.magnitude
        ax.set_xlabel(xlabel)

        for idx, yi in enumerate(y):
            _y = self[yi]
            if std_dev_shading:
                if yi.split('_')[-1] in _aggregations:
                    _yi = '_'.join(yi.split('_')[:-1])
                else:
                    _yi = yi
                try:
                    _y_std = self[_yi + '_std']
                except ValueError:
                    logger.warning('Cannot calculate standard deviation')
                if self.aggregation != 'mean':
                    _y_mean = self[_yi + '_mean']
                else:
                    _y_mean = _y
            label = yi.capitalize().replace('_', ' ')
            if y_unit is not None:
                _y = _y.to(y_unit[idx])
                if std_dev_shading:
                    _y_std = _y_std.to(y_unit[idx])
                    _y_mean = _y_mean.to(y_unit[idx])
            label = ' '.join([label, f'[{_y.units:~P}]'])
            _y = _y.magnitude
            if std_dev_shading:
                _y_std = _y_std.magnitude
                _y_mean = _y_mean.magnitude
            if kwargs.get('label') is None:
                lines = ax.plot(_x, _y, label=label, **kwargs)
            else:
                lines = ax.plot(_x, _y, **kwargs)
            if std_dev_shading:
                color = lines[0].get_color()
                ax.fill_between(
                    _x, _y_mean - _y_std, _y_mean + _y_std, color=color, alpha=0.2,
                )

        ax.legend()

        return ax

    def to_dataframe(
        self, columns: List[str] = None, units: List[str] = None
    ) -> DataFrame:
        """Convert Profile to DataFrame.

        Parameters
        ----------
        columns : optional
            A list of columns to add to the data frame. If None, add all
            loaded columns. Default is None.
        units : optional
            A list of units corresponding to columns add to the data
            frame. Units must be strings, and must be base units. I.e.
            'cm' not '10 cm'. If None, use default, i.e. cgs. Default is
            None.

        Returns
        -------
        DataFrame
        """
        data = dict()
        if columns is None:
            columns = self.loaded_profiles()
        if units is None:
            _units = list()
            for column in columns:
                try:
                    _units.append(self[column].units)
                except AttributeError:
                    _units.append(plonk_units['dimensionless'])
        else:
            _units = list()
            for unit in units:
                u = plonk_units(unit)
                if np.allclose(u.m, 1.0):
                    _units.append(u.units)
                else:
                    raise ValueError(
                        'Units must be strings, and must be base units. '
                        'I.e. "cm" not "10 cm".'
                    )
        if len(_units) != len(columns):
            raise ValueError('units and columns must have same length')
        for column, unit in zip(columns, _units):
            try:
                name = column + f' [{unit:~}]'
                array = self[column].to(unit).magnitude
            except AttributeError:
                name = column
                array = self[column]
            data[name] = array
        return pd.DataFrame(data)

    def particles_to_binned_quantity(self, aggregation: str, array: ndarray) -> ndarray:
        """Calculate binned quantities from particles.

        This takes care of the bin indices and ignoring accreted
        particles (if requested in instantiating the profile).

        Parameters
        ----------
        aggregation
            The aggregation function that acts on particles in the
            radial bin.
        array
            The particle array.
        """
        if aggregation not in _aggregations:
            raise ValueError('Cannot determine aggregation method')

        binned_quantity = np.zeros(self.n_bins)
        for idx, bin_ind in enumerate(self.bin_indicies):
            if bin_ind.size == 0:
                continue
            if aggregation == 'average':
                val = average(
                    array[self._mask][bin_ind],
                    weights=self._weights[self._mask][bin_ind],
                )
            elif aggregation == 'mean':
                val = np.mean(array[self._mask][bin_ind])
            elif aggregation == 'median':
                val = np.median(array[self._mask][bin_ind])
            elif aggregation == 'std':
                val = np.std(array[self._mask][bin_ind])
            elif aggregation == 'sum':
                val = np.sum(array[self._mask][bin_ind])

            binned_quantity[idx] = val.magnitude

        binned_quantity *= val.units

        return binned_quantity


def _generate_profiles(n_dust: int = 0):
    @Profile.profile_property
    def mass(prof) -> ndarray:
        """Mass profile."""
        M = prof.snap['mass']
        return prof.particles_to_binned_quantity('sum', M)

    @Profile.profile_property
    def surface_density(prof) -> ndarray:
        """Surface density profile.

        Units are [mass / length ** ndim], which depends on ndim of profile.
        """
        return prof['mass'] / prof['size']

    @Profile.profile_property
    def scale_height(prof) -> ndarray:
        """Scale height profile."""
        z = prof.snap['z']
        return prof.particles_to_binned_quantity('std', z)

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
        angular_momentum_magnitude = prof['angular_momentum_mag']

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
        G = (1 * plonk_units.newtonian_constant_of_gravitation).to_base_units()
        return (
            prof['sound_speed']
            * prof['keplerian_frequency']
            / (np.pi * G * prof['density'])
        )

    if n_dust > 0:

        @Profile.profile_property
        def gas_mass(prof) -> ndarray:
            """Gas mass profile."""
            M = prof.snap['gas_mass']
            return prof.particles_to_binned_quantity('sum', M)

        @Profile.profile_property
        def gas_surface_density(prof) -> ndarray:
            """Gas surface density profile.

            Units are [mass / length ** ndim], which depends on ndim of profile.
            """
            return prof['gas_mass'] / prof['size']

        for idx in range(n_dust):

            def dust_mass(idx, prof) -> ndarray:
                """Dust mass profile."""
                M = prof.snap[f'dust_mass_{idx+1:03}']
                return prof.particles_to_binned_quantity('sum', M)

            def dust_surface_density(idx, prof) -> ndarray:
                """Dust surface density profile.

                Units are [mass / length ** ndim], which depends on ndim of profile.
                """
                return prof[f'dust_mass_{idx+1:03}'] / prof['size']

            Profile._profile_functions[f'dust_mass_{idx+1:03}'] = partial(
                dust_mass, idx
            )

            Profile._profile_functions[f'dust_surface_density_{idx+1:03}'] = partial(
                dust_surface_density, idx
            )


@is_documented_by(Profile)
def load_profile(
    snap: SnapLike,
    ndim: int = 2,
    cmin: Any = None,
    cmax: Any = None,
    n_bins: int = 100,
    aggregation: str = 'average',
    spacing: str = 'linear',
    coordinate: str = 'x',
    ignore_accreted: bool = True,
) -> Profile:
    logger.debug(f'Loading profile: {snap.file_path.name}')
    return Profile(
        snap=snap,
        ndim=ndim,
        cmin=cmin,
        cmax=cmax,
        n_bins=n_bins,
        aggregation=aggregation,
        spacing=spacing,
        coordinate=coordinate,
        ignore_accreted=ignore_accreted,
    )


def _check_aggregation(method: str) -> str:
    if method in _aggregations:
        return method
    raise ValueError(
        f'Cannot determine aggregation method: choose from {_aggregations}'
    )


def _check_spacing(spacing: str) -> str:
    if spacing.lower() in ('lin', 'linear'):
        return 'linear'
    if spacing.lower() in ('log', 'logarithm', 'logarithmic'):
        return 'log'
    raise ValueError('Cannot determine spacing')
