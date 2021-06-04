"""Profiles to reduce data to 1-dimension."""

# Heavily inspired by pynbody (https://pynbody.github.io/).

from __future__ import annotations

from bisect import bisect
from copy import copy
from typing import (
    TYPE_CHECKING,
    Any,
    Callable,
    Dict,
    List,
    Optional,
    Sequence,
    Tuple,
    Union,
)

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from numpy import ndarray
from pandas import DataFrame
from scipy.interpolate import interp1d

from .._logging import logger
from .._units import Quantity
from .._units import units as plonk_units
from ..utils.math import average
from ..utils.snap import dust_array_names, vector_array_names
from ..utils.strings import is_documented_by, pretty_array_name
from .extra import extra_profiles

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
        '10 au', or a quantity with units, e.g. plonk.units('10 au').
        Defaults to minimum on the particles.
    cmax : optional
        The maximum coordinate for binning. Can be a string, e.g.
        '10 au', or a quantity with units, e.g. plonk.units('10 au').
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
    >>> prof = plonk.load_profile(snap=snap, cmin='10 au', cmax='300 au')
    >>> prof = plonk.load_profile(snap=snap, spacing='log')

    To access a profile.

    >>> prof['surface_density']
    >>> prof['scale_height']

    To set a new profile.

    >>> prof['aspect_ratio'] = prof['scale_height'] / prof['radius']

    Alternatively use the add_profile decorator.

    >>> @prof.add_profile
    ... def mass(prof):
    ...     M = prof.snap['mass']
    ...     return prof.particles_to_binned_quantity('sum', M)

    Plot one or many quantities on the profile.

    >>> prof.plot('radius', 'density')
    >>> prof.plot('radius', ['angular_momentum_x', 'angular_momentum_y'])

    Plot a quantity on the profile with units.

    >>> units = {'position': 'au', 'surface_density'='g/cm^2'}
    >>> prof.plot('radius', 'surface_density', units=units)
    """

    _profile_aliases: Dict[str, str] = {}

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

        self._profiles: Dict[str, Quantity] = {}
        self._profile_functions: Dict[str, Callable] = {}
        self._default_units = copy(self.snap._default_units)

        self._weights = self.snap['mass']
        self._mask = _setup_particle_mask(snap, ignore_accreted)
        self._x = _calculate_x(snap, self._mask, ndim, coordinate)
        self.range = _set_range(self._x, cmin, cmax)
        self.n_bins = n_bins

        self.bin_edges, self['size'] = _setup_bins(ndim, spacing, self.range, n_bins)
        self.bin_centers = 0.5 * (self.bin_edges[:-1] + self.bin_edges[1:])
        self._x = self._x.to(self.bin_edges.units)
        self._particle_bin = np.digitize(self._x.magnitude, self.bin_edges.magnitude)
        self.bin_indices = _set_particle_bin_indices(self._particle_bin, n_bins)

        if ndim == 1:
            self._coordinate = coordinate
        else:
            self._coordinate = 'radius'
        self._profiles[self._coordinate] = self.bin_centers
        self._profiles['number'] = np.histogram(
            self._x.magnitude, self.bin_edges.magnitude
        )[0] * plonk_units('dimensionless')

        # Add pre-defined profiles
        try:
            num_separate_dust = len(self.snap.num_particles_of_type['dust'])
        except KeyError:
            num_separate_dust = 0
        num_mixture_dust = self.snap.num_dust_species - num_separate_dust
        extra_profiles(self, num_separate_dust, num_mixture_dust)

    def add_profile(self, fn: Callable) -> Callable:
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
        self._profile_functions[fn.__name__] = fn
        return fn

    def add_alias(self, name: str, alias: str) -> None:
        """Add alias to array.

        Parameters
        ----------
        name
            The name of the array.
        alias
            The alias to reference the array.
        """
        self._profile_aliases[alias] = name

    def loaded_profiles(self) -> List[str]:
        """Return a listing of loaded profiles."""
        return sorted(self._profiles.keys())

    def available_profiles(self) -> List[str]:
        """Return a listing of available profiles."""
        loaded = list(self.loaded_profiles())
        available = list(self._profile_functions.keys())
        snap_arrays = _1d_arrays(list(self.snap.available_arrays(verbose=True)))
        return sorted(set(loaded + available + snap_arrays))

    @property
    def default_units(self) -> Dict[str, Any]:
        """Profile default units."""
        return {
            key: self._default_units[key] for key in sorted(self._default_units.keys())
        }

    def set_units(self, **kwargs) -> Profile:
        """Set default unit for profiles.

        Parameters
        ----------
        kwargs
            Keyword arguments with keys as the profile name, e.g.
            'pressure', and with values as the unit as a string, e.g.
            'pascal'.

        Examples
        --------
        Set multiple default units.

        >>> profile.set_units(pressure='pascal', density='g/cm^3')
        """
        for key, val in kwargs.items():
            defaults = list(self.default_units) + list(self.snap.default_units)
            if key not in defaults:
                logger.info(f'adding profile {key} to default_units dict')
            self._default_units[key] = val

        return self

    def base_profile_name(self, name: str) -> str:
        """Get the base profile name from a string.

        For example, 'velocity_x' returns 'velocity', 'density' returns
        'density', 'dust_fraction_001' returns 'dust_fraction', 'x'
        returns 'position'.

        Parameters
        ----------
        name
            The name as a string

        Returns
        -------
        str
            The base name.
        """
        try:
            return self.snap.base_array_name(name)
        except ValueError:
            if name == self._coordinate:
                return 'position'
            name_root = '_'.join(name.split('_')[:-1])
            name_suffix = name.split('_')[-1]
            if name_root == '' and name_suffix in ('x', 'y', 'z'):
                return 'position'
            if name_root in self.snap._array_aliases:
                return self.snap._array_aliases[name_root]
            if name_suffix in ('x', 'y', 'z', 'mag'):
                return name_root
            if _str_is_int(name_suffix):
                return name_root
            return name

    def plot(
        self,
        x: str,
        y: Union[str, List[str]],
        units: Dict[str, Union[str, List[str]]] = None,
        std: str = None,
        label: Union[str, List[str]] = None,
        ax: Any = None,
        ax_kwargs={},
        **kwargs,
    ) -> Any:
        """Plot profile.

        Parameters
        ----------
        x
            The x axis to plot as a string.
        y
            The y axis to plot. Can be string or a list of strings.
        units
            The units of the plot as a dictionary. The keys correspond
            to quantities such as 'position', 'density', 'velocity', and
            so on. The values are strings representing units, e.g.
            'g/cm^3' for density.
        std : optional
            Add standard deviation on profile. Can be 'shading' or
            'errorbar'.
        label : optional
            A label for the plot. Can be a string or a list of strings,
            one per y.
        ax : optional
            A matplotlib Axes object to plot to.
        ax_kwargs
            Keyword arguments to pass to matplotlib Axes.
        **kwargs
            Keyword arguments to pass to Axes plot method.

        Returns
        -------
        ax
            The matplotlib Axes object.
        """
        if std is not None and std not in ('shading', 'errorbar'):
            raise ValueError('std must be "shading" or "errorbar"')

        ynames = _yname_from_yinput(y, self)
        labels = _labels(label, ynames)

        xunit = _get_unit(self, x, units)
        yunits = [_get_unit(self, y, units) for y in ynames]

        xdata = self[x].to(xunit)

        if ax is None:
            _, ax = plt.subplots()

        for yname, yunit, label in zip(ynames, yunits, labels):
            ydata = self[yname].to(yunit)
            if label is None:
                label = f'{pretty_array_name(yname)}'
                if ydata.units != plonk_units.dimensionless:
                    label = f'{pretty_array_name(yname)}'
                    label += f' [{ydata.units:~P}]'
            [line] = ax.plot(xdata.magnitude, ydata.magnitude, label=label, **kwargs)
            if std:
                color = line.get_color()
                _std_plot(self, xdata, ydata, yname, yunit, std, color, ax)

        ax.set_xlabel(f'{pretty_array_name(x)} [{xdata.units:~P}]')
        ax.legend()
        ax.set(**ax_kwargs)

        return ax

    def to_function(self, profile: str, **kwargs) -> Callable:
        """Create function via interpolation.

        The function is of the coordinate of the profile, e.g.
        'radius', and returns values of the selected profile, e.g.
        'scale_height'. The function is generated from the
        scipy.interpolate function interp1d.

        Parameters
        ----------
        profile
            The profile function to create as a string, e.g.
            'scale_height'.

        Returns
        -------
        Callable
            The function.

        Examples
        --------
        Select all particles within a scale height in a disc.

        >>> scale_height = prof.to_function('scale_height')
        >>> subsnap = snap[np.abs(snap['z']) < scale_height(snap['R'])]
        """
        coord = self.bin_centers
        prof = self[profile]

        def fn(x):
            nonlocal coord, prof
            _coord = coord.to(x.units).magnitude
            _prof = prof.magnitude
            y = interp1d(_coord, _prof, fill_value='extrapolate', **kwargs)(x.magnitude)
            return y * self[profile].units

        return fn

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
                    _units.append(plonk_units('dimensionless'))
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

    def particles_to_binned_quantity(
        self, aggregation: str, array: Quantity
    ) -> Quantity:
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

        _array = array[self._mask]
        binned_quantity = np.zeros(self.n_bins) * _array.units
        for idx, bin_ind in enumerate(self.bin_indices):
            if bin_ind.size == 0:
                continue
            if aggregation == 'average':
                val = average(
                    _array[bin_ind],
                    weights=self._weights[self._mask][bin_ind],
                )
            elif aggregation == 'mean':
                val = np.mean(_array[bin_ind])
            elif aggregation == 'median':
                val = np.median(_array[bin_ind])
            elif aggregation == 'std':
                val = np.std(_array[bin_ind])
            elif aggregation == 'sum':
                val = np.sum(_array[bin_ind])

            binned_quantity[idx] = val

        return binned_quantity

    def __getitem__(self, name: str) -> Quantity:
        """Return the profile of a given kind."""
        name_root = '_'.join(name.split('_')[:-1])
        name_suffix = name.split('_')[-1]

        if name in self._profile_aliases:
            name = self._profile_aliases[name]
        if name in self._profiles:
            return self._profiles[name]
        if name in self._profile_functions:
            self._profiles[name] = self._profile_functions[name](self)
            return self._profiles[name]

        if name_suffix in _aggregations:
            aggregation = name_suffix
            array_name = name_root
        else:
            aggregation = self.aggregation
            array_name = name
        try:
            array: Quantity = self.snap[array_name]
        except ValueError as e:
            logger.error(e)
            raise ValueError(f'array "{array_name}" not available on snap')
        if array.ndim == 1:
            self._profiles[name] = self.particles_to_binned_quantity(aggregation, array)
            return self._profiles[name]
        raise ValueError(
            'Requested profile has array dimension > 1.\nTo access x-, y-, or '
            'z-components, or magnitude of vector quantities,\ntry, for '
            'example, prof["velocity_x"] or prof["momentum_mag"].\nTo '
            'access dust profiles, try, for example, prof["stopping_time_001"]'
        )

    def __setitem__(self, name: str, item: Quantity):
        """Set the profile directly."""
        if not isinstance(item, Quantity):
            raise ValueError('"item" must be pint Quantity')
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


def _setup_particle_mask(snap: SnapLike, ignore_accreted: bool) -> ndarray:
    if ignore_accreted is False:
        return np.ones(len(snap), dtype=bool)
    h: Quantity = snap['h']
    return h > 0


def _calculate_x(snap: SnapLike, mask: ndarray, ndim: int, coordinate: str) -> Quantity:
    pos: Quantity = snap['position']
    pos = pos[mask]
    if ndim == 1:
        if coordinate == 'x':
            return pos[:, 0]
        if coordinate == 'y':
            return pos[:, 1]
        if coordinate == 'z':
            return pos[:, 2]
        raise ValueError('coordinate must be "x", "y", or "z" for ndim==1')
    if ndim == 2:
        return np.sqrt(pos[:, 0] ** 2 + pos[:, 1] ** 2)
    if ndim == 3:
        return np.sqrt(pos[:, 0] ** 2 + pos[:, 1] ** 2 + pos[:, 2] ** 2)
    raise ValueError('Unknown ndim: cannot calculate x array')


def _set_range(x: Quantity, cmin: Any, cmax: Any) -> Tuple[float, float]:
    if cmin is None:
        rmin = x.min()
    else:
        rmin = Quantity(cmin)
        if not rmin.dimensionality == Quantity('cm').dimensionality:
            raise ValueError('must specify cmin units, e.g. cmin="10 au"')
    if cmax is None:
        rmax = np.percentile(x.magnitude, 99, axis=0) * x.units
    else:
        rmax = Quantity(cmax)
        if not rmax.dimensionality == Quantity('cm').dimensionality:
            raise ValueError('must specify cmin units, e.g. cmax="100 au"')

    return rmin, rmax


def _setup_bins(ndim: int, spacing: str, xrange: Quantity, n_bins: int) -> Quantity:
    bin_edges = _bin_edges(spacing, xrange, n_bins)
    if ndim == 1:
        bin_sizes = bin_edges[1:] - bin_edges[:-1]
    elif ndim == 2:
        bin_sizes = np.pi * (bin_edges[1:] ** 2 - bin_edges[:-1] ** 2)
    elif ndim == 3:
        bin_sizes = 4 / 3 * np.pi * (bin_edges[1:] ** 3 - bin_edges[:-1] ** 3)
    return bin_edges, bin_sizes


def _bin_edges(spacing: str, xrange: Quantity, n_bins: int) -> Quantity:
    if spacing == 'linear':
        bin_edges = (
            np.linspace(xrange[0].magnitude, xrange[1].magnitude, n_bins + 1)
            * xrange[0].units
        )
    elif spacing == 'log':
        bin_edges = (
            np.logspace(
                np.log10(xrange[0].magnitude),
                np.log10(xrange[1].magnitude),
                n_bins + 1,
            )
            * xrange[0].units
        )
    else:
        raise ValueError('Cannot determine spacing to setup bins')
    return bin_edges


def _set_particle_bin_indices(particle_bin: ndarray, n_bins: int) -> List[ndarray]:
    sortind = particle_bin.argsort()
    sort_pind = particle_bin[sortind]
    binind = list()
    prev_index = bisect(sort_pind, 0)
    for i in range(n_bins):
        new_index = bisect(sort_pind, i + 1)
        binind.append(np.sort(sortind[prev_index:new_index]))
        prev_index = new_index
    return binind


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


def _1d_arrays(arrays: list) -> list:
    _arrays = list()
    for array in arrays:
        if (array + '_x' in arrays) or (array + '_001' in arrays):
            pass
        else:
            _arrays.append(array)
    return _arrays


def _get_unit(profile, name, units):
    if name is None:
        return None
    base_name = profile.base_profile_name(name)
    if units is not None:
        if name in units:
            return 1 * plonk_units(units[name])
        if base_name in units:
            return 1 * plonk_units(units[base_name])
    if name in profile.default_units:
        return 1 * plonk_units(profile.default_units[name])
    if base_name in profile.default_units:
        return 1 * plonk_units(profile.default_units[base_name])
    return 1 * profile[name].units


def _yname_from_yinput(y, profile):
    if isinstance(y, str):
        if y not in profile.available_profiles():
            if y + '_001' in profile.available_profiles():
                return dust_array_names(
                    name=y, num_dust_species=profile.snap.num_dust_species
                )
            elif y + '_x' in profile.available_profiles():
                return vector_array_names(name=y)
            return [y]
        return [y]
    return y


def _labels(label, ynames):
    labels: Sequence[Optional[str]]
    if label is not None:
        if isinstance(label, str):
            return [label for _ in ynames]
        return label
    else:
        return [None for _ in ynames]


def _std_plot(profile, xdata, ydata, yname, yunit, std, color, ax):
    if yname.split('_')[-1] in _aggregations:
        _yname = '_'.join(yname.split('_')[:-1])
    else:
        _yname = yname
    try:
        y_std = profile[_yname + '_std']
    except ValueError:
        logger.warning('Cannot calculate standard deviation')
        return
    if profile.aggregation in ('std', 'sum'):
        y_mean = profile[_yname + '_mean']
    else:
        y_mean = ydata
    y_std = y_std.to(yunit).magnitude
    y_mean = y_mean.to(yunit).magnitude
    xdata = xdata.magnitude
    if std == 'shading':
        ax.fill_between(
            xdata,
            y_mean - y_std,
            y_mean + y_std,
            color=color,
            alpha=0.2,
        )
    elif std == 'errorbar':
        ax.errorbar(xdata, y_mean, yerr=y_std, linestyle='', color=color, alpha=0.5)


def _str_is_int(string: str) -> bool:
    try:
        int(string)
        return True
    except ValueError:
        return False
