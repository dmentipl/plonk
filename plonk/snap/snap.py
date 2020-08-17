"""Snap, SubSnap, Sinks classes for snapshot files.

The Snap class contains all information related to a smoothed particle
hydrodynamics simulation snapshot file. The SubSnap class is for
accessing a subset of particles in a Snap.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Callable, Dict, List, Tuple, Union

import h5py
import numpy as np
import pandas as pd
from numpy import ndarray
from pandas import DataFrame
from scipy.spatial import cKDTree
from scipy.spatial.transform import Rotation

from .. import visualize
from .._logging import logger
from .._units import Quantity
from .._units import units as plonk_units
from ..utils.kernels import kernel_names, kernel_radius
from ..utils.math import norm
from . import context


class Snap:
    """Smoothed particle hydrodynamics Snap object.

    Snapshot files contain the state of the simulation at a point in
    time. Typical minimum data from a smoothed particle hydrodynamics
    simulation include the particle positions and smoothing length, from
    which the density field can be reconstructed, as well as the
    particle type. In addition, the particle velocities are required to
    restart the simulation.

    Other data stored in the snapshot file include equation of state,
    dust, and magnetic field information, as well as numerical
    quantities related to time-stepping.

    Examples
    --------
    Use load_snap to generate a Snap object.

    >>> snap = plonk.load_snap('file_name.h5')

    To access arrays on the particles.

    >>> snap['position']
    >>> snap['density']

    To access sink arrays.

    >>> snap.sinks['position']
    >>> snap.sinks['spin']

    To access a subset of particles as a SubSnap.

    >>> subsnap = snap[:100]
    >>> subsnap = snap[snap['x'] > 0]
    >>> subsnap = snap['gas']

    To set a new array directly.

    >>> snap['my_r'] = np.sqrt(snap['x'] ** 2 + snap['y'] ** 2)

    Alternatively, define a function.

    >>> @snap.add_array()
    ... def my_radius(snap):
    ...     return np.hypot(snap['x'], snap['y'])

    Or, use an existing one function.

    >>> @snap.add_array()
    ... def my_R(snap):
    ...     return plonk.analysis.particles.radial_distance(snap)
    """

    _array_aliases = {
        'B': 'magnetic_field',
        'c_s': 'sound_speed',
        'dt': 'timestep',
        'e': 'eccentricity',
        'eps': 'dust_fraction',
        'h': 'smoothing_length',
        'j': 'specific_angular_momentum',
        'L': 'angular_momentum',
        'm': 'mass',
        'omega_k': 'keplerian_frequency',
        'p': 'momentum',
        'P': 'pressure',
        'pos': 'position',
        'phi': 'azimuthal_angle',
        'R': 'radius_cylindrical',
        'r': 'radius_spherical',
        'rho': 'density',
        'rho_d': 'dust_density',
        'rho_g': 'gas_density',
        'St': 'stokes_number',
        'T': 'temperature',
        'theta': 'polar_angle',
        't_s': 'stopping_time',
        'u': 'internal_energy',
        'v': 'velocity',
        'v_R': 'velocity_radial_cylindrical',
        'v_r': 'velocity_radial_spherical',
        'v_phi': 'angular_velocity',
        'velocity_R': 'radial_velocity_cylindrical',
        'velocity_r': 'radial_velocity_spherical',
        'velocity_phi': 'angular_velocity',
    }

    _vector_arrays = {
        'magnetic_field',
        'position',
        'spin',
        'velocity',
    }

    _dust_arrays = {
        'stopping_time',
        'dust_fraction',
        'dust_to_gas_ratio',
    }

    particle_type = {
        'gas': 1,
        'dust': 2,
        'boundary': 3,
        'star': 4,
        'darkmatter': 5,
        'bulge': 6,
    }

    @staticmethod
    def add_alias(name: str, alias: str) -> None:
        """Add alias to array.

        Parameters
        ----------
        name
            The name of the array.
        alias
            The alias to reference the array.
        """
        Snap._array_aliases[alias] = name

    def __init__(self):

        self.data_source = None
        self.file_path = None
        self.code_units = {}
        self._properties = {}
        self._array_code_units = {}
        self._array_registry: Dict[str, Callable] = {}
        self._sink_registry: Dict[str, Callable] = {}
        self._cache_arrays = True
        self._arrays = {}
        self._sinks = {}
        self._file_pointer = None
        self._num_particles = -1
        self._num_particles_of_type = -1
        self._num_sinks = -1
        self._num_dust_species = -1
        self.rotation = None
        self.translation = None
        self._tree = None

    def close_file(self):
        """Close access to underlying file."""
        self._file_pointer.close()

    def reopen_file(self):
        """Re-open access to the underlying file."""
        self._file_pointer = h5py.File(self.file_path, mode='r')

    def add_array(self, vector: bool = False, dust: bool = False) -> Callable:
        """Decorate function to add array to Snap.

        This function decorates a function that returns an array. The
        name of the function is the string with which to reference the
        array.

        The function being decorated should return a Pint Quantity
        array, not a unitless numpy array.

        Parameters
        ----------
        vector
            A bool to represent if the array is a vector in space that
            should have rotations applied to it. If True the rotation
            should be applied. If False the rotation cannot be applied.
            Default is False.
        dust
            A bool to represent if the array is a dust array, in that
            it has one column per dust species. Default is False.

        Returns
        -------
        Callable
            The decorator which returns the array.
        """

        def _add_array(fn):
            self._array_registry[fn.__name__] = fn
            if vector is True:
                self._vector_arrays.add(fn.__name__)
            if dust is True:
                self._dust_arrays.add(fn.__name__)
            return fn

        return _add_array

    def add_unit(self, name: str, unit: str) -> Snap:
        """Add missing code unit to array.

        Add code units to an array from file that has not had units set
        automatically.

        Parameters
        ----------
        name
            The name of the array
        unit
            A unit string representing the array code unit, e.g.
            '1.234 kg * m / s'.
        """
        self._array_code_units[name] = plonk_units(unit)
        if name in self._arrays:
            del self[name]

        return self

    def loaded_arrays(self) -> List[str]:
        """Return a list of loaded arrays.

        Returns
        -------
        List
            A list of names of loaded particle arrays.
        """
        return sorted(self._arrays.keys())

    def _available_arrays(
        self, sinks: bool = False, verbose: bool = False, aliases: bool = False
    ) -> List[str]:
        """Return a list of available arrays.

        Parameters
        ----------
        sinks
            If True, return available sink arrays. Default is
            False.
        verbose
            Also display suffixed arrays, e.g. 'position_x',
            'position_y', etc. Default is False
        aliases
            If True, return array aliases. Default is False.

        Returns
        -------
        List
            A list of names of available arrays.
        """
        if sinks:
            loaded = self.sinks.loaded_arrays()
            registered = list(self._sink_registry.keys())
        else:
            loaded = self.loaded_arrays()
            registered = list(self._array_registry.keys())
        if verbose:
            for arr in self._vector_arrays:
                if arr in registered:
                    registered += [f'{arr}_x', f'{arr}_y', f'{arr}_z', f'{arr}_mag']
            for arr in self._dust_arrays:
                if arr in registered:
                    registered += [
                        f'{arr}_{n:03}' for n in range(1, self.num_dust_species + 1)
                    ]
                    registered.append(f'{arr}_tot')

        if aliases:
            extra = [
                key
                for key, val in self._array_aliases.items()
                if val[0] in self.loaded_arrays() or val in self._array_registry
            ]
            if verbose:
                _extra = list()
                for arr in extra:
                    if self.base_array_name(arr) in self._vector_arrays:
                        for suffix in ['x', 'y', 'z', 'mag']:
                            _extra.append(arr + f'_{suffix}')
                    if self.base_array_name(arr) in self._dust_arrays:
                        suffixes = [f'{n+1:03}' for n in range(self.num_dust_species)]
                        suffixes += ['tot']
                        for suffix in suffixes:
                            _extra.append(arr + f'_{suffix}')
                extra += _extra
            return sorted(set(extra), key=lambda x: x.lower())

        return sorted(set(loaded + registered))

    def available_arrays(
        self, verbose: bool = False, aliases: bool = False
    ) -> List[str]:
        """Return a list of available particle arrays.

        Parameters
        ----------
        verbose
            Also display suffixed arrays, e.g. 'position_x',
            'position_y', etc. Default is False
        aliases
            If True, return array aliases. Default is False.

        Returns
        -------
        List
            A list of names of available arrays.
        """
        return self._available_arrays(sinks=False, verbose=verbose, aliases=aliases)

    def bulk_load(self, arrays: List[str] = None) -> Snap:
        """Load arrays into memory in bulk.

        Parameters
        ----------
        arrays
            A list of arrays to load as strings. If None, then load all
            available arrays.
        """
        if arrays is None:
            _arrays = self.available_arrays()
        else:
            _arrays = arrays
        with self.context(cache=True):
            for array in _arrays:
                try:
                    self[array]
                except ValueError as e:
                    logger.warning(f'Cannot load {array}\n{e}')

        return self

    def bulk_unload(self, arrays: List[str] = None) -> Snap:
        """Un-load arrays from memory in bulk.

        Parameters
        ----------
        arrays
            A list of arrays to load as strings. If None, then unload
            all loaded arrays.
        """
        if arrays is None:
            _arrays = self.loaded_arrays()
        else:
            _arrays = arrays
        for array in _arrays:
            try:
                del self[array]
            except KeyError:
                logger.warning(f'Cannot un-load {array}')

        return self

    @property
    def properties(self) -> Dict[str, Any]:
        """Snap properties."""
        return {key: self._properties[key] for key in sorted(self._properties.keys())}

    @property
    def sinks(self) -> Sinks:
        """Sink particle arrays."""
        return Sinks(self)

    @property
    def num_particles(self) -> int:
        """Return number of particles."""
        if self._num_particles == -1:
            self._num_particles = len(self['type'])
        return self._num_particles

    @property
    def num_particles_of_type(self) -> Dict[str, Any]:
        """Return number of particles per type."""
        if self._num_particles_of_type == -1:
            int_to_name = {idx: name for name, idx in self.particle_type.items()}
            d = {}
            ptype: Quantity = self['type']
            stype: Quantity = self['sub_type']
            for idx, num in enumerate(np.bincount(ptype.magnitude)):
                if num > 0:
                    if idx == self.particle_type['dust']:
                        # Dust particle sub-type skips zero: 1, 2, 3 ...
                        d['dust'] = list(
                            np.bincount(stype[ptype.magnitude == idx].magnitude)[1:]
                        )
                    elif idx == self.particle_type['boundary']:
                        # Boundary particle sub-type: 0 (gas), 1, 2, 3... (dust)
                        d['boundary'] = list(np.bincount(stype[ptype == idx].magnitude))
                    else:
                        d[int_to_name[idx]] = num
            self._num_particles_of_type = d
        return self._num_particles_of_type

    @property
    def num_sinks(self) -> int:
        """Return number of sinks."""
        if self._num_sinks == -1:
            try:
                self._num_sinks = len(self._sink_registry['mass'](self))
            except KeyError:
                self._num_sinks = 0
        return self._num_sinks

    @property
    def num_dust_species(self) -> int:
        """Return number of dust species."""
        if self._num_dust_species == -1:
            self._num_dust_species = len(self._properties.get('grain_size', []))
        return self._num_dust_species

    @property
    def cache_arrays(self) -> bool:
        """Cache arrays in memory for faster access."""
        return self._cache_arrays

    @cache_arrays.setter
    def cache_arrays(self, value):
        self._cache_arrays = value

    def reset(
        self, arrays: bool = False, rotation: bool = True, translation: bool = True
    ) -> Snap:
        """Reset Snap.

        Reset rotation and translations transformations on the Snap to
        initial (on-file) values. In addition, unload cached arrays.

        Parameters
        ----------
        arrays
            Set to True to unload arrays from memory. Default is False.
        rotation
            Set to True to reset rotation. Default is True.
        translation
            Set to True to reset translation. Default is True.

        Returns
        -------
        Snap
            The reset Snap. Note that the reset operation is in-place.
        """
        if any((arrays, rotation, translation)):
            for arr in self.loaded_arrays():
                del self._arrays[arr]
            for arr in self.sinks.loaded_arrays():
                del self._sinks[arr]
        else:
            logger.warning('Select something to reset')

        if rotation:
            self.rotation = None
        if translation:
            self.translation = None

        return self

    def rotate(
        self,
        rotation: Rotation = None,
        axis: Union[ndarray, List[float], Tuple[float, float, float]] = None,
        angle: float = None,
    ) -> Snap:
        """Rotate snapshot.

        The rotation can be defined by a scipy Rotation object, or a
        combination of a vector, specifying a rotation axis, and an
        angle, specifying the rotation in radians.

        Parameters
        ----------
        rotation
            The rotation as a scipy.spatial.transform.Rotation object.
        axis
            An array specifying a rotation axis, like (x, y, z).
        angle
            A float specifying the rotation in radians.

        Returns
        -------
        Snap
            The rotated Snap. Note that the rotation operation is
            in-place.

        Examples
        --------
        Rotate a Snap by π/3 around [1, 1, 0].

        >>> axis = (1, 1, 0)
        >>> angle = np.pi / 3
        >>> snap.rotate(axis=axis, angle=angle)
        """
        logger.debug(f'Rotating snapshot: {self.file_path.name}')
        if rotation is not None:
            if axis is not None or angle is not None:
                logger.warning('ignoring axis and angle as rotation is passed in')
            _rotation = rotation
        else:
            _rotation = axis / norm(axis) * angle
        if isinstance(_rotation, (list, tuple, ndarray)):
            _rotation = Rotation.from_rotvec(_rotation)
        for arr in self._vector_arrays:
            if arr in self.loaded_arrays():
                array_m, array_u = self._arrays[arr].magnitude, self._arrays[arr].units
                self._arrays[arr] = _rotation.apply(array_m) * array_u
            if arr in self.sinks.loaded_arrays():
                array_m, array_u = self._sinks[arr].magnitude, self._sinks[arr].units
                self._sinks[arr] = _rotation.apply(array_m) * array_u
        for arr in self._vector_arrays:
            if arr in self.loaded_arrays():
                del self._arrays[arr]
            if arr in self.sinks.loaded_arrays():
                del self._sinks[arr]

        if self.rotation is None:
            self.rotation = _rotation
        else:
            rot = _rotation * self.rotation
            self.rotation = rot

        return self

    def translate(
        self, translation: Union[Quantity, ndarray], unit: str = None
    ) -> Snap:
        """Translate snapshot.

        I.e. shift the snapshot origin.

        Parameters
        ----------
        translation
            The translation as an array like (x, y, z), optionally with
            units. If no units are specified you must specify the units
            parameter below.
        unit
            The length unit for the translation. E.g. 'au'.

        Returns
        -------
        Snap
            The translated Snap. Note that the translation operation is
            in-place.
        """
        logger.debug(f'Translating snapshot: {self.file_path.name}')
        if isinstance(translation, (list, tuple)):
            translation = np.array(translation)
        if translation.shape != (3,):
            raise ValueError('translation must be like (x, y, z)')
        if isinstance(translation, Quantity):
            if unit is not None:
                logger.warning('units argument ignored as translation has units')
        else:
            if unit is None:
                raise ValueError(
                    'translation must have units, or you must specify units argument'
                )
            translation *= plonk_units(unit)
        if 'position' in self.loaded_arrays():
            self._arrays['position'] += translation
        if 'position' in self.sinks.loaded_arrays():
            self._sinks['position'] += translation

        if self.translation is None:
            self.translation = translation
        else:
            self.translation += translation

        return self

    def particle_indices(
        self, particle_type: str, squeeze: bool = False,
    ) -> Union[ndarray, List[ndarray]]:
        """Particle indices of a particular type.

        Parameters
        ----------
        particle_type
            The particle type as a string.
        squeeze
            If True return all subtypes in a single array. Default is
            False.

        Returns
        -------
        ndarray or list of ndarray
            If particle has no subtypes then returns an array of
            indices of that type. Otherwise return a list of arrays of
            indices, one for each subtype. However, if squeeze
            is True, return a single array.
        """
        if particle_type == 'dust' and not squeeze:
            # Dust particle sub-type skips zero: 1, 2, 3 ...
            return [
                np.flatnonzero(
                    (self['type'] == self.particle_type['dust'])
                    & (self['sub_type'] == idx + 1)
                )
                for idx in range(self.num_dust_species)
            ]
        if particle_type == 'boundary' and not squeeze:
            # Boundary particle sub-type: 0 (gas), 1, 2, 3... (dust)
            return [
                np.flatnonzero(
                    (self['type'] == self.particle_type['boundary'])
                    & (self['sub_type'] == idx)
                )
                for idx in range(self.num_dust_species + 1)
            ]
        return np.flatnonzero(self['type'] == self.particle_type[particle_type])

    def subsnaps_by_type(
        self, split_subtypes: bool = True
    ) -> Dict[str, Union[SubSnap, List[SubSnap]]]:
        """Return particle-type subsnaps as a dict of SubSnaps.

        Parameters
        ----------
        split_subtypes : optional
            If True, split particle types into subtypes.

        Returns
        -------
        Dict
            A dict of all SubSnaps.
        """
        if not split_subtypes:
            raise NotImplementedError('Combined sub-types not implemented yet')
        subsnaps_by_type: Dict[str, Any] = dict()
        for key in self.particle_type:
            try:
                subsnap = self[key]
                if isinstance(subsnap, list):
                    if len(subsnap[0]) > 0:
                        subsnaps_by_type[key] = self[key]
                else:
                    if len(subsnap) > 0:
                        subsnaps_by_type[key] = self[key]
            except ValueError:
                pass

        return subsnaps_by_type

    def subsnaps_as_list(self, split_subtypes: bool = True) -> List[SubSnap]:
        """Return particle-type subsnaps as a list of SubSnaps.

        Parameters
        ----------
        split_subtypes : optional
            If True, split particle types into subtypes.

        Returns
        -------
        List[SubSnap]
            A list of all SubSnaps.
        """
        if not split_subtypes:
            raise NotImplementedError('Combined sub-types not implemented yet')
        subsnaps: List[SubSnap] = list()
        for val in self.subsnaps_by_type().values():
            if isinstance(val, list):
                subsnaps = subsnaps + val
            else:
                subsnaps.append(val)
        return subsnaps

    def set_kernel(self, kernel: str) -> Snap:
        """Set kernel.

        Parameters
        ----------
        kernel
            The kernel name as a string.

        Returns
        -------
        Snap
            The Snap.
        """
        if kernel not in kernel_names:
            raise ValueError(f'Kernel must be in {kernel_names}')
        self._properties['kernel'] = kernel

        return self

    def set_gravitational_parameter(self, sink_idx: Union[int, List[int]]) -> Snap:
        """Set standard gravitational parameter.

        Calculate the standard gravitational parameter (G M) given a
        sink index, or list of sink indices for multiple systems, etc,
        and set snap.properties['gravitational_parameter']. This is
        required to calculate orbital quantities such as eccentricity
        or Stokes number.

        Parameters
        ----------
        sink_idx
            The sink index or list of indices.

        Returns
        -------
        Snap
            The Snap.
        """
        G = plonk_units.newtonian_constant_of_gravitation
        if isinstance(sink_idx, (int, list)):
            M = self.sinks[sink_idx]['mass']
        else:
            raise ValueError('Cannot determine gravitational parameter')
        G = (G * M).to_base_units()
        self._properties['gravitational_parameter'] = G

        return self

    def set_molecular_weight(self, molecular_weight: float) -> Snap:
        """Set molecular weight.

        Set the molecular weight of the gas in gram / mole in
        snap.properties['molecular_weight']. This is required to
        calculate temperature.

        Parameters
        ----------
        molecular_weight
            The molecular weight in units of gram / mole. E.g. Phantom
            uses 2.381 for molecular hydrogen with solar metallicity.

        Returns
        -------
        Snap
            The Snap.
        """
        self._properties['molecular_weight'] = molecular_weight
        return self

    @property
    def tree(self) -> cKDTree:
        """Particle neighbour kd-tree.

        Trees are represented by scipy cKDTree objects.
        """
        if self._tree is None:
            pos: Quantity = self['position']
            self._tree = cKDTree(pos.magnitude)
        return self._tree

    def neighbours(self, indices: Union[ndarray, List[int]]) -> ndarray:
        """Get neighbours of particles.

        Parameters
        ----------
        indices
            The particle indices.

        Returns
        -------
        ndarray of list
            An array of neighbours lists.
        """
        kernel = self._properties.get('kernel')
        if kernel is None:
            raise ValueError(
                'To calculate particle neighbours, first set the kernel\n'
                'via snap.set_kernel.'
            )

        subsnap = self[indices]
        position: Quantity = subsnap['position']
        smoothing_length: Quantity = subsnap['smoothing_length']
        r_kern = kernel_radius[kernel]

        neighbours = self.tree.query_ball_point(
            position.magnitude, r_kern * smoothing_length.magnitude, n_jobs=-1,
        )

        return neighbours

    def write_extra_arrays(
        self, arrays: List[str], filename: Union[str, Path] = None
    ) -> Snap:
        """Write extra arrays to file.

        Parameters
        ----------
        arrays
            A list of strings with array names.
        filename : optional
            A filename to write to.

        Returns
        -------
        Snap
            The Snap.
        """
        if filename is None:
            filename = f'{self.file_path.stem}_extra.h5'
        f = h5py.File(filename, mode='w')
        for array in arrays:
            arr_with_units: Quantity = self[array]
            units = self.array_code_unit(array)
            arr = (arr_with_units / units).magnitude
            dset = f.create_dataset(
                array, arr.shape, dtype=arr.dtype, compression='gzip',
            )
            dset[:] = arr
        f.close()

        return self

    def read_extra_arrays(self, filename: Union[str, Path] = None) -> Snap:
        """Read extra arrays from file.

        Parameters
        ----------
        filename : optional
            A filename to read from.

        Returns
        -------
        Snap
            The Snap.
        """
        if filename is None:
            filename = f'{self.file_path.stem}_extra.h5'
        f = h5py.File(filename, mode='r')
        for array in f:
            self[array] = f[array][()] * plonk_units['dimensionless']
        f.close()

        return self

    def to_dataframe(
        self, columns: Union[Tuple[str, ...], List[str]] = None, units: List[str] = None
    ) -> DataFrame:
        """Convert Snap to DataFrame.

        Parameters
        ----------
        columns : optional
            A list of columns to add to the data frame. Default is
            None.
        units : optional
            A list of units corresponding to columns add to the data
            frame. Units must be strings, and must be base units. I.e.
            'cm' not '10 cm'. If None, use default, i.e. cgs. Default is
            None.

        Returns
        -------
        DataFrame
        """
        d = dict()
        if columns is None:
            columns = self.loaded_arrays()
        cols = list(columns)
        if units is None:
            _units = list()
            for column in cols:
                try:
                    arr: Quantity = self[column]
                    _units.append(arr.units)
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
        if len(_units) != len(cols):
            raise ValueError('units and columns must have same length')
        for column, unit in zip(cols, _units):
            name = column
            array: Quantity = self[column]
            try:
                suffix = f' [{unit:~}]'
                array = array.to(unit).magnitude
            except AttributeError:
                suffix = ''
            if array.ndim == 1:
                d[name + suffix] = array
            if array.ndim == 2:
                for idx in range(array.shape[1]):
                    d[f'{name}.{idx+1}' + suffix] = array[:, idx]
        return pd.DataFrame(d)

    def family(self, name: str, squeeze: bool = False) -> Union[SubSnap, List[SubSnap]]:
        """Get a SubSnap of a particle family by name.

        Parameters
        ----------
        name
            A string representing the name of the family of particles,
            e.g. 'gas'.
        squeeze
            Squeeze sub-types. If False and the particle family has
            sub-types then return a list of SubSnaps of each sub-type.
            Otherwise return a SubSnap with all particle of that type.

        Returns
        -------
        SubSnap or List[SubSnap]
            If the particle family has no sub-types then return a
            SubSnap. Otherwise, if split_subtype is True, return a list
            of SubSnaps.
        """
        if name in self.particle_type:
            ind = self.particle_indices(particle_type=name, squeeze=squeeze)
            if isinstance(ind, list):
                return [SubSnap(self, _ind) for _ind in ind]
            return SubSnap(self, ind)
        raise ValueError('Family not available')

    def array(self, name: str, sinks: bool = False) -> Quantity:
        """Get a particle (or sink) array.

        Parameters
        ----------
        name
            A string representing the name of the particle array.
        sinks
            Whether or not to reference particle or sinks arrays.

        Returns
        -------
        Quantity
        """
        base_name = self.base_array_name(name)
        suffix = self._array_suffix(name)
        array = self._get_array(base_name, sinks)
        if suffix == '':
            return array
        transform, _slice, kwargs = self._array_transform(
            base_name=base_name, suffix=suffix
        )
        return transform(array, **kwargs)[_slice]

    def array_code_unit(self, arr: str) -> Quantity:
        """Get array code units.

        Parameters
        ----------
        arr
            The string representing the quantity.

        Returns
        -------
        Quantity
            The Pint unit quantity, or the float 1.0 if no unit found.
        """
        base_name = self.base_array_name(arr)
        try:
            unit = self._array_code_units[base_name]
        except KeyError:
            _arr: Quantity = self[base_name]
            dim = _arr.units.dimensionality
            unit = 1.0
            for d in self.code_units:
                unit *= self.code_units[d] ** dim[f'[{d}]']
        return unit

    def array_in_code_units(self, name: str) -> ndarray:
        """Get array in code units.

        Parameters
        ----------
        snap
            The Snap or SubSnap.
        name
            The array name.

        Returns
        -------
        ndarray
            The array on the particles in code units.
        """
        return (self[name] / self.array_code_unit(name)).magnitude

    def base_array_name(self, name: str) -> str:
        """Get the base array name from a string.

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
            The base array name.
        """
        if name in self.available_arrays():
            return name
        if name in self.sinks.available_arrays():
            return name
        if name in self._array_aliases:
            return self._array_aliases[name]
        if name in self._arrays:
            return name
        if name in self._sinks:
            return name

        name_root = '_'.join(name.split('_')[:-1])
        name_suffix = name.split('_')[-1]
        if name_root == '' and name_suffix in ('x', 'y', 'z'):
            return 'position'
        if name_root in self._array_aliases:
            return self._array_aliases[name_root]
        if name_root in self._vector_arrays and name_suffix in ('x', 'y', 'z', 'mag'):
            return name_root
        if name_root in self._dust_arrays:
            if _str_is_int(name_suffix) or name_suffix == 'tot':
                return name_root

        raise ValueError('Unknown array')

    def _array_suffix(self, name: str) -> str:
        if (
            name in self.available_arrays()
            or name in self.sinks.available_arrays()
            or name in self._array_aliases
        ):
            return ''
        if name in self._arrays:
            return ''
        if name in self._sinks:
            return ''

        name_root = '_'.join(name.split('_')[:-1])
        name_suffix = name.split('_')[-1]
        if name_root == '' and name_suffix in ('x', 'y', 'z'):
            return name_suffix
        name_root = self.base_array_name(name_root)
        if name_root in self._vector_arrays:
            if name_suffix in ('x', 'y', 'z', 'mag'):
                return name_suffix
        if name_root in self._dust_arrays:
            if _str_is_int(name_suffix):
                return name_suffix
            if name_suffix == 'tot':
                return name_suffix

        raise ValueError('Unknown array')

    def _array_transform(
        self, base_name: str, suffix: str
    ) -> Tuple[Callable, tuple, Dict[str, Any]]:
        def nothing(x):
            return x

        if base_name in self._vector_arrays:
            if suffix == 'x':
                return nothing, (..., 0), {}
            if suffix == 'y':
                return nothing, (..., 1), {}
            if suffix == 'z':
                return nothing, (..., 2), {}
            if suffix == 'mag':
                return norm, (), {'axis': 1}
        if base_name in self._dust_arrays:
            if _str_is_int(suffix):
                return nothing, (..., int(suffix) - 1), {}
            if suffix == 'tot':
                return np.sum, (), {'axis': 1}

        raise ValueError('Unknown array')

    def _get_array_from_registry(self, name: str, sinks: bool = False) -> Quantity:
        if sinks:
            array = self._sink_registry[name](self)
        else:
            array = self._array_registry[name](self)
        if self.rotation is not None and name in self._vector_arrays:
            array_m, array_u = array.magnitude, array.units
            array = self.rotation.apply(array_m) * array_u
        if self.translation is not None and name == 'position':
            array += self.translation
        return array

    def _get_array(self, name: str, sinks: bool = False) -> Quantity:
        """Get an array by name."""
        if sinks:
            array_dict = self._sinks
        else:
            array_dict = self._arrays
        if name in array_dict:
            return array_dict[name]
        if name in self._array_registry or name in self._sink_registry:
            array = self._get_array_from_registry(name, sinks)
            if self.cache_arrays:
                if sinks:
                    self._sinks[name] = array
                else:
                    self._arrays[name] = array
            return array
        raise ValueError('Array not available')

    def _getitem(
        self, inp: Union[str, ndarray, int, slice], sinks: bool = False,
    ) -> Union[Quantity, SubSnap, List[SubSnap]]:
        """Return an array, or family, or subset."""
        if isinstance(inp, str):
            if inp in self.particle_type:
                return self.family(name=inp)
            return self.array(name=inp, sinks=sinks)
        if sinks:
            raise ValueError('Cannot return sinks as SubSnap')
        return SubSnap(self, inp)

    def __getitem__(self, inp: Union[str, ndarray, int, slice]):
        """Return an array, or family, or subset."""
        return self._getitem(inp, sinks=False)

    def __setitem__(self, name: str, item: Quantity):
        """Set a particle array."""
        if not isinstance(item, Quantity):
            raise ValueError('"item" must be Pint Quantity')
        if item.shape[0] != len(self):
            raise ValueError('Length of array does not match particle number')
        if name in self.loaded_arrays():
            raise ValueError(
                'Attempting to overwrite existing array. To do so, first delete the '
                'array\nwith del snap["array"], then try again.'
            )
        if name in self._available_arrays() or name in self._array_aliases.keys():
            raise ValueError(
                'Attempting to set array already available. '
                'See snap.available_arrays().'
            )
        self._arrays[name] = item

    def __delitem__(self, name):
        """Delete an array from memory."""
        del self._arrays[name]

    def _ipython_key_completions_(self) -> List[str]:
        """Tab completion for IPython __getitem__ method."""
        return self.available_arrays(verbose=True)

    def __len__(self) -> int:
        """Length as number of particles."""
        return self.num_particles

    def __repr__(self) -> str:
        """Dunder repr method."""
        return self.__str__()

    def __str__(self) -> str:
        """Dunder str method."""
        return f'<plonk.Snap "{self.file_path.name}">'

    image = visualize.image
    plot = visualize.plot
    vector = visualize.vector
    context = context.context


class SubSnap(Snap):
    """A Snap subset of particles.

    A SubSnap can be generated from a Snap via an index array, a
    particle mask, or a string. SubSnaps can be used like a Snap,
    including: accessing arrays, plotting, finding neighbours, etc.

    Parameters
    ----------
    base
        The base Snap or SubSnap
    indices
        A (N,) array of particle indices to include in the SubSnap.

    Examples
    --------
    Generate a SubSnap directly.

    >>> subsnap = SubSnap(base=snap, indices=[0, 1, 2, 3])

    You can generate a SubSnap from a Snap object. For example, generate
    a SubSnap of the gas particles on a Snap.

    >>> subsnap = snap['gas']

    Generate a SubSnap of particles with a mask.

    >>> subsnap = snap[snap['x'] > 0]

    Generate a SubSnap of particles from indices.

    >>> subsnap = snap[:100]
    >>> subsnap = snap[[0, 9, 99]]
    """

    def __init__(self, base: Snap, indices: Union[ndarray, slice, list, int, tuple]):
        super().__init__()

        self.base = base

        ind = _input_indices_array(inp=indices, max_slice=len(base))
        if len(ind) == 0:
            raise ValueError('SubSnap has no particles')
        self._indices = ind

        # Attributes different to Snap
        self._num_particles = len(self._indices)
        self._num_particles_of_type = -1
        self._num_dust_species = -1
        self._tree = None

        # Attributes same as Snap
        self.data_source = self.base.data_source
        self.file_path = self.base.file_path
        self.code_units = self.base.code_units
        self._properties = self.base._properties
        self._array_code_units = self.base._array_code_units
        self._array_registry = self.base._array_registry
        self._sink_registry = self.base._sink_registry
        self._cache_arrays = self.base._cache_arrays
        self._arrays = self.base._arrays
        self._sinks = self.base._sinks
        self._file_pointer = self.base._file_pointer
        self.rotation = self.base.rotation
        self.translation = self.base.translation

    @property
    def indices(self) -> ndarray:
        """Particle indices."""
        return self._indices

    @property
    def sinks(self) -> Sinks:
        """Sink particle arrays."""
        return self.base.sinks

    def __repr__(self) -> str:
        """Dunder repr method."""
        return self.__str__()

    def __str__(self) -> str:
        """Dunder str method."""
        return f'<plonk.SubSnap "{self.file_path.name}">'

    def _get_array(self, name: str, sinks: bool = False) -> Quantity:
        return self.base._get_array(name, sinks)[self.indices]


SnapLike = Union[Snap, SubSnap]


class Sinks:
    """Sink particles in a Snap.

    A Sinks object is generated from a Snap.

    Parameters
    ----------
    base
        The base Snap.
    indices : optional
        Indices to specify a subset of sink particles.

    Examples
    --------
    Generate a Sinks object directly.

    >>> sinks = Sinks(base=snap)

    Generate a Sinks object from a Snap object.

    >>> sinks = snap.sinks

    Choose a subset of sink particles.

    >>> sinks = snap.sinks[[0, 1]]
    >>> star = snap.sinks[0]
    >>> planets = snap.sinks[1:4]
    """

    def __init__(
        self, base: Snap, indices: Union[ndarray, slice, list, int, tuple] = None
    ):
        self.base = base

        if indices is None:
            indices = np.arange(base.num_sinks)
        ind = _input_indices_array(inp=indices, max_slice=base.num_sinks)
        if len(ind) == 0:
            raise ValueError('Sinks has no particles')
        self._indices = ind

        # Attributes same as Snap
        self.file_path = self.base.file_path

    @property
    def indices(self) -> ndarray:
        """Sink particle indices."""
        return self._indices

    def available_arrays(
        self, verbose: bool = False, aliases: bool = False
    ) -> List[str]:
        """Return a list of available sink arrays.

        Parameters
        ----------
        verbose
            Also display suffixed arrays, e.g. 'position_x',
            'position_y', etc. Default is False
        aliases
            If True, return array aliases. Default is False.

        Returns
        -------
        List
            A list of names of available arrays.
        """
        return self.base._available_arrays(sinks=True, verbose=verbose, aliases=aliases)

    def loaded_arrays(self) -> List[str]:
        """Return a list of loaded arrays.

        Returns
        -------
        List
            A list of names of loaded sink arrays.
        """
        return sorted(self.base._sinks.keys())

    def array(self, name: str) -> Quantity:
        """Get an array.

        Parameters
        ----------
        name
            A string representing the name of the particle array.

        Returns
        -------
        Quantity
        """
        return np.squeeze(self.base.array(name, sinks=True)[self.indices])[()]

    def __setitem__(self, name: str, item: Quantity):
        """Set an array."""
        if not isinstance(item, Quantity):
            raise ValueError('"item" must be an array with units, i.e. a Pint Quantity')
        if item.shape[0] != len(self):
            raise ValueError('Length of array does not match particle number')
        self.base._sinks[name] = item

    def __getitem__(self, inp):
        """Return an array or subset."""
        if isinstance(inp, str):
            return self.array(name=inp)
        ind = _input_indices_array(inp=inp, max_slice=len(self))
        if ind is not None:
            return Sinks(self.base, ind)
        raise ValueError('Cannot determine item to return')

    def __len__(self) -> int:
        """Length as number of particles."""
        return len(self.indices)

    def __repr__(self) -> str:
        """Dunder repr method."""
        return self.__str__()

    def __str__(self) -> str:
        """Dunder str method."""
        return f'<plonk.Sinks "{self.base.file_path.name}">'

    def _ipython_key_completions_(self) -> List[str]:
        """Tab completion for IPython __getitem__ method."""
        return self.available_arrays(verbose=True)

    plot = visualize.plot


def _str_is_int(string: str) -> bool:
    try:
        int(string)
        return True
    except ValueError:
        return False


def _input_indices_array(
    inp: Union[ndarray, slice, list, int, tuple], max_slice: int
) -> Union[ndarray, List[int]]:
    """Take array, slice, int, list, tuple and return indices array."""
    if isinstance(inp, ndarray):
        if np.issubdtype(np.bool, inp.dtype):
            return np.flatnonzero(inp)
        if np.issubdtype(np.int, inp.dtype):
            return inp
    if isinstance(inp, (list, tuple)):
        if isinstance(inp[0], int):
            return np.array(inp)
    if isinstance(inp, int):
        return np.array([inp])
    if isinstance(inp, slice):
        i1 = inp.start if inp.start is not None else 0
        i2 = inp.stop if inp.stop is not None else max_slice
        return np.arange(i1, i2, inp.step)
    return []
