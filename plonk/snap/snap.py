"""Snap, SubSnap, Sinks classes for snapshot files.

The Snap class contains all information related to a smoothed particle
hydrodynamics simulation snapshot file. The SubSnap class is for
accessing a subset of particles in a Snap.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Callable, Dict, List, Set, Tuple, Union

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

    _array_name_mapper = {
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
        'vel': 'velocity',
        'v_R': 'velocity_radial_cylindrical',
        'v_r': 'velocity_radial_spherical',
        'v_phi': 'angular_velocity',
        'velocity_R': 'radial_velocity_cylindrical',
        'velocity_r': 'radial_velocity_spherical',
        'velocity_phi': 'angular_velocity',
        'xyz': 'position',
    }

    _array_split_mapper = {
        'x': ('position', 0),
        'y': ('position', 1),
        'z': ('position', 2),
        'v_x': ('velocity', 0),
        'v_y': ('velocity', 1),
        'v_z': ('velocity', 2),
        'p_x': ('momentum', 0),
        'p_y': ('momentum', 1),
        'p_z': ('momentum', 2),
        'L_x': ('angular_momentum', 0),
        'L_y': ('angular_momentum', 1),
        'L_z': ('angular_momentum', 2),
        'j_x': ('specific_angular_momentum', 0),
        'j_y': ('specific_angular_momentum', 1),
        'j_z': ('specific_angular_momentum', 2),
        's_x': ('spin', 0),
        's_y': ('spin', 1),
        's_z': ('spin', 2),
        'B_x': ('magnetic_field', 0),
        'B_y': ('magnetic_field', 1),
        'B_z': ('magnetic_field', 2),
    }

    _vector_arrays = {
        'magnetic_field',
        'position',
        'spin',
        'velocity',
    }

    _vector_component_arrays: Set[str] = set()

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
        Snap._array_name_mapper[alias] = name

    def __init__(self):

        self.data_source = None
        self.file_path = None
        self.code_units = {}
        self._properties = {}
        self._array_units = {}
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
        self._families = {key: None for key in self.particle_type}
        self.rotation = None
        self.translation = None
        self._tree = None

    def close_file(self):
        """Close access to underlying file."""
        self._file_pointer.close()

    def reopen_file(self):
        """Re-open access to the underlying file."""
        self._file_pointer = h5py.File(self.file_path, mode='r')

    def add_array(self, rotatable: bool = None, dust: bool = False) -> Callable:
        """Decorate function to add array to Snap.

        This function decorates a function that returns an array. The
        name of the function is the string with which to reference the
        array.

        Parameters
        ----------
        rotatable
            A bool to represent if the array should have rotations
            applied to it, i.e. it is a vector arrray. If True the
            rotation should be applied. If False the rotation cannot be
            applied. If None no rotation is required. Default is None.
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
            if rotatable is True:
                self._vector_arrays.add(fn.__name__)
            elif rotatable is False:
                self._vector_component_arrays.add(fn.__name__)
            if dust is True:
                self._dust_arrays.add(fn.__name__)
            return fn

        return _add_array

    def loaded_arrays(self):
        """Return a tuple of loaded arrays.

        Returns
        -------
        A tuple of names of loaded particle arrays.
        """
        return tuple(sorted(self._arrays.keys()))

    def _available_arrays(
        self, sinks: bool = False, verbose: bool = False, aliases: bool = False
    ):
        """Return a tuple of available arrays.

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
        A tuple of names of available arrays.
        """
        if sinks:
            loaded = self.sinks.loaded_arrays()
            registered = list(self._sink_registry.keys())
        else:
            loaded = self.loaded_arrays()
            registered = list(sorted(self._array_registry.keys()))
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
            extra = tuple(
                key
                for key, val in self._array_split_mapper.items()
                if val[0] in self.loaded_arrays() or val[0] in self._array_registry
            )
            extra += tuple(
                key
                for key, val in self._array_name_mapper.items()
                if val[0] in self.loaded_arrays() or val in self._array_registry
            )
            return tuple(sorted(set(extra), key=lambda x: x.lower()))

        return tuple(sorted(set(loaded + tuple(registered))))

    def available_arrays(self, verbose: bool = False, aliases: bool = False):
        """Return a tuple of available particle arrays.

        Parameters
        ----------
        verbose
            Also display suffixed arrays, e.g. 'position_x',
            'position_y', etc. Default is False
        aliases
            If True, return array aliases. Default is False.

        Returns
        -------
        A tuple of names of available arrays.
        """
        return self._available_arrays(sinks=False, verbose=verbose, aliases=aliases)

    def bulk_load(self, arrays: List[str] = None):
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
        for array in _arrays:
            try:
                self[array]
            except ValueError as e:
                logger.warning(f'Cannot load {array}\n{e}')

        return self

    @property
    def properties(self):
        """Snap properties."""
        return {key: self._properties[key] for key in sorted(self._properties.keys())}

    @property
    def sinks(self):
        """Sink particle arrays."""
        return Sinks(self)

    @property
    def num_particles(self):
        """Return number of particles."""
        if self._num_particles == -1:
            self._num_particles = len(self['type'])
        return self._num_particles

    @property
    def num_particles_of_type(self):
        """Return number of particles per type."""
        if self._num_particles_of_type == -1:
            int_to_name = {idx: name for name, idx in self.particle_type.items()}
            d = {}
            for idx, num in enumerate(np.bincount(self['type'].magnitude)):
                if num > 0:
                    if idx == self.particle_type['dust']:
                        # Dust particle sub-type skips zero: 1, 2, 3 ...
                        d['dust'] = list(
                            np.bincount(
                                self[self['type'] == idx]['sub_type'].magnitude
                            )[1:]
                        )
                    elif idx == self.particle_type['boundary']:
                        # Boundary particle sub-type: 0 (gas), 1, 2, 3... (dust)
                        d['boundary'] = list(
                            np.bincount(self[self['type'] == idx]['sub_type'].magnitude)
                        )
                    else:
                        d[int_to_name[idx]] = num
            self._num_particles_of_type = d
        return self._num_particles_of_type

    @property
    def num_sinks(self):
        """Return number of sinks."""
        if self._num_sinks == -1:
            try:
                self._num_sinks = len(self._sink_registry['mass'](self))
            except KeyError:
                self._num_sinks = 0
        return self._num_sinks

    @property
    def num_dust_species(self):
        """Return number of dust species."""
        if self._num_dust_species == -1:
            self._num_dust_species = len(self._properties.get('grain_size', []))
        return self._num_dust_species

    @property
    def cache_arrays(self):
        """Cache arrays in memory for faster access."""
        return self._cache_arrays

    @cache_arrays.setter
    def cache_arrays(self, value):
        if value is False:
            self._arrays = {}
        self._cache_arrays = value

    def unset(self, rotation: bool = True, translation: bool = True):
        """Unset transformations.

        Unset rotation and translations transformations on the Snap.

        Parameters
        ----------
        rotation
            Set to True to unset rotation. Default is True.
        translation
            Set to True to unset translation. Default is True.
        """
        if any((rotation, translation)):
            for arr in self.loaded_arrays():
                del self._arrays[arr]
            for arr in self.sinks.loaded_arrays():
                del self._sinks[arr]
        else:
            raise ValueError('Select something to unset')

        if rotation:
            self.rotation = None
        if translation:
            self.translation = None

        return self

    def rotate(
        self,
        rotation: Rotation = None,
        axis: Union[ndarray, List, Tuple[float, float, float]] = None,
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
        for arr in self._vector_component_arrays:
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
        self, particle_type: str, squeeze_subtype: bool = False,
    ) -> Union[ndarray, List[ndarray]]:
        """Particle indices of a particular type.

        Parameters
        ----------
        particle_type
            The particle type as a string.
        squeeze_subtype
            If True return all subtypes in a single array. Default is
            False.

        Returns
        -------
        ndarray or list of ndarray
            If particle has no subtypes then returns an array of
            indices of that type. Otherwise return a list of arrays of
            indices, one for each subtype. However, if squeeze_subtype
            is True, return a single array.
        """
        if particle_type == 'dust' and not squeeze_subtype:
            # Dust particle sub-type skips zero: 1, 2, 3 ...
            return [
                np.flatnonzero(
                    (self['type'] == self.particle_type['dust'])
                    & (self['sub_type'] == idx + 1)
                )
                for idx in range(self.num_dust_species)
            ]
        if particle_type == 'boundary' and not squeeze_subtype:
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

    def set_kernel(self, kernel: str):
        """Set kernel.

        Parameters
        ----------
        kernel
            The kernel name as a string.
        """
        if kernel not in kernel_names:
            raise ValueError(f'Kernel must be in {kernel_names}')
        self._properties['kernel'] = kernel

    def set_gravitational_parameter(self, sink_idx: Union[int, List[int]]):
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
        """
        G = plonk_units.newtonian_constant_of_gravitation
        if isinstance(sink_idx, (int, list)):
            M = self.sinks[sink_idx]['mass']
        else:
            raise ValueError('Cannot determine gravitational parameter')
        G = (G * M).to_base_units()
        self._properties['gravitational_parameter'] = G

    def set_molecular_weight(self, molecular_weight: float):
        """Set molecular weight.

        Set the molecular weight of the gas in gram / mole in
        snap.properties['molecular_weight']. This is required to
        calculate temperature.

        Parameters
        ----------
        molecular_weight
            The molecular weight in units of gram / mole. E.g. Phantom
            uses 2.381 for molecular hydrogen with solar metallicity.
        """
        self._properties['molecular_weight'] = molecular_weight

    @property
    def tree(self):
        """Particle neighbour kd-tree.

        Trees are represented by scipy cKDTree objects.
        """
        if self._tree is None:
            self._tree = cKDTree(self['position'].magnitude)
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

    def write_extra_arrays(self, arrays: List[str], filename: Union[str, Path] = None):
        """Write extra arrays to file.

        Parameters
        ----------
        arrays
            A list of strings with array names.
        filename : optional
            A filename to write to.
        """
        if filename is None:
            filename = f'{self.file_path.stem}_extra.h5'
        f = h5py.File(filename, mode='w')
        for array in arrays:
            arr_with_units: Quantity = self[array]
            units = self.get_array_code_unit(array)
            arr = (arr_with_units / units).magnitude
            dset = f.create_dataset(
                array, arr.shape, dtype=arr.dtype, compression='gzip',
            )
            dset[:] = arr
        f.close()

    def read_extra_arrays(self, filename: Union[str, Path] = None):
        """Read extra arrays from file.

        Parameters
        ----------
        filename : optional
            A filename to read from.
        """
        if filename is None:
            filename = f'{self.file_path.stem}_extra.h5'
        f = h5py.File(filename, mode='r')
        for array in f:
            self[array] = f[array][()]
        f.close()

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

    def _get_family_subsnap(self, name: str):
        """Get a family by name."""
        if name in self._families:
            if self._families[name] is None:
                ind = self.particle_indices(name)
                if len(ind) == 0:
                    raise ValueError(f'No {name} particles available')
                self._families[name] = ind
            ind = self._families[name]
            if isinstance(ind, list):
                return [SubSnap(self, _ind) for _ind in ind]
            return SubSnap(self, ind)
        raise ValueError('Family not available')

    def get_array_code_unit(self, arr: str) -> Any:
        """Get array code units.

        Parameters
        ----------
        arr
            The string representing the quantity.

        Returns
        -------
        unit
            The Pint unit quantity, or the float 1.0 if no unit found.
        """
        arr_root = '_'.join(arr.split('_')[:-1])
        if arr in self._array_split_mapper:
            arr = self._array_split_mapper[arr][0]
        elif arr in self._array_name_mapper:
            arr = self._array_name_mapper[arr]
        elif arr_root in self._vector_arrays | self._dust_arrays:
            arr = arr_root
        try:
            unit = self._array_units[arr]
        except KeyError:
            _arr: Quantity = self[arr]
            dim = _arr.units.dimensionality
            unit = 1.0
            for d in ['length', 'mass', 'time']:
                unit *= self.code_units[d] ** dim[f'[{d}]']
        return unit

    def get_canonical_array_name(self, name: str) -> str:
        """Get the canonical array name from a string.

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
            The canonical name.
        """
        if name in self.available_arrays():
            return name
        if name in self._array_name_mapper:
            return self._array_name_mapper[name]
        if name in self._array_split_mapper:
            return self._array_split_mapper[name][0]
        if name in self._arrays:
            return name

        name_root = '_'.join(name.split('_')[:-1])
        if name_root not in self.available_arrays():
            raise ValueError('Unknown array')
        name_suffix = name.split('_')[-1]
        if name_root in self._vector_arrays:
            if name_suffix in ('x', 'y', 'z', 'mag'):
                return name_root
        if name_root in self._dust_arrays:
            if _str_is_int(name_suffix):
                return name_root
            if name_suffix == 'tot':
                return name_root

        raise ValueError('Unknown array')

    def _get_array_from_registry(self, name: str, sinks: bool = False):
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

    def _get_array(self, name: str, sinks: bool = False) -> ndarray:
        """Get an array by name."""
        index = None

        if name in self._available_arrays(sinks):
            pass
        elif name in self._array_name_mapper.keys():
            name = self._array_name_mapper[name]
        elif name in self._array_split_mapper.keys():
            name, index = self._array_split_mapper[name]
        else:
            raise ValueError('Array not available')

        if sinks:
            array_dict = self._sinks
        else:
            array_dict = self._arrays
        if name in array_dict:
            if index is None:
                return array_dict[name]
            return array_dict[name][:, index]
        if name in self._array_registry or name in self._sink_registry:
            array = self._get_array_from_registry(name, sinks)
            if self.cache_arrays:
                if sinks:
                    self._sinks[name] = array
                else:
                    self._arrays[name] = array
            if index is None:
                return array
            return array[:, index]
        raise ValueError('Array not available')

    def _getitem_from_str(self, inp: str, sinks: bool = False) -> ndarray:
        """Return item from string."""
        inp_root = '_'.join(inp.split('_')[:-1])
        inp_suffix = inp.split('_')[-1]

        if inp in self._families:
            return self._get_family_subsnap(inp)
        if inp in self._available_arrays(sinks):
            return self._get_array(inp, sinks)
        if inp in self._array_name_mapper.keys():
            return self._get_array(inp, sinks)
        if inp in self._array_split_mapper.keys():
            return self._get_array(inp, sinks)
        if inp_root in self._vector_arrays:
            if inp_suffix == 'x':
                return self._get_array(inp_root, sinks)[:, 0]
            if inp_suffix == 'y':
                return self._get_array(inp_root, sinks)[:, 1]
            if inp_suffix == 'z':
                return self._get_array(inp_root, sinks)[:, 2]
            if inp_suffix == 'mag':
                return norm(self._get_array(inp_root, sinks), axis=1)
        if inp_root in self._dust_arrays:
            if _str_is_int(inp_suffix):
                return self._get_array(inp_root)[:, int(inp_suffix) - 1]
            if inp_suffix == 'tot':
                return self._get_array(inp_root).sum(axis=1)

        raise ValueError('Cannot determine item to return.')

    def _getitem(
        self, inp: Union[str, ndarray, int, slice], sinks: bool = False,
    ) -> Union[ndarray, SubSnap]:
        """Return an array, or family, or subset."""
        if isinstance(inp, str):
            return self._getitem_from_str(inp, sinks)
        if sinks:
            raise ValueError('Cannot return sinks as SubSnap')
        ind = _input_indices_array(inp=inp, max_slice=len(self))
        if ind is not None:
            return SubSnap(self, ind)
        raise ValueError('Cannot determine item to return')

    def __getitem__(
        self, inp: Union[str, ndarray, int, slice]
    ) -> Union[ndarray, SubSnap]:
        """Return an array, or family, or subset."""
        return self._getitem(inp, sinks=False)

    def __setitem__(self, name: str, item: ndarray):
        """Set an array."""
        if not isinstance(item, (ndarray, Quantity)):
            raise ValueError('"item" must be ndarray or Pint Quantity')
        if item.shape[0] != len(self):
            raise ValueError('Length of array does not match particle number')
        if name in self.loaded_arrays():
            raise ValueError(
                'Attempting to overwrite existing array. To do so, first delete the '
                'array\nwith del snap["array"], then try again.'
            )
        if (
            name in self._available_arrays()
            or name in self._array_split_mapper.keys()
            or name in self._array_name_mapper.keys()
        ):
            raise ValueError(
                'Attempting to set array already available. '
                'See snap.available_arrays().'
            )
        self._arrays[name] = item

    def __delitem__(self, name):
        """Delete an array from memory."""
        del self._arrays[name]

    def _ipython_key_completions_(self):
        """Tab completion for IPython __getitem__ method."""
        return self.available_arrays(verbose=True)

    def __len__(self):
        """Length as number of particles."""
        return self.num_particles

    def __repr__(self):
        """Dunder repr method."""
        return self.__str__()

    def __str__(self):
        """Dunder str method."""
        return f'<plonk.Snap "{self.file_path.name}">'

    image = visualize.image
    plot = visualize.plot
    vector = visualize.vector


class SubSnap(Snap):
    """A Snap subset of particles.

    A SubSnap can be generated from a Snap via an index array, a
    particle mask, or a string. SubSnaps can be used like a Snap,
    including: accessing arrays, plotting, finding neighbours, etc.

    Parameters
    ----------
    base
        The base Snap.
    indices
        A (N,) array of particle indices to include in the SubSnap.

    Examples
    --------
    Generate a SubSnap directly.

    >>> subsnap = SubSnap(snap=snap, indices=[0, 1, 2, 3])

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
        if ind is None:
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
        self._array_units = self.base._array_units
        self._array_registry = self.base._array_registry
        self._sink_registry = self.base._sink_registry
        self._cache_arrays = self.base._cache_arrays
        self._arrays = self.base._arrays
        self._sinks = self.base._sinks
        self._file_pointer = self.base._file_pointer
        self.rotation = self.base.rotation
        self.translation = self.base.translation

    @property
    def indices(self):
        """Particle indices."""
        return self._indices

    @property
    def sinks(self):
        """Sink particle arrays."""
        return self.base.sinks

    def __repr__(self):
        """Dunder repr method."""
        return self.__str__()

    def __str__(self):
        """Dunder str method."""
        return f'<plonk.SubSnap "{self.file_path.name}">'

    def _get_array(self, name: str, sinks: bool = False) -> ndarray:
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

    >>> sinks = Sinks(snap=snap)

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
        self._getitem_from_str = base._getitem_from_str

        if indices is None:
            indices = np.arange(base.num_sinks)
        ind = _input_indices_array(inp=indices, max_slice=base.num_sinks)
        if ind is None:
            raise ValueError('Sinks has no particles')
        self._indices = ind

        # Attributes same as Snap
        self.file_path = self.base.file_path

    @property
    def indices(self):
        """Sink particle indices."""
        return self._indices

    def available_arrays(self, verbose: bool = False, aliases: bool = False):
        """Return a tuple of available sink arrays.

        Parameters
        ----------
        verbose
            Also display suffixed arrays, e.g. 'position_x',
            'position_y', etc. Default is False
        aliases
            If True, return array aliases. Default is False.

        Returns
        -------
        A tuple of names of available arrays.
        """
        return self.base._available_arrays(sinks=True, verbose=verbose, aliases=aliases)

    def loaded_arrays(self):
        """Return a tuple of loaded arrays.

        Returns
        -------
        A tuple of names of loaded sink arrays.
        """
        return tuple(sorted(self.base._sinks.keys()))

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
            return np.squeeze(self._getitem_from_str(inp, sinks=True)[self.indices])[()]
        ind = _input_indices_array(inp=inp, max_slice=len(self))
        if ind is not None:
            return Sinks(self.base, ind)
        raise ValueError('Cannot determine item to return')

    def __len__(self):
        """Length as number of particles."""
        return len(self.indices)

    def __repr__(self):
        """Dunder repr method."""
        return self.__str__()

    def __str__(self):
        """Dunder str method."""
        return f'<plonk.Sinks "{self.base.file_path.name}">'

    def _ipython_key_completions_(self):
        """Tab completion for IPython __getitem__ method."""
        return self.available_arrays(verbose=True)

    plot = visualize.plot


def _str_is_int(string):
    try:
        int(string)
        return True
    except ValueError:
        return False


def _input_indices_array(inp: Union[ndarray, slice, list, int, tuple], max_slice: int):
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
        i1, i2, step = inp.start, inp.stop, inp.step
        if i1 is None:
            i1 = 0
        if i2 is None:
            i2 = max_slice
        if step is not None:
            return np.arange(i1, i2, step)
        return np.arange(i1, i2)
    return None
