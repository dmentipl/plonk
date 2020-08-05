"""Snap, SubSnap, Sinks classes for snapshot files.

The Snap class contains all information related to a smoothed particle
hydrodynamics simulation snapshot file. The SubSnap class is for
accessing a subset of particles in a Snap.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Callable, Dict, List, Set, Tuple, Union, cast

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


class _SinkUtility:
    def __init__(self, fn):
        self.fn = fn

    def __getitem__(self, inp):
        return self.fn(inp, sinks=True)

    def __repr__(self):
        """Dunder repr method."""
        return self.__str__()

    def __str__(self):
        """Dunder str method."""
        return f'<plonk.snap sinks>'


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
        self.units = {}
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

    def loaded_arrays(self, sinks: bool = False):
        """Return a tuple of loaded arrays.

        Parameters
        ----------
        sinks
            If True, return loaded sink arrays.

        Returns
        -------
        A tuple of names of loaded arrays.
        """
        if sinks:
            return tuple(sorted(self._sinks.keys()))
        return tuple(sorted(self._arrays.keys()))

    def _available_arrays(
        self, sinks: bool = False, all: bool = False, aliases: bool = False
    ):
        """Return a tuple of available arrays.

        Parameters
        ----------
        sinks
            If True, return available sink arrays. Default is
            False.
        all
            TODO. Default is False
        aliases
            If True, return array aliases. Default is False.

        Returns
        -------
        A tuple of names of available arrays.
        """
        if sinks:
            loaded = self.loaded_arrays(sinks)
            registered = list(self._sink_registry.keys())
        else:
            loaded = self.loaded_arrays()
            registered = list(sorted(self._array_registry.keys()))
            if all:
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

    def available_arrays(self, all: bool = False, aliases: bool = False):
        """Return a tuple of available particle arrays.

        Parameters
        ----------
        all
            TODO. Default is False
        aliases
            If True, return array aliases. Default is False.

        Returns
        -------
        A tuple of names of available arrays.
        """
        return self._available_arrays(sinks=False, all=all, aliases=aliases)

    def available_sink_arrays(self, aliases: bool = False):
        """Return a tuple of available sink arrays.

        Parameters
        ----------
        aliases
            If True, return array aliases. Default is False.

        Returns
        -------
        A tuple of names of available arrays.
        """
        return self._available_arrays(sinks=True, aliases=aliases)

    @property
    def properties(self):
        """Snap properties."""
        return {key: self._properties[key] for key in sorted(self._properties.keys())}

    @property
    def sinks(self):
        """Sink particle arrays."""
        return _SinkUtility(self._getitem)

    @property
    def num_particles(self):
        """Return number of particles."""
        if self._num_particles == -1:
            self._num_particles = len(self._array_registry['type'](self))
        return self._num_particles

    @property
    def num_particles_of_type(self):
        """Return number of particles per type."""
        if self._num_particles_of_type == -1:
            int_to_name = {idx: name for name, idx in self.particle_type.items()}
            d = {}
            for idx, num in enumerate(np.bincount(self['type'])):
                if num > 0:
                    if idx == self.particle_type['dust']:
                        # Dust particle sub-type skips zero: 1, 2, 3 ...
                        d['dust'] = list(
                            np.bincount(self[self['type'] == idx]['sub_type'])[1:]
                        )
                    elif idx == self.particle_type['boundary']:
                        # Boundary particle sub-type: 0 (gas), 1, 2, 3... (dust)
                        d['boundary'] = list(
                            np.bincount(self[self['type'] == idx]['sub_type'])
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

    def unset(self, rotation: bool = False, translation: bool = False):
        """Unset.

        Unset some transformations on the Snap data.

        Parameters
        ----------
        rotation
            Set to True to unset rotation. Default is False.
        translation
            Set to True to unset translation. Default is False.
        """
        if any((rotation, translation)):
            for arr in self.loaded_arrays():
                del self._arrays[arr]
            for arr in self.loaded_arrays(sinks=True):
                del self._sinks[arr]
        else:
            raise ValueError('Select something to unset')

        if rotation:
            self.rotation = None
        if translation:
            self.translation = None

        return self

    def rotate(self, rotation: Union[ndarray, Rotation]) -> Snap:
        """Rotate snapshot.

        Parameters
        ----------
        rotation
            The rotation as a scipy.spatial.transform.Rotation object
            or ndarray that can be converted to a Rotation object via
            Rotation.from_rotvec.

        Returns
        -------
        Snap
            The rotated Snap. Note that the rotation operation is
            in-place.

        Examples
        --------
        Rotate a Snap by π/3 around [1, 1, 0].

        >>> rot = np.array([1, 1, 0])
        >>> rot = rot * np.pi / 3 * np.linalg.norm(rot)
        >>> snap.rotate(rot)
        """
        logger.debug(f'Rotating snapshot: {self.file_path.name}')
        if isinstance(rotation, (list, tuple, ndarray)):
            rotation = Rotation.from_rotvec(rotation)
        for arr in self._vector_arrays:
            if arr in self.loaded_arrays():
                array_m, array_u = self._arrays[arr].magnitude, self._arrays[arr].units
                self._arrays[arr] = rotation.apply(array_m) * array_u
            if arr in self.loaded_arrays(sinks=True):
                array_m, array_u = self._sinks[arr].magnitude, self._sinks[arr].units
                self._sinks[arr] = rotation.apply(array_m) * array_u
        for arr in self._vector_component_arrays:
            if arr in self.loaded_arrays():
                del self._arrays[arr]
            if arr in self.loaded_arrays(sinks=True):
                del self._sinks[arr]

        if self.rotation is None:
            self.rotation = rotation
        else:
            rot = rotation * self.rotation
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
            else:
                translation *= plonk_units(unit)
        if 'position' in self.loaded_arrays():
            self._arrays['position'] += translation
        if 'position' in self.loaded_arrays(sinks=True):
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
            M = self.sinks['mass'][sink_idx]
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
        """Particle neighbour kd-trees.

        Trees are represented by scipy cKDTree objects. There is one
        tree per particle type stored in a dictionary. Particles with
        subtypes, e.g. dust, have one tree per subtype.
        """
        if self._tree is None:
            self._tree = dict()
            for key in self.particle_type:
                ind = self.particle_indices(key)
                if len(ind) == 0:
                    self._tree[key] = None
                elif isinstance(ind, list):
                    self._tree[key] = list()
                    for _ind in ind:
                        self._tree[key].append(cKDTree(self['position'][_ind]))
                else:
                    self._tree[key] = cKDTree(self['position'][ind])
        return self._tree

    def get_neighbours(self, idx):
        """Get neighbours of a particle.

        Parameters
        ----------
        idx
            The particle index.

        Returns
        -------
        list
            The list of neighbours relative to the particle
            type/sub-type.
        """
        kernel = self._properties.get('kernel')
        if kernel is None:
            raise ValueError(
                'To calculate particle neighbours, first set the kernel\n'
                'via snap.set_kernel.'
            )
        r_kern = kernel_radius[kernel]
        int_to_str_type = {val: key for key, val in self.particle_type.items()}
        particle_type = int_to_str_type[self['type'][idx]]
        tree = self.tree[particle_type]
        if particle_type == 'dust':
            sub_type = self['sub_type'][idx]
            tree = tree[sub_type]
        if particle_type == 'boundary':
            sub_type = self['sub_type'][idx]
            tree = tree[sub_type]
        neighbours = tree.query_ball_point(
            self['position'][idx], r_kern * self['smoothing_length'][idx], n_jobs=-1,
        )
        return neighbours

    def get_many_neighbours(self, indices):
        """Get neighbours of more than one particle.

        Assumes all particles are of same type.

        Parameters
        ----------
        indices
            The particle indices.

        Returns
        -------
        ndarray of list
            An array of neighbours lists relative to the particle
            type/sub-type.
        """
        kernel = self._properties.get('kernel')
        if kernel is None:
            raise ValueError(
                'To calculate particle neighbours, first set the kernel\n'
                'via snap.set_kernel.'
            )
        r_kern = kernel_radius[kernel]
        int_to_str_type = {val: key for key, val in self.particle_type.items()}
        particle_type = int_to_str_type[self['type'][indices[0]]]
        tree = self.tree[particle_type]
        if particle_type == 'dust':
            sub_type = self['sub_type'][indices[0]]
            tree = tree[sub_type]
        if particle_type == 'boundary':
            sub_type = self['sub_type'][indices[0]]
            tree = tree[sub_type]
        neighbours = tree.query_ball_point(
            self['position'][indices],
            r_kern * self['smoothing_length'][indices],
            n_jobs=-1,
        )
        return neighbours

    def get_all_neighbours(self):
        """Get all particle neighbours.

        Returns
        -------
        ndarray of lists
            An 1d array where the index is the particle, and the item
            is a list of neighbours.

        Notes
        -----
        This has a large memory requirement. It is likely to crash if
        the Snap has more than approximately 1 million particles
        depending on the available memory.
        """
        kernel = self._properties.get('kernel')
        if kernel is None:
            raise ValueError(
                'To calculate particle neighbours, first set the kernel\n'
                'via snap.set_kernel.'
            )
        r_kern = kernel_radius[kernel]
        logger.info('Finding neighbours... may take some time...', end='', flush=True)

        neighbours = np.zeros(len(self), dtype=object)
        for key, tree in self.tree.items():
            ind = self.particle_indices(key)
            if len(ind) == 0:
                continue
            if isinstance(ind, list):
                for idx, _ind in enumerate(ind):
                    h = self['smoothing_length'][_ind]
                    pos = self['position'][_ind]
                    neighbours[_ind] = self.tree[key][idx].query_ball_point(
                        pos, r_kern * h, n_jobs=-1
                    )
            else:
                h = self['smoothing_length'][ind]
                pos = self['position'][ind]
                neighbours[ind] = self.tree[key].query_ball_point(
                    pos, r_kern * h, n_jobs=-1
                )

        _neighbours = np.zeros(len(self), dtype=object)
        ind = {key: self.particle_indices(key) for key in self.particle_type}
        for idx, neigh in enumerate(neighbours):
            if self['type'][idx] == self.particle_type['gas']:
                _neighbours[idx] = ind['gas'][neigh]
            elif self['type'][idx] == self.particle_type['dust']:
                _neighbours[idx] = ind['dust'][self['sub_type'][idx] - 1][neigh]
            elif self['type'][idx] == self.particle_type['boundary']:
                _neighbours[idx] = ind['boundary'][self['sub_type'][idx] - 1][neigh]
            elif self['type'][idx] == self.particle_type['star']:
                _neighbours[idx] = ind['star'][neigh]
            elif self['type'][idx] == self.particle_type['darkmatter']:
                _neighbours[idx] = ind['darkmatter'][neigh]
            elif self['type'][idx] == self.particle_type['bulge']:
                _neighbours[idx] = ind['bulge'][neigh]

        logger.info(' Done!', flush=True)

        return _neighbours

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
            arr: ndarray = self[array]
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
        self, columns: Union[Tuple[str, ...], List[str]] = None
    ) -> DataFrame:
        """Convert Snap to DataFrame.

        Parameters
        ----------
        columns : optional
            A list of columns to add to the data frame. Default is
            None.

        Returns
        -------
        DataFrame
        """
        d = dict()
        if columns is None:
            columns = self.loaded_arrays()
        cols = list(columns)
        for col in cols:
            arr = self[col]
            arr = cast(ndarray, arr)
            if arr.ndim == 2:
                for idx in range(arr.shape[1]):
                    d[f'{col}.{idx+1}'] = arr[:, idx]
            else:
                d[col] = arr
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
                unit *= self.units[d] ** dim[f'[{d}]']
        return unit

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
        if isinstance(inp, ndarray):
            if np.issubdtype(np.bool, inp.dtype):
                return SubSnap(self, np.flatnonzero(inp))
            if np.issubdtype(np.int, inp.dtype):
                return SubSnap(self, inp)
        if isinstance(inp, int):
            return SubSnap(self, np.array([inp]))
        if isinstance(inp, slice):
            i1, i2, step = inp.start, inp.stop, inp.step
            if i1 is None:
                i1 = 0
            if i2 is None:
                i2 = len(self)
            if step is not None:
                return SubSnap(self, np.arange(i1, i2, step))
            return SubSnap(self, np.arange(i1, i2))
        raise ValueError('Cannot determine item to return.')

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
        return self.available_arrays()

    def __len__(self):
        """Length as number of particles."""
        return self.num_particles

    def __repr__(self):
        """Dunder repr method."""
        return self.__str__()

    def __str__(self):
        """Dunder str method."""
        return f'<plonk.Snap "{self.file_path.name}">'

    plot = visualize.plot
    particle_plot = visualize.particle_plot


class SubSnap(Snap):
    """A Snap subset of particles.

    The sub-snap is generated via an index array.

    Parameters
    ----------
    base
        The base snapshot.
    indices
        A (N,) array of particle indices to include in the sub-snap.
    """

    def __init__(self, base: Snap, indices: ndarray):
        super().__init__()

        # Attributes different to Snap
        self.base = base
        self._indices = indices
        self._num_particles = len(indices)
        self._num_particles_of_type = -1
        self._num_dust_species = -1
        self._tree = None

        # Attributes same as Snap
        self.data_source = self.base.data_source
        self.file_path = self.base.file_path
        self.units = self.base.units
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
        return self.base.__str__().replace('Snap', 'SubSnap')

    def _get_array(self, name: str, sinks: bool = False) -> ndarray:
        return self.base._get_array(name, sinks)[self.indices]


SnapLike = Union[Snap, SubSnap]


def _str_is_int(string):
    try:
        int(string)
        return True
    except ValueError:
        return False
