"""Snap, SubSnap, Sinks classes for snapshot files.

The Snap class contains all information related to a smoothed particle
hydrodynamics simulation snapshot file. The SubSnap class is for
accessing a subset of particles in a Snap.
"""

from __future__ import annotations

from typing import Any, Callable, Dict, List, Optional, Set, Tuple, Union, cast

import numpy as np
import pandas as pd
from numpy import ndarray
from pandas import DataFrame
from scipy.spatial.transform import Rotation

from .. import Quantity
from .. import units as plonk_units
from ..analysis import particles
from ..utils import norm


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

    To set a new array.

    >>> snap['r'] = np.sqrt(snap['x'] ** 2 + snap['y'] ** 2)

    Alternatively, define a function.

    >>> @plonk.Snap.add_array()
    ... def radius(snap) -> ndarray:
    ...     radius = np.hypot(snap['x'], snap['y'])
    ...     return radius

    Possibly with units.

    >>> @plonk.Snap.add_array(unit='length')
    ... def radius(snap) -> ndarray:
    ...     radius = np.hypot(snap['x'], snap['y'])
    ...     return radius

    Or, use an existing one.

    >>> snap['R'] = plonk.analysis.particles.radial_distance(snap)

    Set physical units. Arrays are now Pint quantities.

    >>> snap.physical_units()
    """

    _array_registry: Dict[str, Callable] = {}
    _sink_registry: Dict[str, Callable] = {}

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
        'v_R': 'radial_velocity_cylindrical',
        'v_r': 'radial_velocity_spherical',
        'v_phi': 'angular_velocity',
        'xyz': 'position',
    }

    _array_split_mapper = {
        'x': ('position', 0),
        'y': ('position', 1),
        'z': ('position', 2),
        'vx': ('velocity', 0),
        'vy': ('velocity', 1),
        'vz': ('velocity', 2),
        'px': ('momentum', 0),
        'py': ('momentum', 1),
        'pz': ('momentum', 2),
        'Lx': ('angular_momentum', 0),
        'Ly': ('angular_momentum', 1),
        'Lz': ('angular_momentum', 2),
        'jx': ('specific_angular_momentum', 0),
        'jy': ('specific_angular_momentum', 1),
        'jz': ('specific_angular_momentum', 2),
        'sx': ('spin', 0),
        'sy': ('spin', 1),
        'sz': ('spin', 2),
        'Bx': ('magnetic_field', 0),
        'By': ('magnetic_field', 1),
        'Bz': ('magnetic_field', 2),
    }

    _array_units = {
        'alpha_viscosity_numerical': 'dimensionless',
        'density': 'density',
        'differential_velocity': 'velocity',
        'dust_fraction': 'dimensionless',
        'dust_type': 'dimensionless',
        'gravitational_potential': 'energy',
        'internal_energy': 'energy',
        'magnetic_field': 'magnetic_field',
        'mass': 'mass',
        'position': 'length',
        'pressure': 'pressure',
        'smoothing_length': 'length',
        'sound_speed': 'velocity',
        'spin': 'angular_momentum',
        'stopping_time': 'time',
        'timestep': 'time',
        'type': 'dimensionless',
        'velocity': 'velocity',
        'velocity_divergence': 'frequency',
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
    def add_array(
        unit: str = None, rotatable: bool = None, dust: bool = False
    ) -> Callable:
        """Decorate function to add array to Snap.

        This function decorates a function that returns an array. The
        name of the function is the string with which to reference the
        array.

        Parameters
        ----------
        unit
            A string to represent the units of the array. E.g. 'length'
            for a 'radius' array. Default is None.
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
            Snap._array_registry[fn.__name__] = fn
            Snap._array_units[fn.__name__] = unit
            if rotatable is True:
                Snap._vector_arrays.add(fn.__name__)
            elif rotatable is False:
                Snap._vector_component_arrays.add(fn.__name__)
            if dust is True:
                Snap._dust_arrays.add(fn.__name__)
            return fn

        return _add_array

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
        self.properties = {}
        self.units = {}
        self._arrays = {}
        self._sinks = {}
        self._file_pointer = None
        self._num_particles = -1
        self._num_sinks = -1
        self._families = {key: None for key in Snap.particle_type.keys()}
        self.rotation = None
        self.translation = None
        self._physical_units = False
        self._extra_quantities = False

    def close_file(self):
        """Close access to underlying file."""
        self._file_pointer.close()

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

    def available_arrays(self, sinks: bool = False, aliases: bool = False):
        """Return a tuple of available arrays.

        Parameters
        ----------
        sinks
            If True, return available sink arrays. Default is
            False.
        aliases
            If True, return array aliases. Default is False.

        Returns
        -------
        A tuple of names of available arrays.

        Notes
        -----
        If an array is not available, it may be available after calling
        the extra_quantities method.
        """
        if sinks:
            loaded = self.loaded_arrays(sinks)
            registered = tuple(self._sink_registry.keys())
        else:
            loaded = self.loaded_arrays()
            registered = tuple(sorted(self._array_registry.keys()))

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

        return tuple(sorted(set(loaded + registered)))

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
    def num_sinks(self):
        """Return number of sinks."""
        if self._num_sinks == -1:
            try:
                self._num_sinks = len(self._sink_registry['mass'](self))
            except KeyError:
                self._num_sinks = 0
        return self._num_sinks

    def extra_quantities(self):
        """Make extra quantities available."""
        if self._extra_quantities:
            raise ValueError('Extra quantities already available')
        n_dust = len(self.properties.get('grain_size', []))
        dust = n_dust > 0
        dust_method = self.properties.get('dust_method')
        extra_quantities(dust=dust, dust_method=dust_method)
        self._extra_quantities = True
        return self

    def add_unit(self, name: str, unit: Any, unit_str: str):
        """Define a unit on an array.

        Parameters
        ----------
        name
            The name of the array.
        unit
            The Pint units Quantity.
        unit_str
            The unit string. See units attribute for units.

        Examples
        --------
        New array 'arr' with dimension 'length' and units 'cm'.

        >>> snap.add_unit('arr', plonk.units('cm'), 'length')
        """
        if name in self._array_units:
            raise ValueError('Array unit already defined on Snap')
        if unit_str not in self.units:
            self.units[unit_str] = unit
        self._array_units[name] = unit_str

    def unset(
        self, units: bool = False, rotation: bool = False, translation: bool = False,
    ):
        """Unset.

        Unset some transformations on the Snap data.

        Parameters
        ----------
        units
            Set to True to unset physical units. Default is False.
        rotation
            Set to True to unset rotation. Default is False.
        translation
            Set to True to unset translation. Default is False.
        """
        if any((units, rotation, translation)):
            for arr in self.loaded_arrays():
                del self._arrays[arr]
            for arr in self.loaded_arrays(sinks=True):
                del self._sinks[arr]
        else:
            raise ValueError('Select something to unset')

        if units:
            self._physical_units = None
        if rotation:
            self.rotation = None
        if translation:
            self.translation = None

        return self

    def physical_units(self) -> Snap:
        """Set physical units.

        Returns
        -------
        Snap
        """
        if self._physical_units:
            raise ValueError(
                'Physical units already set: snap.unset(units=True) to unset.'
            )
        for arr in self.loaded_arrays():
            self._arrays[arr] = self._arrays[arr] * self.get_array_unit(arr)
        for arr in self.loaded_arrays(sinks=True):
            self._sinks[arr] = self._sinks[arr] * self.get_array_unit(arr)
        self._physical_units = True

        return self

    def rotate(self, rotation: Rotation) -> Snap:
        """Rotate snapshot.

        Parameters
        ----------
        rotation
            The rotation as a scipy.spatial.transform.Rotation object.

        Returns
        -------
        Snap
            The rotated Snap. Note that the rotation operation is
            in-place.
        """
        for arr in self._vector_arrays:
            if arr in self.loaded_arrays():
                self._arrays[arr] = rotation.apply(self._arrays[arr])
            if arr in self.loaded_arrays(sinks=True):
                self._sinks[arr] = rotation.apply(self._sinks[arr])
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

    def translate(self, translation: ndarray) -> Snap:
        """Translate snapshot.

        I.e. shift the snapshot origin.

        Parameters
        ----------
        translation
            The translation as a (3,) ndarray like (x, y, z).

        Returns
        -------
        Snap
            The translated Snap. Note that the translation operation is
            in-place.
        """
        translation = np.array(translation)
        if translation.shape != (3,):
            raise ValueError('translation must be like (x, y, z)')
        if 'position' in self.loaded_arrays():
            self._arrays['position'] += translation
        if 'position' in self.loaded_arrays(sinks=True):
            self._sinks['position'] += translation

        if self.translation is None:
            self.translation = translation
        else:
            self.translation += translation

        return self

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
        if self._physical_units:
            G = (G * M).to_base_units()
        else:
            G = (G * self.units['mass'] * M).to_base_units()
        self.properties['gravitational_parameter'] = G

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
        self.properties['molecular_weight'] = molecular_weight

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

    def _get_family_indices(self, name: str):
        """Get a family by name."""
        if name in self._families:
            if self._families[name] is None:
                self._families[name] = np.flatnonzero(
                    self['type'] == Snap.particle_type[name]
                )
            return self._families[name]
        else:
            raise ValueError('Family not available')

    def get_array_unit(self, arr: str) -> Any:
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
        if arr in self._array_split_mapper:
            arr = self._array_split_mapper[arr][0]
        elif arr in self._array_name_mapper:
            arr = self._array_name_mapper[arr]
        try:
            unit = self.units[self._array_units[arr]]
        except KeyError:
            unit = 1.0
        return unit

    def _get_array_from_registry(self, name: str, sinks: bool = False):
        if sinks:
            array = Snap._sink_registry[name](self)
            array_dict = self._sinks
        else:
            array = Snap._array_registry[name](self)
            array_dict = self._arrays
        if self.rotation is not None and name in self._vector_arrays:
            array = self.rotation.apply(array)
        if self.translation is not None and name == 'position':
            array += self.translation
        if self._physical_units and not isinstance(array, Quantity):
            unit = self.get_array_unit(name)
            array_dict[name] = unit * array
        else:
            array_dict[name] = array

    def _get_array(self, name: str, sinks: bool = False) -> ndarray:
        """Get an array by name."""
        index = None

        if name in self.available_arrays(sinks):
            name = name
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
        elif name in Snap._array_registry or name in Snap._sink_registry:
            self._get_array_from_registry(name, sinks)
            if index is None:
                return array_dict[name]
            return array_dict[name][:, index]
        else:
            raise ValueError('Array not available')

    def _getitem_from_str(self, inp: str, sinks: bool = False) -> ndarray:
        """Return item from string."""
        inp_root = '_'.join(inp.split('_')[:-1])
        inp_suffix = inp.split('_')[-1]

        if inp in self._families:
            return SubSnap(self, self._get_family_indices(inp))
        elif inp in self.available_arrays(sinks):
            return self._get_array(inp, sinks)
        elif inp in self._array_name_mapper.keys():
            return self._get_array(inp, sinks)
        elif inp in self._array_split_mapper.keys():
            return self._get_array(inp, sinks)
        elif inp_root in self._vector_arrays:
            if inp_suffix == 'x':
                return self._get_array(inp_root, sinks)[:, 0]
            elif inp_suffix == 'y':
                return self._get_array(inp_root, sinks)[:, 1]
            elif inp_suffix == 'z':
                return self._get_array(inp_root, sinks)[:, 2]
            elif inp_suffix == 'magnitude':
                return norm(self._get_array(inp_root, sinks), axis=1)
        elif inp_root in self._dust_arrays:
            if _str_is_int(inp_suffix):
                return self._get_array(inp_root)[:, int(inp_suffix) - 1]
            elif inp_suffix == 'total':
                return self._get_array(inp_root).sum(axis=1)

        if self._extra_quantities:
            raise ValueError('Cannot determine item to return.')
        else:
            raise ValueError(
                'Cannot determine item to return. Extra quantities are available via\n'
                'snap.extra_quantities().'
            )

    def _getitem(
        self, inp: Union[str, ndarray, int, slice], sinks: bool = False,
    ) -> Union[ndarray, SubSnap]:
        """Return an array, or family, or subset."""
        if isinstance(inp, str):
            return self._getitem_from_str(inp, sinks)
        elif sinks:
            raise ValueError('Cannot return sinks as SubSnap')
        elif isinstance(inp, ndarray):
            if np.issubdtype(np.bool, inp.dtype):
                return SubSnap(self, np.flatnonzero(inp))
            elif np.issubdtype(np.int, inp.dtype):
                return SubSnap(self, inp)
        elif isinstance(inp, int):
            return SubSnap(self, np.array([inp]))
        elif isinstance(inp, slice):
            i1, i2, step = inp.start, inp.stop, inp.step
            if i1 is None:
                i1 = 0
            if i2 is None:
                i2 = len(self)
            if step is not None:
                return SubSnap(self, np.arange(i1, i2, step))
            return SubSnap(self, np.arange(i1, i2))
        if self._extra_quantities:
            raise ValueError('Cannot determine item to return.')
        else:
            raise ValueError(
                'Cannot determine item to return. Extra quantities are available via\n'
                'snap.extra_quantities().'
            )

    def __getitem__(
        self, inp: Union[str, ndarray, int, slice]
    ) -> Union[ndarray, SubSnap]:
        """Return an array, or family, or subset."""
        return self._getitem(inp, sinks=False)

    def __setitem__(self, name: str, item: ndarray):
        """Set an array."""
        if not isinstance(item, ndarray):
            raise ValueError('"item" must be ndarray')
        if item.shape[0] != len(self):
            raise ValueError('Length of array does not match particle number')
        if name in self.loaded_arrays():
            raise ValueError(
                'Attempting to overwrite existing array. To do so, first delete the '
                'array\nwith del snap["array"], then try again.'
            )
        elif (
            name in self.available_arrays()
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

    def __len__(self):
        """Length as number of particles."""
        return self.num_particles

    def __repr__(self):
        """Dunder repr method."""
        return self.__str__()

    def __str__(self):
        """Dunder str method."""
        return f'<plonk.Snap "{self.file_path.name}">'


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

        self.base = base
        self._indices = indices
        self._num_particles = len(indices)

        self.data_source = self.base.data_source
        self.file_path = self.base.file_path
        self.properties = self.base.properties
        self.units = self.base.units
        self._file_pointer = self.base._file_pointer
        self._arrays = self.base._arrays
        self._sinks = self.base._sinks
        self.rotation = self.base.rotation
        self.translation = self.base.translation
        self._physical_units = self.base._physical_units
        self._extra_quantities = self.base._extra_quantities

    def __repr__(self):
        """Dunder repr method."""
        return self.__str__()

    def __str__(self):
        """Dunder str method."""
        return f'<plonk.SubSnap>'

    def _get_array(self, name: str, sinks: bool = False) -> ndarray:
        return self.base._get_array(name, sinks)[self._indices]


SnapLike = Union[Snap, SubSnap]


def get_array_from_input(
    snap: SnapLike, inp: Union[str, ndarray], default: str = None
) -> ndarray:
    """Get array on Snap.

    Parameters
    ----------
    snap
        The Snap or SubSnap.
    inp
        The input as a string or ndarray. If a string return
        snap[inp], otherwise return inp as a ndarray.
    default
        The default array as a string resolved as snap[default].

    Returns
    -------
    ndarray
        The array on the particles.
    """
    if isinstance(inp, ndarray):
        return inp
    elif isinstance(inp, str):
        return get_array_in_code_units(snap, inp)
    elif default is not None:
        return get_array_in_code_units(snap, default)
    raise ValueError('Cannot determine array to return')


def get_array_in_code_units(snap: SnapLike, name: str) -> ndarray:
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
    arr = snap[name]
    if isinstance(arr, Quantity):
        return (arr / snap.get_array_unit(name)).magnitude
    return arr


def extra_quantities(
    dust: bool = False, dust_method: Optional[str] = None,
):
    """Make extra quantities available.

    Parameters
    ----------
    dust
        Whether to add dust quantities.
    dust_method
        The dust method.
    """

    @Snap.add_array(unit='momentum', rotatable=True)
    def momentum(snap) -> ndarray:
        """Momentum."""
        return particles.momentum(snap=snap)

    @Snap.add_array(unit='angular_momentum', rotatable=True)
    def angular_momentum(snap) -> ndarray:
        """Angular momentum."""
        origin = snap.translation if snap.translation is not None else (0.0, 0.0, 0.0)
        return particles.angular_momentum(snap=snap, origin=origin)

    @Snap.add_array(unit='specific_angular_momentum', rotatable=True)
    def specific_angular_momentum(snap) -> ndarray:
        """Specific angular momentum."""
        origin = snap.translation if snap.translation is not None else (0.0, 0.0, 0.0)
        return particles.specific_angular_momentum(snap=snap, origin=origin)

    @Snap.add_array(unit='energy')
    def kinetic_energy(snap) -> ndarray:
        """Kinetic energy."""
        return particles.kinetic_energy(snap=snap)

    @Snap.add_array(unit='specific_energy')
    def specific_kinetic_energy(snap) -> ndarray:
        """Specific kinetic energy."""
        return particles.specific_kinetic_energy(snap=snap)

    @Snap.add_array(unit='frequency')
    def keplerian_frequency(snap) -> ndarray:
        """Keplerian orbital frequency."""
        gravitational_parameter = snap.properties.get('gravitational_parameter')
        if gravitational_parameter is None:
            raise ValueError(
                'To get Keplerian frequency, first set the gravitational parameter\n'
                'via snap.set_gravitational_parameter.'
            )
        origin = snap.translation if snap.translation is not None else (0.0, 0.0, 0.0)
        return particles.keplerian_frequency(
            snap=snap, gravitational_parameter=gravitational_parameter, origin=origin
        )

    @Snap.add_array(unit='length')
    def semi_major_axis(snap) -> ndarray:
        """Semi-major axis."""
        gravitational_parameter = snap.properties.get('gravitational_parameter')
        if gravitational_parameter is None:
            raise ValueError(
                'To get semi-major axis, first set the gravitational parameter\n'
                'via snap.set_gravitational_parameter.'
            )
        origin = snap.translation if snap.translation is not None else (0.0, 0.0, 0.0)
        return particles.semi_major_axis(
            snap=snap, gravitational_parameter=gravitational_parameter, origin=origin
        )

    @Snap.add_array(unit='dimensionless')
    def eccentricity(snap) -> ndarray:
        """Eccentricity."""
        gravitational_parameter = snap.properties.get('gravitational_parameter')
        if gravitational_parameter is None:
            raise ValueError(
                'To get eccentricity, first set the gravitational parameter\n'
                'via snap.set_gravitational_parameter.'
            )
        origin = snap.translation if snap.translation is not None else (0.0, 0.0, 0.0)
        return particles.eccentricity(
            snap=snap, gravitational_parameter=gravitational_parameter, origin=origin
        )

    @Snap.add_array(unit='radian', rotatable=False)
    def inclination(snap) -> ndarray:
        """Inclination."""
        return particles.inclination(snap=snap)

    @Snap.add_array(unit='length', rotatable=False)
    def radius_cylindrical(snap) -> ndarray:
        """Cylindrical radius."""
        return particles.radial_distance(snap=snap, coordinates='cylindrical')

    @Snap.add_array(unit='length', rotatable=False)
    def radius_spherical(snap) -> ndarray:
        """Spherical radius."""
        return particles.radial_distance(snap=snap, coordinates='spherical')

    @Snap.add_array(unit='radian', rotatable=False)
    def azimuthal_angle(snap) -> ndarray:
        """Azimuthal angle."""
        return particles.azimuthal_angle(snap=snap)

    @Snap.add_array(unit='radian', rotatable=False)
    def polar_angle(snap) -> ndarray:
        """Polar angle."""
        return particles.polar_angle(snap=snap)

    @Snap.add_array(unit='velocity', rotatable=False)
    def radial_velocity_cylindrical(snap) -> ndarray:
        """Cylindrical radial velocity."""
        return particles.radial_velocity(snap=snap, coordinates='cylindrical')

    @Snap.add_array(unit='velocity', rotatable=False)
    def radial_velocity_spherical(snap) -> ndarray:
        """Spherical radial velocity."""
        return particles.radial_velocity(snap=snap, coordinates='spherical')

    @Snap.add_array(unit='frequency', rotatable=False)
    def angular_velocity(snap) -> ndarray:
        """Angular velocity."""
        return particles.angular_velocity(snap=snap)

    @Snap.add_array(unit='temperature')
    def temperature(snap) -> ndarray:
        """Temperature."""
        molecular_weight = snap.properties.get('molecular_weight')
        if molecular_weight is None:
            raise ValueError(
                'To get temperature, first set the molecular weight parameter\n'
                'via snap.set_molecular_weight method.'
            )
        return particles.temperature(snap=snap, molecular_weight=molecular_weight)

    if dust:

        @Snap.add_array(unit='dimensionless')
        def gas_fraction(snap) -> ndarray:
            """Gas fraction."""
            return particles.gas_fraction(snap=snap)

        @Snap.add_array(unit='mass')
        def gas_mass(snap) -> ndarray:
            """Gas mass."""
            return particles.gas_mass(snap=snap)

        @Snap.add_array(unit='density')
        def gas_density(snap) -> ndarray:
            """Gas density."""
            return particles.gas_density(snap=snap)

        if dust_method == 'dust as separate sets of particles':

            @Snap.add_array(unit='dimensionless', dust=True)
            def dust_fraction(snap) -> ndarray:
                """Dust fraction."""
                return particles.dust_fraction(snap=snap)

        @Snap.add_array(unit='mass', dust=True)
        def dust_mass(snap) -> ndarray:
            """Dust mass."""
            return particles.dust_mass(snap=snap)

        @Snap.add_array(unit='density', dust=True)
        def dust_density(snap) -> ndarray:
            """Dust density."""
            return particles.dust_density(snap=snap)

        @Snap.add_array(unit='dimensionless', dust=True)
        def stokes_number(snap) -> ndarray:
            """Stokes number."""
            gravitational_parameter = snap.properties.get('gravitational_parameter')
            if gravitational_parameter is None:
                raise ValueError(
                    'To get eccentricity, first set the gravitational parameter\n'
                    'via snap.set_gravitational_parameter.'
                )
            origin = (
                snap.translation if snap.translation is not None else (0.0, 0.0, 0.0)
            )
            return particles.stokes_number(
                snap=snap,
                gravitational_parameter=gravitational_parameter,
                origin=origin,
            )


def _str_is_int(string):
    try:
        int(string)
        return True
    except ValueError:
        return False
