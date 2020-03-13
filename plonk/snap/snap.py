"""Snap, SubSnap, Sinks classes for snapshot files.

The Snap class contains all information related to a smoothed particle
hydrodynamics simulation snapshot file. The SubSnap class is for
accessing a subset of particles in a Snap.
"""

from __future__ import annotations

from typing import Any, Callable, Dict, List, Optional, Tuple, Union, cast

import numpy as np
import pandas as pd
from numpy import ndarray
from pandas import DataFrame
from scipy.spatial.transform import Rotation

from .. import Quantity


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
        'xyz': 'position',
        'pos': 'position',
        'vxyz': 'velocity',
        'vel': 'velocity',
        'h': 'smooth',
        'm': 'mass',
        'rho': 'density',
        'Bxyz': 'magfield',
        'spinxyz': 'spin',
    }

    _array_split_mapper = {
        'x': ('position', 0),
        'y': ('position', 1),
        'z': ('position', 2),
        'vx': ('velocity', 0),
        'vy': ('velocity', 1),
        'vz': ('velocity', 2),
        'velx': ('velocity', 0),
        'vely': ('velocity', 1),
        'velz': ('velocity', 2),
        'Bx': ('magfield', 0),
        'By': ('magfield', 1),
        'Bz': ('magfield', 2),
        'sx': ('spin', 0),
        'sy': ('spin', 1),
        'sz': ('spin', 2),
    }

    _array_units = {
        'alpha': 'dimensionless',
        'density': 'density',
        'divv': 'frequency',
        'dt': 'time',
        'dustfrac': 'dimensionless',
        'dust_type': 'dimensionless',
        'magfield': 'magnetic_field',
        'mass': 'mass',
        'position': 'length',
        'poten': 'energy',
        'smooth': 'length',
        'spin': 'angular_momentum',
        'tstop': 'time',
        'type': 'dimensionless',
        'velocity': 'velocity',
    }

    _particle_type = {
        'gas': 1,
        'dust': 2,
        'boundary': 3,
        'star': 4,
        'darkmatter': 5,
        'bulge': 6,
    }

    @staticmethod
    def add_array(unit: str = None) -> Callable:
        """Decorate function to add array to Snap.

        This function decorates a function that returns an array. The
        name of the function is the string with which to reference the
        array.

        Parameters
        ----------
        unit
            A string to represent the units of the array. E.g. 'length'
            for a 'radius' array.

        Returns
        -------
        Callable
            The decorator which returns the array.
        """

        def _add_array(fn):
            Snap._array_units[fn.__name__] = unit
            Snap._array_registry[fn.__name__] = fn
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

        self.properties = {}
        self.units = {}
        self._arrays = {}
        self._sinks = {}
        self._file_pointer = None
        self._num_particles = -1
        self._num_sinks = -1
        self._families = {key: None for key in Snap._particle_type.keys()}
        self.rotation = None
        self.translation = None
        self._physical_units = False

    def close_file(self):
        """Close access to underlying file."""
        self._file_pointer.close()

    def loaded_arrays(self, sinks: bool = False):
        """Return a list of loaded arrays."""
        if sinks:
            return tuple(sorted(self._sinks.keys()))
        return tuple(sorted(self._arrays.keys()))

    def available_arrays(self, sinks: bool = False):
        """Return a list of available arrays."""
        if sinks:
            loaded = self.loaded_arrays(sinks)
            registered = tuple(sorted(self._sink_registry.keys()))
        else:
            loaded = self.loaded_arrays()
            registered = tuple(sorted(self._array_registry.keys()))
        return tuple(sorted(set(loaded + registered)))

    @property
    def sinks(self):
        """Sink particle arrays."""
        return _SinkUtility(self._get_array)

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

    def physical_units(self, unset: bool = False) -> Snap:
        """Set/unset physical units.

        Uses Pint.

        Returns
        -------
        Snap
        """
        if self._physical_units:
            raise ValueError(
                'Physical units already set: snap.physical_units(unset=True) to unset.'
            )
        if unset:
            for arr in self.loaded_arrays():
                del self._arrays[arr]
            for arr in self.loaded_arrays(sinks=True):
                del self._sinks[arr]
            self._physical_units = False
        else:
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
        for arr in self._rotation_required():
            if arr in self.loaded_arrays():
                self._arrays[arr] = rotation.apply(self._arrays[arr])
            if arr in self.loaded_arrays(sinks=True):
                self._sinks[arr] = rotation.apply(self._sinks[arr])

        self.rotation = rotation

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
        if len(translation) != 3:
            raise ValueError('translation must be like (x, y, z)')
        if 'position' in self.loaded_arrays():
            self._arrays['position'] += translation
        if 'position' in self.loaded_arrays(sinks=True):
            self._sinks['position'] += translation

        self.translation = translation

        return self

    def to_dataframe(self, columns: Union[Tuple[str, ...], List[str]]) -> DataFrame:
        """Convert Snap to DataFrame.

        Parameters
        ----------
        columns
            A list of columns to add to the data frame.

        Returns
        -------
        DataFrame
        """
        d = dict()
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

    def _rotation_required(self):
        return set([val[0] for val in self._array_split_mapper.values()])

    def _get_family_indices(self, name: str):
        """Get a family by name."""
        if name in self._families:
            if self._families[name] is None:
                self._families[name] = np.flatnonzero(
                    self['type'] == Snap._particle_type[name]
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
        if self.rotation is not None and name in self._rotation_required():
            array = self.rotation.apply(array)
        if self.translation is not None and name == 'position':
            array += self.translation
        if self._physical_units and not isinstance(array, Quantity):
            unit = self.get_array_unit(name)
            array_dict[name] = unit * array
        else:
            array_dict[name] = array

    def _get_array(
        self, name: str, index: Optional[int] = None, sinks: bool = False
    ) -> ndarray:
        """Get an array by name."""
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

    def __getitem__(
        self, inp: Union[str, ndarray, int, slice]
    ) -> Union[ndarray, SubSnap]:
        """Return an array, or family, or subset."""
        if isinstance(inp, str):
            if inp in self._families:
                return SubSnap(self, self._get_family_indices(inp))
            elif inp in self.available_arrays():
                return self._get_array(inp)
            elif inp in self._array_name_mapper.keys():
                return self._get_array(inp)
            elif inp in self._array_split_mapper.keys():
                return self._get_array(inp)
        elif isinstance(inp, ndarray):
            if np.issubdtype(np.bool, inp.dtype):
                return SubSnap(self, np.flatnonzero(inp))
            elif np.issubdtype(np.int, inp.dtype):
                return SubSnap(self, inp)
        elif isinstance(inp, int):
            raise NotImplementedError
        elif isinstance(inp, slice):
            i1, i2, step = inp.start, inp.stop, inp.step
            if step is not None:
                return SubSnap(self, np.arange(i1, i2, step))
            return SubSnap(self, np.arange(i1, i2))
        raise ValueError('Cannot determine item to return')

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
        return f'<plonk.Snap>'


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
        self.properties = self.base.properties
        self._file_pointer = self.base._file_pointer
        self._indices = indices
        self._num_particles = len(indices)
        self._arrays = self.base._arrays

    def __repr__(self):
        """Dunder repr method."""
        return self.__str__()

    def __str__(self):
        """Dunder str method."""
        return f'<plonk.SubSnap>'

    def _get_array(
        self, name: str, index: Optional[int] = None, sinks: bool = False
    ) -> ndarray:
        return self.base._get_array(name, index, sinks)[self._indices]


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
