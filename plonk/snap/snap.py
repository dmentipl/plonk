"""Snap, SubSnap, Sinks classes for snapshot files.

The Snap class contains all information related to a smoothed particle
hydrodynamics simulation snapshot file. The SubSnap class is for
accessing a subset of particles in a Snap.
"""

from __future__ import annotations

from typing import Callable, Dict, List, Optional, Tuple, Union, cast

import numpy as np
import pandas as pd
from numpy import ndarray
from pandas import DataFrame
from scipy.spatial.transform import Rotation


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

    >>> @plonk.Snap.add_array
    ... def radius(snap) -> ndarray:
    ...     radius = np.hypot(snap['x'], snap['y'])
    ...     return radius
    """

    _array_registry: Dict[str, Callable] = {}

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

    _particle_id = {
        'gas': 1,
        'dust': 2,
        'boundary': 3,
        'star': 4,
        'darkmatter': 5,
        'bulge': 6,
    }

    @staticmethod
    def add_array(fn: Callable) -> Callable:
        """Decorate function to add array to Snap.

        Parameters
        ----------
        fn
            A function that returns the array. The name of the function
            is the string with which to reference the array.

        Returns
        -------
        Callable
            The function which returns the array.
        """
        Snap._array_registry[fn.__name__] = fn
        return fn

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
        self.sinks = Sinks()
        self._arrays = {}
        self._file_pointer = None
        self._num_particles = 0
        self._families = {key: None for key in Snap._particle_id.keys()}
        self._rotation = None

    def close_file(self):
        """Close access to underlying file."""
        self._file_pointer.close()

    def loaded_arrays(self):
        """Return a list of loaded arrays."""
        return tuple(sorted(self._arrays.keys()))

    def available_arrays(self):
        """Return a list of available arrays."""
        loaded = self.loaded_arrays()
        array_reg = tuple(sorted(self._array_registry.keys()))
        return tuple(sorted(set(loaded + array_reg)))

    @property
    def num_particles(self):
        """Return number of particles."""
        if self._num_particles == 0:
            self._num_particles = self['id'].size
        return self._num_particles

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

        self._rotation = rotation
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
                    self['id'] == Snap._particle_id[name]
                )
            return self._families[name]
        else:
            raise ValueError('Family not available')

    def _get_array(self, name: str, index: Optional[int] = None) -> ndarray:
        """Get an array by name."""
        if name in self._arrays:
            if index is None:
                return self._arrays[name]
            return self._arrays[name][:, index]
        elif name in Snap._array_registry:
            if self._rotation is not None and name in self._rotation_required():
                self._arrays[name] = self._rotation.apply(
                    Snap._array_registry[name](self)
                )
            else:
                self._arrays[name] = Snap._array_registry[name](self)
            if index is None:
                return self._arrays[name]
            return self._arrays[name][:, index]
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
                return self._get_array(self._array_name_mapper[inp])
            elif inp in self._array_split_mapper.keys():
                return self._get_array(*self._array_split_mapper[inp])
            elif inp in self._arrays:
                return self._arrays[inp]
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

    def __setitem__(self, name: str, item: ndarray, force: bool = False):
        """Set an array."""
        if not isinstance(item, ndarray):
            raise ValueError('"item" must be ndarray')
        if item.shape[0] != len(self):
            raise ValueError('Length of array does not match particle number')
        if force:
            self._arrays[name] = item
        else:
            if (
                name in self.available_arrays()
                or name in self._array_split_mapper.keys()
                or name in self._array_name_mapper.keys()
            ):
                raise ValueError(
                    'Attempting to overwrite array. Set force=True to continue.'
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

    def _get_array(self, name: str, index: Optional[int] = None):
        """Get an array by name."""
        if name in self.base._arrays:
            if index is None:
                return self.base._arrays[name][self._indices]
            return self.base._arrays[name][:, index][self._indices]
        elif name in Snap._array_registry:
            self.base._arrays[name] = Snap._array_registry[name](self)
            return self.base._arrays[name][self._indices]
        else:
            raise ValueError('Array not available')


class Sinks:
    """Sink particles in a Snap."""

    _array_name_mapper = {
        'xyz': 'position',
        'pos': 'position',
        'vxyz': 'velocity',
        'vel': 'velocity',
        'h': 'smooth',
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
        'sx': ('spin', 0),
        'sy': ('spin', 1),
        'sz': ('spin', 2),
        'spinx': ('spin', 0),
        'spiny': ('spin', 1),
        'spinz': ('spin', 2),
    }

    def __init__(self):
        self._data = None

    def add_sinks(self, structured_array: ndarray) -> None:
        """Add sinks via structured array.

        Parameters
        ----------
        structured_array
            A structured ndarray with labels such as 'position',
            'velocity', and so on, representing quantities on the sink
            particles.
        """
        self._data = structured_array

    @property
    def columns(self) -> Tuple[str, ...]:
        """Available sink quantities."""
        return self._data.dtype.names

    def __getitem__(self, inp: Union[str, int, slice, List[int]]) -> ndarray:
        """Return an array."""
        if isinstance(inp, (int, slice, list)):
            return self._data[inp]
        elif isinstance(inp, str):
            if inp in self.columns:
                return self._data[inp]
            elif inp in self._array_name_mapper:
                return self._data[self._array_name_mapper[inp]]
            elif inp in self._array_split_mapper:
                array, index = self._array_split_mapper[inp]
                return self._data[array][:, index]
        raise ValueError('Cannot determine quantity to return')

    def __len__(self):
        """Dunder len method."""
        return len(self._data)

    def __repr__(self):
        """Dunder repr method."""
        return self.__str__()

    def __str__(self):
        """Dunder str method."""
        return f'<plonk.snap.Sinks>'
