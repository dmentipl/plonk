"""Dump class for dump files.

The Dump class contains all information related to a smoothed particle
hydrodynamics simulation dump file.
"""

from typing import Callable, Dict, Union

from numpy import ndarray


class Dump:
    """Smoothed particle hydrodynamics dump object.

    Dump files contain the state of the simulation at a point in time.
    Typical minimum data from a smoothed particle hydrodynamics
    simulation include the particle positions and smoothing length, from
    which the density field can be reconstructed, as well as the
    particle type. In addition, the particle velocities are required to
    restart the simulation.

    Other data stored in the dump file include equation of state, dust,
    and magnetic field information, as well as numerical quantities
    related to time-stepping.
    """

    _array_registry: Dict[str, Callable] = {}

    _array_name_mapper = {
        'xyz': 'position',
        'pos': 'position',
        'vxyz': 'velocity',
        'vel': 'velocity',
        'h': 'smooth',
        'rho': 'density',
    }

    @staticmethod
    def add_array(fn):
        """Add array to Dump."""
        Dump._array_registry[fn.__name__] = fn
        return fn

    def __init__(self):

        self.families = {}
        self.properties = {}
        self._arrays = {}
        self._file_pointer = {}
        self._num_particles = 0

    def loaded_arrays(self):
        """Return a list of loaded arrays."""
        return tuple(sorted(self._arrays.keys()))

    def available_arrays(self):
        """Return a list of available arrays."""
        return tuple(sorted(self._array_registry.keys()))

    def _get_array(self, name: str):
        """Get an array by name."""
        if name in self._arrays:
            return self._arrays[name]

        elif name in Dump._array_registry:
            self._arrays[name] = Dump._array_registry[name](self)
            return self._arrays[name]

        else:
            raise ValueError('Array not available')

    @property
    def num_particles(self):
        """Return number of particles."""
        if self._num_particles == 0:
            self._num_particles = self['smooth'].size
        return self._num_particles

    def __getitem__(self, name: Union[str, ndarray]) -> ndarray:
        """Return an array, or family, or subset."""
        if isinstance(name, str):
            if name in self.families:
                raise NotImplementedError('')
            elif name in self.available_arrays():
                return self._get_array(name)
            elif name in self._array_name_mapper.keys():
                return self._get_array(self._array_name_mapper[name])
            else:
                raise ValueError('Cannot determine item to return')
        elif isinstance(name, ndarray):
            raise NotImplementedError('')

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
        return f'<plonk.Dump>'
