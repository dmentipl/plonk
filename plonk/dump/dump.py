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

    def __getitem__(self, name: Union[str, ndarray]) -> ndarray:
        """Return an array, or family, or subset."""
        if isinstance(name, str):
            if name in self.families:
                raise NotImplementedError('')
            elif name in self.available_arrays():
                return self._get_array(name)
            else:
                raise ValueError('Cannot determine item to return')
        elif isinstance(name, ndarray):
            raise NotImplementedError('')

    def __delitem__(self, name):
        """Delete an array from memory."""
        del self._arrays[name]

    def __repr__(self):
        """Dunder repr method."""
        return self.__str__()

    def __str__(self):
        """Dunder str method."""
        return f'<plonk.Dump>'


########################################################################################
########################################################################################
########################################################################################
########################################################################################

# @property
# def mass(self):
#     """Particle masses."""
#     return self._mass_from_itype()

# @property
# def density(self):
#     """Density on particles from smoothing length."""
#     return self._density_from_smoothing_length(self.header['hfact'])

# @property
# def header(self):
#     """File header.

#     Quantities such as units, number and mass of particles,
#     and numerical parameters.
#     """
#     return self._header

# def _density_from_smoothing_length(self, hfact=1.2):
#     """Calculate density from particle mass and smoothing length."""
#     if self.particles._can_compute_density:
#         return self.mass * (hfact / np.abs(self.particles.arrays['h'])) ** 3
#     else:
#         raise ValueError(f'Cannot compute density on {self.particles}')

# def _mass_from_itype(self):
#     return self.header['massoftype'][self.particles.arrays['itype'][:] - 1]

# def _load_arrays(self, array):
#     """Load arrays into memory."""
#     _array = '_' + array
#     setattr(self, _array, self._read_arrays(array))
#     setattr(self, _array + '_loaded', True)
