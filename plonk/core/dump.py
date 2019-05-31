"""
Dump class for dump files.

The Dump class contains all information related to a smoothed particle
hydrodynamics simulation dump file.
"""

import collections
from pathlib import Path

import h5py
import numpy as np

from .constants import constants
from .particles import Arrays
from .units import Units

FileTypes = collections.namedtuple('FileTypes', 'filetype extension')
FILE_TYPES = [FileTypes(filetype='HDF5', extension='h5')]


class DumpFile:
    def __init__(self, filename):

        if not isinstance(filename, str) and not isinstance(filename, Path):
            raise TypeError('filename must be str or pathlib.Path')

        path = Path(filename)
        self.file_path = path.expanduser().resolve()
        self.file_name = path.name
        self.file_extension = path.suffix[1:]

        for ft in FILE_TYPES:
            if self.file_extension == ft.extension:
                self.file_type = ft.filetype

        self._open_file()

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        return str(self.file_handle)

    def _open_file(self):

        if not self.file_path.is_file():
            raise FileNotFoundError('Cannot find dump file')

        if self.file_type not in [ft.filetype for ft in FILE_TYPES]:
            raise ValueError('Unknown file type')
        else:
            if self.file_type == 'HDF5':
                self.file_handle = h5py.File(self.file_path, mode='r')

    def _close_file(self):
        self.file_handle.close()


class Dump(DumpFile):
    """
    Smoothed particle hydrodynamics dump file object.

    Dump files contain the state of the simulation at a point in time.
    Typical minimum data from a smoothed particle hydrodynamics
    simulation include the particle positions and smoothing length, from
    which the density field can be reconstructed, as well as the
    particle type. In addition, the particle velocities are required to
    restart the simulation.

    Other data stored in the dump file include equation of state, dust,
    and magnetic field information, as well as numerical quantities
    related to time-stepping.

    Parameters
    ----------
    filename : str
        Path to dump file.

    cache_arrays : bool, optional (False)
        Load arrays into memory, otherwise read from file.

    Examples
    --------
    Reading a dump file into a Dump object.

    >>> file_name = 'dumpfile.ext'
    >>> dump = plonk.Dump(file_name)

    Accessing the particle arrays object, available particle arrays, and
    particle positions (in two ways).

    >>> dump.particles
    >>> dump.particles.fields
    >>> dump.particles.xyz[:]
    >>> dump.particles.arrays['xyz']

    Accessing the sink arrays object, array data types, and sink spin.

    >>> dump.sinks
    >>> dump.sinks.dtype
    >>> dump.sinks.spinxyz[:]
    >>> dump.sinks.arrays['spinxyz']

    Accessing the dump header dictionary, dump simulation time, and
    particle mass for each type.

    >>> dump.header
    >>> dump.header['time']
    >>> dump.header['massoftype']
    """

    def __init__(self, filename, cache_arrays=None):
        super().__init__(filename)

        self._header = {
            key: val[()] for key, val in self.file_handle['header'].items()
        }

        self.units = Units(
            ulength=self.header['udist'],
            utime=self.header['utime'],
            umass=self.header['umass'],
        )

        self.is_full_dump = self._determine_if_full_dump()

        if cache_arrays is None:
            self._cache_arrays = False
        else:
            if not isinstance(cache_arrays, bool):
                raise TypeError('cache_array must be bool')
            self._cache_arrays = cache_arrays

        self.particles = Arrays(
            arrays_label='particles',
            file_handle=self.file_handle,
            particle_type='fluid',
            cache_arrays=cache_arrays,
        )

        self.sinks = Arrays(arrays_label='sinks', file_handle=self.file_handle)

    @property
    def mass(self):
        """Particle masses"""
        return self._mass_from_itype()

    @property
    def density(self):
        """Density on particles from smoothing length"""
        return self._density_from_smoothing_length(self.header['hfact'])

    @property
    def header(self):
        """
        File header, e.g. units, number of particles, numerical
        parameters.
        """
        return self._header

    @property
    def available_extra_quantities(self):
        """Available extra quantities to calculate on the particles."""
        return self.particles._available_extra_quantities()

    def extra_quantity(self, quantity, sph_type=None, **kwargs):
        """
        Calculate extra quantity.

        Computes an extra quantity on the dump specified by a string.

        Parameters
        ----------
        quantity : str
            A string specifying the extra quantity to calculate.

        sph_type : str, optional (default 'particles')
            The type of SPH data to compute extra quantity on, e.g.
            'particles' or 'sinks'.

        **kwargs
            Extra arguments to functions to calculate specific quantites.

        Returns
        -------
        numpy.ndarray
            An array of the computed extra quantity.

        Examples
        --------
        Calculating the angular momentum two ways.

        >>> dump.extra_quantity('angular momentum')
        >>> dump.extra_quantity('L')
        """

        if sph_type is None:
            sph_type = 'particles'
        if sph_type not in ('particles', 'sinks'):
            raise ValueError('sph_type must be either "particles" or "sinks"')
        if sph_type == 'particles':
            mass = self.mass
        elif sph_type == 'sinks':
            mass = self.sinks.m[:]

        quantities = self.available_extra_quantities

        if quantity not in [element for tupl in quantities for element in tupl]:
            print(f'{quantity} not available')
            return None

        if quantity in ('e', 'eccentricity'):
            gravitational_constant = constants.gravitational_constant / (
                self.header['udist'] ** 3
                / self.header['umass']
                / self.header['utime'] ** 2
            )
            if self.sinks.number == 1:
                stellar_mass = self.sinks.m[0]
            mu = gravitational_constant * stellar_mass
            kwargs = {'gravitational_parameter': mu}

        return getattr(self, sph_type).extra_quantity(
            quantity, mass=mass, **kwargs
        )

    def _density_from_smoothing_length(self, hfact=1.2):
        """Calculate density from particle mass and smoothing length."""

        if self.particles._can_compute_density:
            return self.mass * (hfact / np.abs(self.particles.h)) ** 3
        else:
            print(f'Cannot compute density on {self.particles}')
            return None

    def _mass_from_itype(self):
        return self.header['massoftype'][self.particles.itype[()] - 1]

    def _load_arrays(self, array):
        """Load arrays into memory."""

        _array = '_' + array
        setattr(self, _array, self._read_arrays(array))
        setattr(self, _array + '_loaded', True)

    def _read_arrays(self, array):
        """Read arrays into structured Numpy array."""

        array_handle = self.file_handle[array]

        dtypes = []
        nvals = None
        for key, val in array_handle.items():
            if val.size > 0:
                if nvals is None:
                    nvals = val.shape[0]
                if val.ndim == 1:
                    dtypes.append((key, val.dtype))
                elif val.ndim > 1:
                    dtypes.append((key, val.dtype, val.shape[1:]))

        struct_array = np.zeros(nvals, dtype=dtypes)
        for key in struct_array.dtype.fields:
            struct_array[key] = array_handle[key][()]

        return struct_array

    def _determine_if_full_dump(self):
        # TODO: works for Phantom HDF dumps, maybe not others
        return 'fulldump' in str(self.header['fileident'])

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        return f'<plonk.Dump: "{self.file_name}">'
