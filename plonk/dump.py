"""
dump.py

Daniel Mentiplay, 2019.
"""

import collections
from pathlib import Path

import h5py
import numpy as np

FileTypes = collections.namedtuple('FileTypes', 'filetype extension')
FILE_TYPES = [FileTypes(filetype='HDF5', extension='h5')]


class DumpFile:
    def __init__(self, filename):

        if not isinstance(filename, str):
            raise TypeError('filename must be str')
        path = Path(filename)
        self._file_path = path.resolve()
        self._file_name = path.name
        self._file_extension = path.suffix[1:]

        for ft in FILE_TYPES:
            if self._file_extension == ft.extension:
                self._file_type = ft.filetype

        self._open_file()

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        return str(self._file_handle)

    def _open_file(self):

        if not self._file_path.is_file():
            raise FileNotFoundError('Cannot find dump file')

        if self._file_type not in [ft.filetype for ft in FILE_TYPES]:
            raise ValueError('Unknown file type')
        else:
            if self._file_type == 'HDF5':
                self._file_handle = h5py.File(self._file_path, mode='r')

    def _close_file(self):
        self._file_handle.close()


class Dump(DumpFile):
    """
    Smoothed particle hydrodynamics dump file object.

    Parameters
    ----------
    filename : str
        Path to dump file.
    cache_arrays : bool (False)
        Load arrays into memory, otherwise read from file.

    Examples
    --------
    Reading a dump file into a Dump object.

    >>> file_name = 'dumpfile.ext'
    >>> dump = plonk.Dump(file_name)
    """

    def __init__(self, filename, cache_arrays=None):
        super().__init__(filename)

        self._header = {
            key: val[()] for key, val in self._file_handle['header'].items()
        }

        if cache_arrays is None:
            self._cache_arrays = False
        else:
            if not isinstance(cache_arrays, bool):
                raise TypeError('cache_array must be bool')
            self._cache_arrays = cache_arrays

        if self._cache_arrays:
            self.cache_arrays()
        else:
            self._particles = None
            self._sinks = None

    @property
    def particles(self):
        """
        Particle arrays, e.g. position, velocity, smoothing length.
        """
        if self._cache_arrays:
            if self._particles is None:
                self._load_arrays('particles')
            return self._particles
        else:
            return self._read_arrays('particles')

    @property
    def sinks(self):
        """
        Sink arrays, e.g. mass, position, velocity, spin, mass accreted.
        """
        if self._cache_arrays:
            if self._sinks is None:
                self._load_arrays('sinks')
            return self._sinks
        else:
            return self._read_arrays('sinks')

    @property
    def header(self):
        """
        File header, e.g. units, number of particles, numerical
        parameters.
        """
        return self._header

    def cache_arrays(self):
        self._load_arrays('particles')
        self._load_arrays('sinks')
        self._cache_arrays = True

    def _load_arrays(self, array):
        """Load arrays into memory."""

        _array = '_' + array
        setattr(self, _array, self._read_arrays(array))
        setattr(self, _array + '_loaded', True)

    def _read_arrays(self, array):
        """Read arrays into structured Numpy array."""

        array_handle = self._file_handle[array]

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

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        return str(self._file_handle)
