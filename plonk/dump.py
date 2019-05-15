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

    Examples
    --------
    Reading a dump file into a Dump object.

    >>> file_name = 'dumpfile.ext'
    >>> dump = plonk.Dump(file_name)
    """

    def __init__(self, filename):
        super().__init__(filename)

        self._particles_handle = {
            key: val for key, val in self._file_handle['particles'].items()
        }

        self._sinks_handle = {
            key: val for key, val in self._file_handle['sinks'].items()
        }

        self._header_handle = {
            key: val for key, val in self._file_handle['header'].items()
        }

        self._particles = None
        self._sinks = None

        self._header = {
            key: val[()] for key, val in self._header_handle.items()
        }

        self._particles_loaded = False
        self._sinks_loaded = False

    def _load_arrays(self, array):
        """Load arrays into memory."""

        _array = '_' + array
        setattr(self, _array, self._read_arrays(array))
        setattr(self, _array + '_loaded', True)

    def _read_arrays(self, array):
        """Read arrays into structured Numpy array."""

        _array = '_' + array
        array_handle = getattr(self, _array + '_handle')

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

    @property
    def particles(self):
        """
        Particle arrays, e.g. position, velocity, smoothing length.
        """
        if not self._particles_loaded:
            self._load_arrays('particles')
        return self._particles

    @property
    def sinks(self):
        """
        Sink arrays, e.g. mass, position, velocity, spin, mass accreted.
        """
        if not self._sinks_loaded:
            self._load_arrays('sinks')
        return self._sinks

    @property
    def header(self):
        """
        File header, e.g. units, number of particles, numerical
        parameters.
        """
        return self._header

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        return str(self._file_handle)
