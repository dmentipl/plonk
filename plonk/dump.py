"""
dump.py

Daniel Mentiplay, 2019.
"""

import os

import h5py
import numpy as np


class DumpFile:
    def __init__(self, filename):

        self._file_path = os.path.abspath(filename)
        self._file_name = filename.split('/')[-1]
        self._file_extension = self._file_name.split('.')[-1]

        self._open_file()

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        return str(self._file_handle)

    def _open_file(self):

        if not os.path.isfile(self._file_path):
            raise FileNotFoundError('Cannot find dump file')

        if not self._file_extension == 'h5':
            raise ValueError('Unknown file type')
        else:
            self._file_handle = h5py.File(self._file_path, mode='r')

    def _close_file(self):
        self._file_handle.close()


class Dump(DumpFile):
    """
    Smoothed particle hydrodynamics dump file object.

    Parameters
    ----------
    name : str
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

        self.header = {key: val[()] for key, val in self._header_handle.items()}

        self._particles_loaded = False
        self._sinks_loaded = False

    def _load_arrays(self, array):

        _array = '_' + array
        array_loaded = getattr(self, _array + '_loaded')
        array_handle = getattr(self, _array + '_handle')

        if not array_loaded:

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

            setattr(self, _array, np.zeros(nvals, dtype=dtypes))
            keys = getattr(self, _array).dtype.fields

            for key in keys:
                getattr(self, _array)[key] = array_handle[key][()]

            setattr(self, _array + '_loaded', True)

    @property
    def particles(self):
        """Particle arrays."""
        self._load_arrays('particles')
        return self._particles

    @property
    def sinks(self):
        """Sink arrays."""
        self._load_arrays('sinks')
        return self._sinks

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        return str(self._file_handle)
