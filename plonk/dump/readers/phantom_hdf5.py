"""Phantom HDF5 dump file."""

from __future__ import annotations

from pathlib import Path
from typing import Callable, Union

import h5py
from numpy import ndarray

from ..dump import Dump


class PhantomHDF5Dump:
    """Phantom HDF5 reader."""

    def __init__(self):
        self.dump = Dump()
        self.dump._file_pointer = None

    def generate_dump_from_file(self, filename: Union[str, Path]) -> Dump:
        """Generate a Dump object from a Phantom HDF5 file.

        Parameters
        ----------
        filename
            The path to the file.

        Returns
        -------
        Dump
            A Dump object.
        """
        self.hdf5_file = _HDF5File(filename)
        self.dump._file_pointer = self.hdf5_file.file_handle

        self._header = {
            key: val[()] for key, val in self.hdf5_file.file_handle['header'].items()
        }

        self._particle_array_handles = {
            field: self.hdf5_file.file_handle[f'particles/{field}']
            for field in self.hdf5_file.file_handle['particles']
        }

        self._sink_array_handles = {
            field: self.hdf5_file.file_handle[f'sinks/{field}']
            for field in self.hdf5_file.file_handle['sinks']
        }

        self._header_to_properties(self._header)
        self._populate_array_registry()

        return self.dump

    def _header_to_properties(self, header: dict):

        self.dump.properties['time'] = header['time']

        self.dump.properties['udist'] = header['udist']
        self.dump.properties['utime'] = header['utime']
        self.dump.properties['umass'] = header['umass']

        self.dump.properties['hfact'] = header['hfact']

        self.dump.properties['ieos'] = header['ieos']
        self.dump.properties['gamma'] = header['gamma']
        self.dump.properties['polyk'] = 2 / 3 * header['RK2']
        self.dump.properties['qfacdisc'] = header['qfacdisc']

        n_dust = header['ndustsmall'] + header['ndustlarge']
        if n_dust > 0:
            self.dump.properties['grain size'] = header['grainsize'][:n_dust]
            self.dump.properties['grain density'] = header['graindens'][:n_dust]

    def _header_to_mass(self):

        return self._header['massoftype'][_get_dataset('itype') - 1]

    def _populate_array_registry(self):

        arrays = list(self.hdf5_file.file_handle['particles'])

        self.dump._array_registry['position'] = _get_dataset('xyz')
        self.dump._array_registry['smooth'] = _get_dataset('h')
        self.dump._array_registry['itype'] = _get_dataset('itype')

        if 'vxyz' in self.dump._file_pointer['particles']:
            self.dump._array_registry['velocity'] = _get_dataset('vxyz')
            arrays.remove('vxyz')

        self.dump._array_registry['mass'] = _mass
        self.dump._array_registry['density'] = _density

        for array in ('xyz', 'itype', 'h'):
            arrays.remove(array)

        for array in arrays:
            self.dump._array_registry[array] = _get_dataset(array)


def _get_dataset(name: str) -> Callable:
    def func(dump: Dump) -> ndarray:
        return dump._file_pointer[f'particles/{name}'][:]

    return func


def _mass(dump: Dump) -> ndarray:
    massoftype = dump._file_pointer['header/massoftype'][:]
    itype = _get_dataset('itype')(dump)
    return massoftype[itype - 1]


def _density(dump: Dump) -> ndarray:
    m = _mass(dump)
    h = _get_dataset('h')(dump)
    hfact = dump.properties['hfact']
    return m * (hfact / h) ** 3


class _HDF5File:
    def __init__(self, filename: Union[str, Path]):
        if isinstance(filename, str):
            path = Path(filename)
        else:
            path = filename

        self.file_path = path.expanduser().resolve()
        self.file_name = path.name
        self.file_extension = path.suffix[1:]
        self.file_type = 'HDF5'
        self.file_handle = self.open_file()

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        return str(self.file_handle)

    def open_file(self):
        if not self.file_path.is_file():
            raise FileNotFoundError('Cannot find dump file')
        return h5py.File(self.file_path, mode='r')

    def close_file(self):
        self.file_handle.close()
