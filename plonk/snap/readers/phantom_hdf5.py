"""Phantom HDF5 snapshot file."""

from __future__ import annotations

from pathlib import Path
from typing import Callable, Union

import h5py
from numpy import ndarray

from ..snap import Snap


class PhantomHDF5Snap:
    """Phantom HDF5 reader."""

    def __init__(self):
        self.snap = Snap()
        self.snap._file_pointer = None

    def generate_snap_from_file(self, filename: Union[str, Path]) -> Snap:
        """Generate a Snap object from a Phantom HDF5 file.

        Parameters
        ----------
        filename
            The path to the file.

        Returns
        -------
        Snap
            A Snap object.
        """
        self.hdf5_file = _HDF5File(filename)
        self.snap._file_pointer = self.hdf5_file.file_handle

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

        return self.snap

    def _header_to_properties(self, header: dict):

        self.snap.properties['time'] = header['time']

        self.snap.properties['udist'] = header['udist']
        self.snap.properties['utime'] = header['utime']
        self.snap.properties['umass'] = header['umass']

        self.snap.properties['hfact'] = header['hfact']

        self.snap.properties['ieos'] = header['ieos']
        self.snap.properties['gamma'] = header['gamma']
        self.snap.properties['polyk'] = 2 / 3 * header['RK2']
        self.snap.properties['qfacdisc'] = header['qfacdisc']

        n_dust = header['ndustsmall'] + header['ndustlarge']
        if n_dust > 0:
            self.snap.properties['grain size'] = header['grainsize'][:n_dust]
            self.snap.properties['grain density'] = header['graindens'][:n_dust]

    def _header_to_mass(self):

        return self._header['massoftype'][_get_dataset('itype') - 1]

    def _populate_array_registry(self):

        arrays = list(self.hdf5_file.file_handle['particles'])

        self.snap._array_registry['position'] = _get_dataset('xyz')
        self.snap._array_registry['smooth'] = _get_dataset('h')
        self.snap._array_registry['itype'] = _get_dataset('itype')

        if 'vxyz' in self.snap._file_pointer['particles']:
            self.snap._array_registry['velocity'] = _get_dataset('vxyz')
            arrays.remove('vxyz')

        self.snap._array_registry['mass'] = _mass
        self.snap._array_registry['density'] = _density

        for array in ('xyz', 'itype', 'h'):
            arrays.remove(array)

        for array in arrays:
            self.snap._array_registry[array] = _get_dataset(array)


def _get_dataset(name: str) -> Callable:
    def func(snap: Snap) -> ndarray:
        return snap._file_pointer[f'particles/{name}'][:]

    return func


def _mass(snap: Snap) -> ndarray:
    massoftype = snap._file_pointer['header/massoftype'][:]
    itype = _get_dataset('itype')(snap)
    return massoftype[itype - 1]


def _density(snap: Snap) -> ndarray:
    m = _mass(snap)
    h = _get_dataset('h')(snap)
    hfact = snap.properties['hfact']
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
            raise FileNotFoundError('Cannot find snapshot file')
        return h5py.File(self.file_path, mode='r')

    def close_file(self):
        self.file_handle.close()
