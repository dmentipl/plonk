"""Phantom HDF5 snapshot file."""

from __future__ import annotations

from pathlib import Path
from typing import Callable, Union

import numpy as np
from numpy import ndarray

from ... import units
from ..snap import Snap
from .hdf5 import HDF5File


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
        self.hdf5_file = HDF5File(filename)
        self.snap._file_pointer = self.hdf5_file.file_handle

        self._header = {
            key: val[()] for key, val in self.hdf5_file.file_handle['header'].items()
        }
        self._header_to_properties(self._header)

        self._populate_array_registry()

        if self._header['nptmass'] > 0:
            self._populate_sink_arrays()

        return self.snap

    def _header_to_properties(self, header: dict):

        self.snap.properties['time'] = header['time']

        self.snap.properties['udist'] = header['udist'] * units.cm
        self.snap.properties['utime'] = header['utime'] * units.s
        self.snap.properties['umass'] = header['umass'] * units.g

        self.snap.properties['hfact'] = header['hfact']

        self.snap.properties['ieos'] = header['ieos']
        self.snap.properties['gamma'] = header['gamma']
        self.snap.properties['polyk'] = 2 / 3 * header['RK2']
        self.snap.properties['qfacdisc'] = header['qfacdisc']

        n_dust = header['ndustsmall'] + header['ndustlarge']
        if n_dust > 0:
            self.snap.properties['grain size'] = header['grainsize'][:n_dust]
            self.snap.properties['grain density'] = header['graindens'][:n_dust]

    def _populate_array_registry(self):

        arrays = list(self.hdf5_file.file_handle['particles'])

        self.snap._array_registry['position'] = _get_dataset('xyz', 'particles')
        self.snap._array_registry['smooth'] = _get_dataset('h', 'particles')

        if 'vxyz' in self.snap._file_pointer['particles']:
            self.snap._array_registry['velocity'] = _get_dataset('vxyz', 'particles')
            arrays.remove('vxyz')

        self.snap._array_registry['id'] = _id
        if self._header['ndustlarge'] > 0:
            self.snap._array_registry['dust_id'] = _dust_id

        self.snap._array_registry['mass'] = _mass
        self.snap._array_registry['density'] = _density

        for array in ('xyz', 'itype', 'h'):
            arrays.remove(array)

        for array in arrays:
            self.snap._array_registry[array] = _get_dataset(array, 'particles')

    def _populate_sink_arrays(self):

        name_map = {
            'xyz': ('position', 'f8', (3,)),
            'vxyz': ('velocity', 'f8', (3,)),
            'm': ('mass', 'f8'),
            'h': ('smooth', 'f8'),
            'hsoft': ('softening', 'f8'),
            'maccreted': ('maccreted', 'f8'),
            'spinxyz': ('spin', 'f8', (3,)),
            'tlast': ('tlast', 'f8'),
        }
        dtype = np.dtype([dt for dt in name_map.values()])
        sinks = np.zeros(self._header['nptmass'], dtype=dtype)

        for name_on_file, array in name_map.items():
            try:
                sinks[array[0]] = self.hdf5_file.file_handle[f'sinks/{name_on_file}'][:]
            except KeyError:
                pass

        self.snap.sinks.add_sinks(sinks)


def _get_dataset(dataset: str, group: str) -> Callable:
    def func(snap: Snap) -> ndarray:
        return snap._file_pointer[f'{group}/{dataset}'][:]

    return func


def _id(snap: Snap) -> ndarray:
    particle_id = np.abs(_get_dataset('itype', 'particles')(snap))
    particle_id[particle_id >= 7] = 2
    return particle_id


def _dust_id(snap: Snap) -> ndarray:
    particle_id = np.abs(_get_dataset('itype', 'particles')(snap))
    dust_id = np.zeros(particle_id.shape, dtype=np.int8)
    dust_id[particle_id >= 7] = particle_id[particle_id >= 7] - 6
    return dust_id


def _mass(snap: Snap) -> ndarray:
    massoftype = snap._file_pointer['header/massoftype'][:]
    particle_id = _get_dataset('itype', 'particles')(snap)
    return massoftype[particle_id - 1]


def _density(snap: Snap) -> ndarray:
    m = _mass(snap)
    h = _get_dataset('h', 'particles')(snap)
    hfact = snap.properties['hfact']
    return m * (hfact / h) ** 3
