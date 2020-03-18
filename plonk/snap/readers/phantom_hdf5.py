"""Phantom HDF5 snapshot file."""

from __future__ import annotations

import pathlib
from pathlib import Path
from typing import Callable, Union

import numpy as np
from numpy import ndarray

from ... import units
from ..snap import Snap
from ..units import generate_units_dictionary
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
        self.snap.filepath = pathlib.Path(filename).expanduser()
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

        length = header['udist'] * units.cm
        time = header['utime'] * units.s
        mass = header['umass'] * units.g
        self.snap.units = generate_units_dictionary(length, mass, time)

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

        # Always read itype, xyz, h
        self.snap._array_registry['type'] = _particle_type
        self.snap._array_registry['position'] = _get_dataset('xyz', 'particles')
        self.snap._array_registry['smoothing_length'] = _get_dataset('h', 'particles')
        arrays.remove('itype')
        arrays.remove('xyz')
        arrays.remove('h')

        # Read arrays if available
        name_map = {
            'vxyz': 'velocity',
            'u': 'internal_energy',
            'dt': 'timestep',
            'divv': 'velocity_divergence',
            'poten': 'gravitational_potential',
            'tstop': 'stopping_time',
            'dustfrac': 'dust_fraction',
            'deltavxyz': 'differential_velocity',
            'Bxyz': 'magnetic_field',
            'divB': 'magnetic_field_divergence',
            'curlBxyz': 'magnetic_field_curl',
            'psi': 'magnetic_field_psi',
            'eta_OR': 'ohmic_resistivity',
            'eta_HE': 'hall_effect',
            'eta_AD': 'ambipolar_diffusion',
            'ne_on_n': 'electron_density',
        }
        for name_on_file, name in name_map.items():
            if name_on_file in self.snap._file_pointer['particles']:
                self.snap._array_registry[name] = _get_dataset(
                    name_on_file, 'particles'
                )
                arrays.remove(name_on_file)

        # Read dust type if there are dust particles
        if self._header['ndustlarge'] > 0:
            self.snap._array_registry['dust_type'] = _dust_particle_type

        # Derived arrays not stored on file
        self.snap._array_registry['mass'] = _mass
        self.snap._array_registry['density'] = _density

        # Read *any* extra arrays
        for array in arrays:
            self.snap._array_registry[array] = _get_dataset(array, 'particles')

    def _populate_sink_arrays(self):

        name_map = {
            'xyz': 'position',
            'vxyz': 'velocity',
            'm': 'mass',
            'h': 'accretion_radius',
            'hsoft': 'softening_radius',
            'maccreted': 'mass_accreted',
            'spinxyz': 'spin',
            'tlast': 'last_injection_time',
        }

        for name_on_file, name in name_map.items():
            self.snap._sink_registry[name] = _get_dataset(name_on_file, 'sinks')


def _get_dataset(dataset: str, group: str) -> Callable:
    def func(snap: Snap) -> ndarray:
        return snap._file_pointer[f'{group}/{dataset}'][()]

    return func


def _particle_type(snap: Snap) -> ndarray:
    idust = _get_dataset('idust', 'header')(snap)
    particle_type = np.abs(_get_dataset('itype', 'particles')(snap))
    particle_type[particle_type >= idust] = 2
    return particle_type


def _dust_particle_type(snap: Snap) -> ndarray:
    idust = _get_dataset('idust', 'header')(snap)
    particle_type = np.abs(_get_dataset('itype', 'particles')(snap))
    dust_type = np.zeros(particle_type.shape, dtype=np.int8)
    dust_type[particle_type >= idust] = (
        particle_type[particle_type >= idust] - idust + 1
    )
    return dust_type


def _mass(snap: Snap) -> ndarray:
    massoftype = snap._file_pointer['header/massoftype'][:]
    particle_type = _get_dataset('itype', 'particles')(snap)
    return massoftype[particle_type - 1]


def _density(snap: Snap) -> ndarray:
    m = _mass(snap)
    h = _get_dataset('h', 'particles')(snap)
    hfact = snap.properties['hfact']
    return m * (hfact / h) ** 3
