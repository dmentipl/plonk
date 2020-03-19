"""Phantom HDF5 snapshot file."""

from __future__ import annotations

import pathlib
from pathlib import Path
from typing import Callable, Union

import numpy as np
from numpy import ndarray

from ... import units as plonk_units
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

        length = header['udist'] * plonk_units('cm')
        time = header['utime'] * plonk_units('s')
        mass = header['umass'] * plonk_units('g')
        magnetic_field = (
            header['umagfd']
            * plonk_units('g ** (1/2) / cm ** (1/2) / s')
            * np.sqrt(plonk_units.magnetic_constant / (4 * np.pi))
        ).to_base_units()
        units = generate_units_dictionary(length, mass, time, magnetic_field)

        prop = dict()
        prop['time'] = header['time'] * units['time']
        prop['smoothing_length_factor'] = header['hfact']

        prop['adiabatic_index'] = header['gamma']
        prop['polytropic_constant'] = 2 / 3 * header['RK2']

        ieos = header['ieos']
        if ieos == 1:
            prop['equation_of_state'] = 'isothermal'
        elif ieos == 2:
            prop['equation_of_state'] = 'adiabatic'
        elif ieos == 3:
            prop['equation_of_state'] = 'locally isothermal disc'
            prop['sound_speed_index'] = header['qfacdisc']

        ndustsmall = header['ndustsmall']
        ndustlarge = header['ndustlarge']
        if ndustsmall > 0 and ndustlarge > 0:
            raise ValueError(
                'Phantom only supports either dust/gas mixtures (aka 1-fluid dust)\n'
                'or dust as separate sets of particles (aka multi-fluid dust).'
            )
        if ndustsmall > 0:
            prop['dust_method'] = 'dust/gas mixture'
        elif ndustlarge > 0:
            prop['dust_method'] = 'dust as separate sets of particles'

        n_dust = ndustsmall + ndustlarge
        if n_dust > 0:
            prop['grain_size'] = header['grainsize'][:n_dust] * units['length']
            prop['grain_density'] = header['graindens'][:n_dust] * units['density']

        self.snap.properties = prop
        self.snap.units = units

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
        self.snap._array_registry['pressure'] = _pressure
        self.snap._array_registry['sound_speed'] = _sound_speed

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
    hfact = snap.properties['smoothing_length_factor']
    rho = np.zeros(m.shape)
    h = np.abs(h)
    rho[h > 0] = m[h > 0] * (hfact / h[h > 0]) ** 3
    return rho


def _pressure(snap: Snap) -> ndarray:
    ieos = snap._file_pointer['header/ieos'][()]
    K = snap.properties['polytropic_constant']
    gamma = snap.properties['adiabatic_index']
    rho = _density(snap)
    if ieos == 1:
        return K * rho
    elif ieos == 2:
        return K * rho ** (gamma - 1)
    elif ieos == 3:
        q = snap.properties['sound_speed_index']
        pos: ndarray = _get_dataset('xyz', 'particles')(snap)
        r_squared = pos[:, 0] ** 2 + pos[:, 1] ** 2 + pos[:, 2] ** 2
        return K * rho * r_squared ** q
    else:
        raise ValueError('Cannot determine equation of state')


def _sound_speed(snap: Snap) -> ndarray:
    ieos = snap._file_pointer['header/ieos'][()]
    gamma = snap.properties['adiabatic_index']
    rho = _density(snap)
    P = _pressure(snap)
    cs = np.zeros(P.shape)
    if ieos in (1, 3):
        cs[rho > 0] = np.sqrt(P[rho > 0] / rho[rho > 0])
    elif ieos == 2:
        cs[rho > 0] = np.sqrt(gamma * P[rho > 0] / rho[rho > 0])
    else:
        raise ValueError('Cannot determine equation of state')
    return cs
