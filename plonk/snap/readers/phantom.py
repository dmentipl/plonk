"""Phantom HDF5 snapshot file."""

from __future__ import annotations

import pathlib
from pathlib import Path
from typing import Callable, Dict, List, Union

import h5py
import numpy as np
from numpy import ndarray

from ... import units as plonk_units
from ..snap import Snap
from ..units import generate_units_dictionary

_bignumber = 1e29

_particle_array_name_map = {
    'abundance': 'abundance',
    'alpha': 'alpha_viscosity_numerical',
    'Bxyz': 'magnetic_field',
    'curlBxyz': 'magnetic_field_curl',
    'deltavxyz': 'differential_velocity',
    'divB': 'magnetic_field_divergence',
    'divv': 'velocity_divergence',
    'dt': 'timestep',
    'dustfrac': 'dust_fraction',
    'eta_AD': 'ambipolar_diffusion_coefficient',
    'eta_HE': 'hall_effect_coefficient',
    'eta_OR': 'ohmic_resistivity_coefficient',
    'graindens': 'grain_density',
    'grainsize': 'grain_size',
    'h': 'smoothing_length',
    'itype': 'type',
    'luminosity': 'luminosity',
    'ne_on_n': 'electron_fraction',
    'poten': 'gravitational_potential',
    'psi': 'magnetic_field_psi',
    'St': 'stokes_number',
    'T': 'temperature',
    'tstop': 'stopping_time',
    'u': 'internal_energy',
    'vrel_on_vfrag': 'relative_on_fragmentation_velocity',
    'vxyz': 'velocity',
    'xyz': 'position',
}

_sink_array_name_map = {
    'xyz': 'position',
    'vxyz': 'velocity',
    'm': 'mass',
    'h': 'accretion_radius',
    'hsoft': 'softening_radius',
    'maccreted': 'mass_accreted',
    'spinxyz': 'spin',
    'tlast': 'last_injection_time',
}


def generate_snap_from_file(filename: Union[str, Path]) -> Snap:
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
    file_path = pathlib.Path(filename).expanduser()
    if not file_path.is_file():
        raise FileNotFoundError('Cannot find snapshot file')
    file_handle = h5py.File(file_path, mode='r')

    snap = Snap()
    snap.data_source = 'Phantom'
    snap.file_path = file_path
    snap._file_pointer = file_handle

    header = {key: val[()] for key, val in file_handle['header'].items()}
    snap.properties, snap.units = _header_to_properties(header)

    arrays = list(file_handle['particles'])
    ndustsmall = header['ndustsmall']
    ndustlarge = header['ndustlarge']
    array_registry = _populate_particle_array_registry(
        arrays=arrays,
        name_map=_particle_array_name_map,
        ndustsmall=ndustsmall,
        ndustlarge=ndustlarge,
    )
    snap._array_registry.update(array_registry)

    if header['nptmass'] > 0:
        sink_registry = _populate_sink_array_registry(name_map=_sink_array_name_map)
        snap._sink_registry.update(sink_registry)

    return snap


def _header_to_properties(header: dict):

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

    return prop, units


def _populate_particle_array_registry(
    arrays: List[str],
    name_map: Dict[str, str],
    ndustsmall: int = 0,
    ndustlarge: int = 0,
):

    array_registry = dict()

    # Always read itype, xyz, h
    array_registry['type'] = _particle_type
    array_registry['position'] = _get_dataset('xyz', 'particles')
    array_registry['smoothing_length'] = _get_dataset('h', 'particles')
    arrays.remove('itype')
    arrays.remove('xyz')
    arrays.remove('h')

    # Handle dust
    if ndustsmall > 0:
        array_registry['dust_fraction'] = _dust_fraction
        arrays.remove('dustfrac')

    elif ndustlarge > 0:
        # Read dust type if there are dust particles
        array_registry['dust_type'] = _dust_particle_type
        array_registry['stopping_time'] = _stopping_time
        array_registry['dust_to_gas_ratio'] = _dust_to_gas_ratio
        arrays.remove('dustfrac')

    # Read arrays if available
    for name_on_file, name in name_map.items():
        if name_on_file in arrays:
            array_registry[name] = _get_dataset(name_on_file, 'particles')
            arrays.remove(name_on_file)

    # Derived arrays not stored on file
    array_registry['mass'] = _mass
    array_registry['density'] = _density
    array_registry['pressure'] = _pressure
    array_registry['sound_speed'] = _sound_speed

    # Read *any* extra arrays
    for array in arrays:
        array_registry[array] = _get_dataset(array, 'particles')

    return array_registry


def _populate_sink_array_registry(name_map: Dict[str, str]):

    sink_registry = dict()

    for name_on_file, name in name_map.items():
        sink_registry[name] = _get_dataset(name_on_file, 'sinks')

    return sink_registry


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
    massoftype = snap._file_pointer['header/massoftype'][()]
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
        return K * rho * r_squared ** (-q)
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


def _stopping_time(snap: Snap) -> ndarray:
    stopping_time = snap._file_pointer['particles/tstop'][()]
    stopping_time[stopping_time == _bignumber] = np.inf
    return stopping_time


def _dust_fraction(snap: Snap) -> ndarray:
    if snap.properties['dust_method'] != 'dust/gas mixture':
        raise ValueError('Dust fraction only available for "dust/gas mixture"')
    dust_fraction = snap._file_pointer['particles/dustfrac'][()]
    return dust_fraction


def _dust_to_gas_ratio(snap: Snap) -> ndarray:
    if snap.properties['dust_method'] != 'dust as separate sets of particles':
        raise ValueError(
            'Dust fraction only available for "dust as separate sets of particles"'
        )
    dust_to_gas_ratio = snap._file_pointer['particles/dustfrac'][()]
    return dust_to_gas_ratio
