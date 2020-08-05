"""Phantom HDF5 snapshot file."""

from __future__ import annotations

import pathlib
from pathlib import Path
from typing import Callable, Dict, List, Union

import h5py
import numpy as np
from numpy import ndarray

from ..._logging import logger
from ..._units import generate_units_array_dictionary, generate_units_dictionary
from ..._units import units as plonk_units
from ..extra import extra_quantities
from ..snap import Snap

igas, iboundary, istar, idarkmatter, ibulge = 1, 3, 4, 5, 6
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
    logger.debug(f'Loading Phantom snapshot: {filename}')
    file_path = pathlib.Path(filename).expanduser()
    if not file_path.is_file():
        raise FileNotFoundError('Cannot find snapshot file')
    file_handle = h5py.File(file_path, mode='r')

    snap = Snap()
    snap.data_source = 'Phantom'
    snap.file_path = file_path
    snap._file_pointer = file_handle

    header = {key: val[()] for key, val in file_handle['header'].items()}
    snap._properties, units = _header_to_properties(header)
    snap._array_units = generate_units_array_dictionary(units)
    snap.units = {
        key: units[key] for key in ['length', 'time', 'mass', 'magnetic_field']
    }

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

    extra_quantities(snap)
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

    gamma = header['gamma']
    ieos = header['ieos']
    if ieos == 1:
        prop['equation_of_state'] = 'isothermal'
        prop['adiabatic_index'] = gamma
    elif ieos == 2:
        prop['equation_of_state'] = 'adiabatic'
        prop['adiabatic_index'] = gamma
    elif ieos == 3:
        prop['equation_of_state'] = 'locally isothermal disc'
        prop['adiabatic_index'] = gamma

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
        array_registry['dust_to_gas_ratio'] = _dust_to_gas_ratio
        arrays.remove('dustfrac')
        # Currently only dust as species has particle sub-types
        array_registry['sub_type'] = _sub_type

    if ndustsmall > 0 or ndustlarge > 0 and 'tstop' in arrays:
        array_registry['stopping_time'] = _stopping_time
        arrays.remove('tstop')

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
        array = snap._file_pointer[f'{group}/{dataset}'][()]
        try:
            unit = snap._array_units[_particle_array_name_map[dataset]]
        except KeyError:
            try:
                unit = snap._array_units[_sink_array_name_map[dataset]]
            except KeyError:
                raise RuntimeError(f'Cannot get unit of {group}/{dataset}')
        return array * unit

    return func


def _particle_type(snap: Snap) -> ndarray:
    # Type        | Phantom                               | Plonk
    # ----------- | ------------------------------------- | -----
    # Gas         |                                     1 | 1
    # Dust        |      idust -> idust + ndustlarge      | 2
    # Boundary    | idustbound -> idustbound + ndustlarge | 3
    # Star        |                                     4 | 4
    # Dark matter |                                     5 | 5
    # Bulge       |                                     6 | 6
    idust = snap._file_pointer['header/idust'][()]
    ndustlarge = snap._file_pointer['header/ndustlarge'][()]
    particle_type = np.abs(_get_dataset('itype', 'particles')(snap))
    particle_type[
        (particle_type >= idust) & (particle_type < idust + ndustlarge)
    ] = snap.particle_type['dust']
    try:
        idustbound = snap._file_pointer['header/idustbound'][()]
        particle_type[
            (particle_type >= idustbound) & (particle_type < idustbound + ndustlarge)
        ] = snap.particle_type['boundary']
    except KeyError:
        if np.any(particle_type >= idust + ndustlarge):
            logger.error('Cannot determine dust boundary particles')
    return np.array(particle_type.magnitude, dtype=int) * plonk_units('dimensionless')


def _sub_type(snap: Snap) -> ndarray:
    #           | Types
    #           | Gas | Dust | Boundary
    # Sub-types |     |      |
    # --------- | --- | ---- | --------
    # Gas       |  0  |  n/a |  0
    # Dust 1    | n/a |   1  |  1
    # Dust 2    | n/a |   2  |  2
    # Dust 3    | n/a |   3  |  3
    # ...
    particle_type = np.abs(_get_dataset('itype', 'particles')(snap))
    sub_type = np.zeros(particle_type.shape, dtype=np.int8)
    sub_type[
        (particle_type == igas)
        | (particle_type == istar)
        | (particle_type == idarkmatter)
        | (particle_type == ibulge)
    ] = 0
    sub_type[particle_type == iboundary] = 0
    idust = snap._file_pointer['header/idust'][()]
    ndustlarge = snap._file_pointer['header/ndustlarge'][()]
    for idx in range(idust, idust + ndustlarge):
        sub_type[particle_type == idx] = idx - idust + 1
    try:
        idustbound = snap._file_pointer['header/idustbound'][()]
        for idx in range(idustbound, idustbound + ndustlarge):
            sub_type[particle_type == idx] = idx - idustbound + 1
    except KeyError:
        if np.any(particle_type >= idust + ndustlarge):
            logger.error('Cannot determine dust boundary particles')
    return sub_type


def _mass(snap: Snap) -> ndarray:
    massoftype = snap._file_pointer['header/massoftype'][()]
    particle_type = np.array(
        np.abs(_get_dataset('itype', 'particles')(snap)).magnitude, dtype=int
    )
    return massoftype[particle_type - 1] * snap._array_units['mass']


def _density(snap: Snap) -> ndarray:
    m = (_mass(snap) / snap._array_units['mass']).magnitude
    h = (
        _get_dataset('h', 'particles')(snap) / snap._array_units['smoothing_length']
    ).magnitude
    hfact = snap.properties['smoothing_length_factor']
    rho = m * (hfact / np.abs(h)) ** 3
    return rho * snap._array_units['density']


def _pressure(snap: Snap) -> ndarray:
    ieos = snap._file_pointer['header/ieos'][()]
    K = 2 / 3 * snap._file_pointer['header/RK2'][()]
    gamma = snap.properties['adiabatic_index']
    rho = _density(snap)
    if ieos == 1:
        # Globally isothermal
        K = K * snap.units['length'] ** 2 * snap.units['time'] ** (-2)
        return K * rho
    if ieos == 2:
        # Adiabatic
        try:
            energy = _get_dataset('u', 'particles')(snap)
            if gamma > 1.0001:
                return (gamma - 1) * energy * rho
            else:
                return 2 / 3 * energy * rho
        except KeyError:
            K = (
                K
                * snap.units['length'] ** (1 - gamma)
                * snap.units['mass'] ** (-1 + 3 * gamma)
                * snap.units['time'] ** (-2)
            )
            return K * rho ** (gamma - 1)
    if ieos == 3:
        # Vertically isothermal (for accretion disc)
        K = K * snap.units['length'] ** 2 * snap.units['time'] ** (-2)
        q = snap._file_pointer['header/qfacdisc'][()]
        pos = _get_dataset('xyz', 'particles')(snap)
        r_squared = pos[:, 0] ** 2 + pos[:, 1] ** 2 + pos[:, 2] ** 2
        r_squared = (r_squared / snap._array_units['position'] ** 2).magnitude
        return K * rho * r_squared ** (-q)
    raise ValueError('Unknown equation of state')


def _sound_speed(snap: Snap) -> ndarray:
    ieos = snap._file_pointer['header/ieos'][()]
    gamma = snap.properties['adiabatic_index']
    rho = _density(snap)
    P = _pressure(snap)
    if ieos in (1, 3):
        return np.sqrt(P / rho)
    if ieos == 2:
        return np.sqrt(gamma * P / rho)
    raise ValueError('Unknown equation of state')


def _stopping_time(snap: Snap) -> ndarray:
    stopping_time = _get_dataset('tstop', 'particles')(snap)
    stopping_time[stopping_time == _bignumber] = np.inf * snap.units['time']
    return stopping_time


def _dust_fraction(snap: Snap) -> ndarray:
    if snap.properties['dust_method'] != 'dust/gas mixture':
        raise ValueError('Dust fraction only available for "dust/gas mixture"')
    return _get_dataset('dustfrac', 'particles')(snap)


def _dust_to_gas_ratio(snap: Snap) -> ndarray:
    if snap.properties['dust_method'] != 'dust as separate sets of particles':
        raise ValueError(
            'Dust fraction only available for "dust as separate sets of particles"'
        )
    return _get_dataset('dustfrac', 'particles')(snap)
