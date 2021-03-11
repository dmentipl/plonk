"""Phantom HDF5 snapshot file."""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING, Any, Callable, Dict, List, Tuple, Union

import h5py
import numpy as np

from ..._logging import logger
from ..._units import Quantity
from ..._units import units as plonk_units

if TYPE_CHECKING:
    from ..snap import Snap

igas, iboundary, istar, idarkmatter, ibulge = 1, 3, 4, 5, 6
bignumber = 1e29
missing_infile_parameters = ['alpha', 'alphaB', 'alphau', 'C_cour', 'C_force', 'tolh']


def add_to_header_from_infile(
    snapfile: Union[str, Path],
    infile: Union[str, Path],
    parameters: List[str] = None,
):
    """Add missing values to Phantom snapshot header.

    Some values can be missing from the Phantom snapshot header after
    conversion from the non-HDF5 file format. These values are
    available in the in-file. This function reads the infile and sets
    the values in the Phantom-HDF5 snapshot file header.

    Parameters
    ----------
    snapfile
        The path to the snapshot file.
    infile
        The path to the in-file.
    parameters : optional
        A list of strings of parameters to set in the snap header from
        the in-file.
    """
    try:
        import phantomconfig as pc
    except ModuleNotFoundError as e:
        print(
            'phantom-config unavailable, '
            'install with python -m pip install phantomconfig'
        )
        raise e
    if parameters is None:
        parameters = missing_infile_parameters
    config = pc.read_config(infile)
    logger.info(f'Reading parameters in {infile}')
    filename = Path(snapfile).expanduser()
    if not filename.is_file():
        raise FileNotFoundError('Cannot find snapshot file')
    file = h5py.File(filename, mode='r+')
    logger.info(f'Opening snapshot file {snapfile}')
    for param in parameters:
        if param in config.config:
            logger.info(f'Updating {param} from in-file')
            file[f'header/{param}'][()] = config.config[param].value
    file.close()


def fileid(snap: Snap) -> str:
    """Return Phantom 'fileid'.

    This contains information on the type of snapshot and which physics
    the simulation contains.

    Parameters
    ----------
    snap
        The Snap object.

    Returns
    -------
    str
        The file id.
    """
    return snap._file_pointer['header/fileident'][()].decode()


def snap_properties_and_units(
    file_pointer: h5py.File,
) -> Tuple[Dict[str, Any], Dict[str, Quantity]]:
    """Generate snap properties and units.

    Parameters
    ----------
    file_pointer
        The h5py file pointer to the snap file.

    Returns
    -------
    properties
        The properties as a dict.
    code_units
        The code units as a dict.
    """
    header = {key: val[()] for key, val in file_pointer['header'].items()}
    length = (header['udist'] * plonk_units('cm')).to_base_units()
    time = (header['utime'] * plonk_units('s')).to_base_units()
    mass = (header['umass'] * plonk_units('g')).to_base_units()
    temperature = plonk_units('K')
    magnetic_field = (
        header['umagfd']
        * plonk_units('g ** (1/2) / cm ** (1/2) / s')
        * np.sqrt(plonk_units.magnetic_constant / (4 * np.pi))
    ).to_base_units()
    current = (mass / time ** 2 / magnetic_field).to_base_units()

    code_units = {
        'length': length,
        'time': time,
        'mass': mass,
        'temperature': temperature,
        'current': current,
    }

    properties = dict()

    properties['time'] = header['time'] * code_units['time']
    properties['smoothing_length_factor'] = header['hfact']

    gamma = header['gamma']
    ieos = header['ieos']
    if ieos == 1:
        properties['equation_of_state'] = 'isothermal'
        properties['adiabatic_index'] = gamma
    elif ieos == 2:
        properties['equation_of_state'] = 'adiabatic'
        properties['adiabatic_index'] = gamma
    elif ieos == 3:
        properties['equation_of_state'] = 'locally isothermal disc'
        properties['adiabatic_index'] = gamma

    ndustsmall = header['ndustsmall']
    ndustlarge = header['ndustlarge']
    if ndustsmall > 0 and ndustlarge > 0:
        raise ValueError(
            'Phantom only supports either dust/gas mixtures (aka 1-fluid dust)\n'
            'or dust as separate sets of particles (aka multi-fluid dust).'
        )
    if ndustsmall > 0:
        properties['dust_method'] = 'dust/gas mixture'
    elif ndustlarge > 0:
        properties['dust_method'] = 'dust as separate sets of particles'

    n_dust = ndustsmall + ndustlarge
    if n_dust > 0:
        properties['grain_size'] = header['grainsize'][:n_dust] * code_units['length']
        properties['grain_density'] = (
            header['graindens'][:n_dust]
            * code_units['mass']
            / code_units['length'] ** 3
        )

    return properties, code_units


def snap_array_registry(
    file_pointer: h5py.File, name_map: Dict[str, str] = None
) -> Dict[str, Callable]:
    """Generate snap array registry.

    Parameters
    ----------
    file_pointer
        The h5py file pointer to the snap file.
    name_map : optional
        A dict to convert from Phantom array names to Plonk names.

    Returns
    -------
    Dict
        The particle array registry.
    """
    # The keys are the names of the arrays and the values are functions that return the
    # array when called with snap as the argument.

    if name_map is None:
        name_map = {}

    header = {key: val[()] for key, val in file_pointer['header'].items()}
    arrays = list(file_pointer['particles'])
    ndustsmall = header['ndustsmall']
    ndustlarge = header['ndustlarge']

    array_registry = dict()

    # Each particle gets an id
    array_registry['id'] = particle_id

    # Always read itype, xyz, h
    array_registry['type'] = particle_type
    array_registry['sub_type'] = sub_type
    array_registry['position'] = get_dataset('xyz', 'particles')
    array_registry['smoothing_length'] = get_dataset('h', 'particles')
    arrays.remove('itype')
    arrays.remove('xyz')
    arrays.remove('h')

    # Handle dust
    if ndustsmall > 0:
        array_registry['dust_fraction'] = dust_fraction
        arrays.remove('dustfrac')
    elif ndustlarge > 0:
        array_registry['dust_to_gas_ratio'] = dust_to_gas_ratio
        arrays.remove('dustfrac')
    if ndustsmall > 0 or ndustlarge > 0 and 'tstop' in arrays:
        array_registry['stopping_time'] = stopping_time
        arrays.remove('tstop')

    # Read arrays if available
    for name_on_file, name in name_map.items():
        if name_on_file in arrays:
            array_registry[name] = get_dataset(name_on_file, 'particles')
            arrays.remove(name_on_file)

    # Derived arrays not stored on file
    array_registry['mass'] = mass
    array_registry['density'] = density
    array_registry['pressure'] = pressure
    array_registry['sound_speed'] = sound_speed

    # Read *any* extra arrays
    for array in arrays:
        array_registry[array] = get_dataset(array, 'particles')

    return array_registry


def snap_sink_registry(
    file_pointer: h5py.File, name_map: Dict[str, str] = None
) -> Dict[str, Callable]:
    """Generate snap sink registry.

    Parameters
    ----------
    file_pointer
        The h5py file pointer to the snap file.
    name_map : optional
        A dict to convert from Phantom array names to Plonk names.

    Returns
    -------
    Dict
        The sink array registry.
    """
    # The keys are the names of the arrays and the values are functions that return the
    # array when called with snap as the argument.

    if name_map is None:
        name_map = {}
    header = {key: val[()] for key, val in file_pointer['header'].items()}

    if header['nptmass'] > 0:
        sinks = list(file_pointer['sinks'])
        sink_registry = dict()

        # Read arrays if available
        for name_on_file, name in name_map.items():
            if name_on_file in sinks:
                sink_registry[name] = get_dataset(name_on_file, 'sinks')
                sinks.remove(name_on_file)

        # Read *any* extra arrays
        for sink in sinks:
            sink_registry[sink] = get_dataset(sink, 'sinks')

        return sink_registry

    return {}


def get_dataset(dataset: str, group: str) -> Callable:
    """Return a function that returns an array from file.

    Parameters
    ----------
    dataset
        The name of the HDF5 dataset, i.e. the array name.
    group
        The name of the HDF5 group. For Phantom, this is one of
        'header', 'particles', or 'sinks'.

    Returns
    -------
    Callable
        The function that when called with a Snap object returns the
        array in the dataset with appropriate units, i.e. as a Pint
        Quantity.
    """

    def func(snap: Snap) -> Quantity:
        array = snap._file_pointer[f'{group}/{dataset}'][()]
        name_map = snap._name_map[group]
        if dataset in name_map:
            name = name_map[dataset]
        else:
            name = dataset
        try:
            unit = snap._array_code_units[name]
        except KeyError:
            logger.error(
                f'Cannot get unit of dataset "{group}/{dataset}" - '
                'assuming dimensionless'
            )
            unit = plonk_units('dimensionless')
        return array * unit

    return func


def particle_id(snap: Snap) -> Quantity:
    """Particle id."""
    num_particles = snap._file_pointer['header/nparttot'][()]
    return np.arange(num_particles) * plonk_units('dimensionless')


def particle_type(snap: Snap) -> Quantity:
    """Particle type.

    Type        | Phantom                               | Plonk
    ----------- | ------------------------------------- | -----
    Gas         |                                     1 | 1
    Dust        |      idust -> idust + ndustlarge      | 2
    Boundary    | idustbound -> idustbound + ndustlarge | 3
    Star        |                                     4 | 4
    Dark matter |                                     5 | 5
    Bulge       |                                     6 | 6
    """
    idust = snap._file_pointer['header/idust'][()]
    ndustlarge = snap._file_pointer['header/ndustlarge'][()]
    particle_type = np.abs(get_dataset('itype', 'particles')(snap))
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


def sub_type(snap: Snap) -> Quantity:
    """Particle sub-type.

              | Types
              | Gas | Dust | Boundary
    Sub-types |     |      |
    --------- | --- | ---- | --------
    Gas       |  0  |  n/a |  0
    Dust 1    | n/a |   1  |  1
    Dust 2    | n/a |   2  |  2
    Dust 3    | n/a |   3  |  3
    ...
    """
    particle_type = np.abs(get_dataset('itype', 'particles')(snap))
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
    return sub_type * plonk_units('dimensionless')


def mass(snap: Snap) -> Quantity:
    """Particle mass."""
    massoftype = snap._file_pointer['header/massoftype'][()]
    particle_type = np.array(
        np.abs(get_dataset('itype', 'particles')(snap)).magnitude, dtype=int
    )
    return massoftype[particle_type - 1] * snap._array_code_units['mass']


def density(snap: Snap) -> Quantity:
    """Density."""
    m = (mass(snap) / snap._array_code_units['mass']).magnitude
    h = (
        get_dataset('h', 'particles')(snap) / snap._array_code_units['smoothing_length']
    ).magnitude
    hfact = snap.properties['smoothing_length_factor']
    rho = m * (hfact / np.abs(h)) ** 3
    return rho * snap._array_code_units['density']


def pressure(snap: Snap) -> Quantity:
    """Pressure."""
    ieos = snap._file_pointer['header/ieos'][()]
    K = 2 / 3 * snap._file_pointer['header/RK2'][()]
    gamma = snap.properties['adiabatic_index']
    rho = density(snap)
    if ieos == 1:
        # Globally isothermal
        K = K * snap.code_units['length'] ** 2 * snap.code_units['time'] ** (-2)
        return K * rho
    if ieos == 2:
        # Adiabatic
        try:
            energy = get_dataset('u', 'particles')(snap)
            if gamma > 1.0001:
                return (gamma - 1) * energy * rho
            else:
                return 2 / 3 * energy * rho
        except KeyError:
            K = (
                K
                * snap.code_units['length'] ** (1 - gamma)
                * snap.code_units['mass'] ** (-1 + 3 * gamma)
                * snap.code_units['time'] ** (-2)
            )
            return K * rho ** (gamma - 1)
    if ieos == 3:
        # Vertically isothermal (for accretion disc)
        K = K * snap.code_units['length'] ** 2 * snap.code_units['time'] ** (-2)
        q = snap._file_pointer['header/qfacdisc'][()]
        pos = get_dataset('xyz', 'particles')(snap)
        r_squared = pos[:, 0] ** 2 + pos[:, 1] ** 2 + pos[:, 2] ** 2
        r_squared = (r_squared / snap._array_code_units['position'] ** 2).magnitude
        return K * rho * r_squared ** (-q)
    raise ValueError('Unknown equation of state')


def sound_speed(snap: Snap) -> Quantity:
    """Sound speed."""
    ieos = snap._file_pointer['header/ieos'][()]
    gamma = snap.properties['adiabatic_index']
    rho = density(snap)
    P = pressure(snap)
    if ieos in (1, 3):
        return np.sqrt(P / rho)
    if ieos == 2:
        return np.sqrt(gamma * P / rho)
    raise ValueError('Unknown equation of state')


def stopping_time(snap: Snap) -> Quantity:
    """Dust stopping time."""
    stopping_time = get_dataset('tstop', 'particles')(snap)
    stopping_time[stopping_time == bignumber] = np.inf * snap.code_units['time']
    return stopping_time


def dust_fraction(snap: Snap) -> Quantity:
    """Dust fraction for mixture method (1-fluid)."""
    if snap.properties['dust_method'] != 'dust/gas mixture':
        raise ValueError('Dust fraction only available for "dust/gas mixture"')
    return get_dataset('dustfrac', 'particles')(snap)


def dust_to_gas_ratio(snap: Snap) -> Quantity:
    """Dust-to-gas ratio for separate particles method (2-fluid)."""
    if snap.properties['dust_method'] != 'dust as separate sets of particles':
        raise ValueError(
            'Dust fraction only available for "dust as separate sets of particles"'
        )
    return get_dataset('dustfrac', 'particles')(snap)
