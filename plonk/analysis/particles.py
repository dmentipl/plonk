"""Calculate quantities on the particles.

The following functions are available:

- angular_momentum
- angular_velocity
- azimuthal_angle
- dust_density
- dust_fraction
- dust_mass
- eccentricity
- gas_density
- gas_fraction
- gas_mass
- inclination
- keplerian_frequency
- kinetic_energy
- momentum
- polar_angle
- radial_distance
- radial_velocity
- semi_major_axis
- specific_angular_momentum
- specific_kinetic_energy
- stokes_number
- temperature
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any, Tuple, Union

import numpy as np
from numpy import ndarray

from .._units import units as plonk_units
from ..utils.math import cross, norm

if TYPE_CHECKING:
    from ..snap.snap import SnapLike

array_units = {
    'angular_momentum': 'angular_momentum',
    'angular_velocity': 'velocity',
    'azimuthal_angle': 'radian',
    'dust_density': 'density',
    'dust_fraction': 'dimensionless',
    'dust_mass': 'mass',
    'eccentricity': 'dimensionless',
    'gas_density': 'density',
    'gas_fraction': 'dimensionless',
    'gas_mass': 'mass',
    'inclination': 'radian',
    'keplerian_frequency': 'frequency',
    'kinetic_energy': 'energy',
    'momentum': 'momentum',
    'polar_angle': 'radian',
    'radial_distance': 'length',
    'radial_velocity': 'velocity',
    'semi_major_axis': 'length',
    'specific_angular_momentum': 'specific_angular_momentum',
    'specific_kinetic_energy': 'specific_energy',
    'stokes_number': 'dimensionless',
    'temperature': 'temperature',
}

array_requires = {
    'angular_momentum': ('mass', 'position', 'velocity'),
    'angular_velocity': ('position', 'velocity'),
    'azimuthal_angle': ('position',),
    'dust_density': ('dust_fraction',),
    'dust_fraction': ('dust_fraction',),
    'dust_mass': ('dust_fraction',),
    'eccentricity': ('position', 'velocity'),
    'gas_density': ('dust_fraction',),
    'gas_fraction': ('dust_fraction',),
    'gas_mass': ('dust_fraction',),
    'inclination': ('mass', 'position', 'velocity'),
    'keplerian_frequency': ('position',),
    'kinetic_energy': ('mass', 'velocity'),
    'momentum': ('mass', 'velocity'),
    'polar_angle': ('position',),
    'radial_distance': ('position',),
    'radial_velocity': ('position', 'velocity'),
    'semi_major_axis': ('position', 'velocity'),
    'specific_angular_momentum': ('position', 'velocity'),
    'specific_kinetic_energy': 'velocity',
    'stokes_number': ('position', 'stopping_time'),
    'temperature': ('sound_speed',),
}

array_rotatable = {
    'angular_momentum': True,
    'angular_velocity': False,
    'azimuthal_angle': False,
    'dust_density': False,
    'dust_fraction': False,
    'dust_mass': False,
    'eccentricity': False,
    'gas_density': False,
    'gas_fraction': False,
    'gas_mass': False,
    'inclination': False,
    'keplerian_frequency': False,
    'kinetic_energy': False,
    'momentum': True,
    'polar_angle': False,
    'radial_distance': False,
    'radial_velocity': False,
    'semi_major_axis': False,
    'specific_angular_momentum': True,
    'specific_kinetic_energy': False,
    'stokes_number': False,
    'temperature': True,
}

array_dust = {
    'angular_momentum': False,
    'angular_velocity': False,
    'azimuthal_angle': False,
    'dust_density': True,
    'dust_fraction': True,
    'dust_mass': True,
    'eccentricity': False,
    'gas_density': False,
    'gas_fraction': False,
    'gas_mass': False,
    'inclination': False,
    'keplerian_frequency': False,
    'kinetic_energy': False,
    'momentum': False,
    'polar_angle': False,
    'radial_distance': False,
    'radial_velocity': False,
    'semi_major_axis': False,
    'specific_angular_momentum': False,
    'specific_kinetic_energy': False,
    'stokes_number': True,
    'temperature': False,
}


def momentum(snap: SnapLike, ignore_accreted: bool = False) -> ndarray:
    """Calculate the momentum.

    Parameters
    ----------
    snap
        The Snap object.
    ignore_accreted : optional
        Ignore accreted particles. Default is False.

    Returns
    -------
    ndarray
        The linear momentum on the particles.
    """
    if ignore_accreted:
        h: ndarray = snap['smoothing_length']
        mass: ndarray = snap['mass'][h > 0]
        vel: ndarray = snap['velocity'][h > 0]
    else:
        mass = snap['mass']
        vel = snap['velocity']

    return mass[:, np.newaxis] * vel


def angular_momentum(
    snap: SnapLike,
    origin: Union[ndarray, Tuple[float, float, float]] = (0.0, 0.0, 0.0),
    ignore_accreted: bool = False,
) -> ndarray:
    """Calculate the angular momentum.

    Parameters
    ----------
    snap
        The Snap object.
    origin : optional
        The origin around which to compute the angular momentum as a
        ndarray or tuple (x, y, z). Default is (0, 0, 0).
    ignore_accreted : optional
        Ignore accreted particles. Default is False.

    Returns
    -------
    ndarray
        The angular momentum on the particles.
    """
    if ignore_accreted:
        h: ndarray = snap['smoothing_length']
        mass: ndarray = snap['mass'][h > 0]
        pos: ndarray = snap['position'][h > 0]
        vel: ndarray = snap['velocity'][h > 0]
    else:
        mass = snap['mass']
        pos = snap['position']
        vel = snap['velocity']

    origin = np.array(origin)
    pos = pos - origin

    return mass[:, np.newaxis] * cross(pos, vel)


def specific_angular_momentum(
    snap: SnapLike,
    origin: Union[ndarray, Tuple[float, float, float]] = (0.0, 0.0, 0.0),
    ignore_accreted: bool = False,
) -> ndarray:
    """Calculate the specific angular momentum.

    Parameters
    ----------
    snap
        The Snap object.
    origin : optional
        The origin around which to compute the angular momentum as a
        ndarray or tuple (x, y, z). Default is (0, 0, 0).
    ignore_accreted : optional
        Ignore accreted particles. Default is False.

    Returns
    -------
    ndarray
        The specific angular momentum on the particles.
    """
    if ignore_accreted:
        h: ndarray = snap['smoothing_length']
        pos: ndarray = snap['position'][h > 0]
        vel: ndarray = snap['velocity'][h > 0]
    else:
        pos = snap['position']
        vel = snap['velocity']

    origin = np.array(origin)
    pos = pos - origin

    return cross(pos, vel)


def kinetic_energy(snap: SnapLike, ignore_accreted: bool = False) -> ndarray:
    """Calculate the kinetic energy.

    Parameters
    ----------
    snap
        The Snap object.
    ignore_accreted : optional
        Ignore accreted particles. Default is False.

    Returns
    -------
    ndarray
        The kinetic energy on the particles.
    """
    if ignore_accreted:
        h: ndarray = snap['smoothing_length']
        mass: ndarray = snap['mass'][h > 0]
        vel: ndarray = snap['velocity'][h > 0]
    else:
        mass = snap['mass']
        vel = snap['velocity']

    return 1 / 2 * mass * norm(vel, axis=1) ** 2


def specific_kinetic_energy(snap: SnapLike, ignore_accreted: bool = False) -> ndarray:
    """Calculate the specific kinetic energy.

    Parameters
    ----------
    snap
        The Snap object.
    ignore_accreted : optional
        Ignore accreted particles. Default is False.

    Returns
    -------
    ndarray
        The specific kinetic energy on the particles.
    """
    if ignore_accreted:
        h: ndarray = snap['smoothing_length']
        vel: ndarray = snap['velocity'][h > 0]
    else:
        vel = snap['velocity']

    return 1 / 2 * norm(vel, axis=1) ** 2


def keplerian_frequency(
    snap: SnapLike,
    gravitational_parameter: Any,
    origin: Union[ndarray, Tuple[float, float, float]] = (0.0, 0.0, 0.0),
    ignore_accreted: bool = False,
) -> ndarray:
    """Calculate the Keplerian orbital frequency.

    The Keplerian orbital frequency of particles around a mass specified
    by gravitational parameter with an optional to specify the position
    of the mass.

    Parameters
    ----------
    snap
        The Snap object.
    gravitational_parameter
        The gravitational parameter (mu = G M) as a Pint quantity.
    origin : optional
        The origin around which to compute the angular momentum as a
        ndarray or tuple (x, y, z). Default is (0, 0, 0).
    ignore_accreted : optional
        Ignore accreted particles. Default is False.

    Returns
    -------
    ndarray
        The Keplerian frequency on the particles.
    """
    if ignore_accreted:
        h: ndarray = snap['smoothing_length']
        pos: ndarray = snap['position'][h > 0]
    else:
        pos = snap['position']

    origin = np.array(origin)
    pos = pos - origin

    mu = gravitational_parameter

    radius = norm(pos, axis=1)
    return np.sqrt(mu / radius ** 3)


def stokes_number(
    snap: SnapLike,
    gravitational_parameter: Any,
    origin: Union[ndarray, Tuple[float, float, float]] = (0.0, 0.0, 0.0),
    ignore_accreted: bool = False,
) -> ndarray:
    """Calculate the Stokes number.

    Parameters
    ----------
    snap
        The Snap object.
    gravitational_parameter
        The gravitational parameter (mu = G M) as a Pint quantity.
    origin : optional
        The origin around which to compute the angular momentum as a
        ndarray or tuple (x, y, z). Default is (0, 0, 0).
    ignore_accreted : optional
        Ignore accreted particles. Default is False.

    Returns
    -------
    ndarray
        The Stokes number on the particles.
    """
    if ignore_accreted:
        h: ndarray = snap['smoothing_length']
        pos: ndarray = snap['position'][h > 0]
        t_s: ndarray = snap['stopping_time'][h > 0]
    else:
        pos = snap['position']
        t_s = snap['stopping_time']

    origin = np.array(origin)
    pos = pos - origin

    mu = gravitational_parameter

    radius = norm(pos, axis=1)
    Omega_k = np.sqrt(mu / radius ** 3)

    Stokes = t_s * Omega_k[:, np.newaxis]
    return Stokes


def semi_major_axis(
    snap: SnapLike,
    gravitational_parameter: Any,
    origin: Union[ndarray, Tuple[float, float, float]] = (0.0, 0.0, 0.0),
    ignore_accreted: bool = False,
) -> ndarray:
    """Calculate the semi-major axis.

    The semi-major axis of particles around a mass specified by
    gravitational parameter with an optional to specify the position of
    the mass.

    Parameters
    ----------
    snap
        The Snap object.
    gravitational_parameter
        The gravitational parameter (mu = G M) as a Pint quantity.
    origin : optional
        The origin around which to compute the angular momentum as a
        ndarray or tuple (x, y, z). Default is (0, 0, 0).
    ignore_accreted : optional
        Ignore accreted particles. Default is False.

    Returns
    -------
    ndarray
        The semi-major axis on the particles.
    """
    if ignore_accreted:
        h: ndarray = snap['smoothing_length']
        pos: ndarray = snap['position'][h > 0]
        vel: ndarray = snap['velocity'][h > 0]
    else:
        pos = snap['position']
        vel = snap['velocity']

    origin = np.array(origin)
    pos = pos - origin

    mu = gravitational_parameter

    radius = norm(pos, axis=1)

    _specific_angular_momentum = cross(pos, vel)
    specific_angular_momentum_magnitude = norm(_specific_angular_momentum, axis=1)

    _specific_kinetic_energy = 1 / 2 * norm(vel, axis=1) ** 2
    specific_potential_energy = -mu / radius
    specific_energy = _specific_kinetic_energy + specific_potential_energy

    term = specific_energy * (specific_angular_momentum_magnitude / mu) ** 2

    _eccentricity = np.sqrt(1 + 2 * term)

    return specific_angular_momentum_magnitude ** 2 / (mu * (1 - _eccentricity ** 2))


def eccentricity(
    snap: SnapLike,
    gravitational_parameter: Any,
    origin: Union[ndarray, Tuple[float, float, float]] = (0.0, 0.0, 0.0),
    ignore_accreted: bool = False,
) -> ndarray:
    """Calculate the eccentricity.

    The eccentricity of particles around a mass specified by
    gravitational parameter with an optional to specify the position of
    the mass.

    Parameters
    ----------
    snap
        The Snap object.
    gravitational_parameter
        The gravitational parameter (mu = G M) as a Pint quantity.
    origin : optional
        The origin around which to compute the angular momentum as a
        ndarray or tuple (x, y, z). Default is (0, 0, 0).
    ignore_accreted : optional
        Ignore accreted particles. Default is False.

    Returns
    -------
    ndarray
        The eccentricity on the particles.
    """
    if ignore_accreted:
        h: ndarray = snap['smoothing_length']
        pos: ndarray = snap['position'][h > 0]
        vel: ndarray = snap['velocity'][h > 0]
    else:
        pos = snap['position']
        vel = snap['velocity']

    origin = np.array(origin)
    pos = pos - origin

    mu = gravitational_parameter

    radius = norm(pos, axis=1)

    _specific_angular_momentum = cross(pos, vel)
    specific_angular_momentum_magnitude = norm(_specific_angular_momentum, axis=1)

    _specific_kinetic_energy = 1 / 2 * norm(vel, axis=1) ** 2
    specific_potential_energy = -mu / radius
    specific_energy = _specific_kinetic_energy + specific_potential_energy

    term = specific_energy * (specific_angular_momentum_magnitude / mu) ** 2
    return np.sqrt(1 + 2 * term)


def inclination(snap: SnapLike, ignore_accreted: bool = False) -> ndarray:
    """Calculate the inclination.

    The inclination is calculated by taking the angle between the
    angular momentum vector and the z-axis, with the angular momentum
    calculated with respect to the center of mass.

    Parameters
    ----------
    snap
        The Snap object.
    ignore_accreted : optional
        Ignore accreted particles. Default is False.

    Returns
    -------
    ndarray
        The inclination on the particles.
    """
    if ignore_accreted:
        h: ndarray = snap['smoothing_length']
        mass: ndarray = snap['mass'][h > 0]
        pos: ndarray = snap['position'][h > 0]
        vel: ndarray = snap['velocity'][h > 0]
    else:
        mass = snap['mass']
        pos = snap['position']
        vel = snap['velocity']

    origin = (mass[:, np.newaxis] * pos).sum(axis=0) / mass.sum()
    pos = pos - origin

    _specific_angular_momentum = cross(pos, vel)

    return np.arccos(
        _specific_angular_momentum[:, 2] / norm(_specific_angular_momentum, axis=1)
    )


def gas_fraction(snap: SnapLike, ignore_accreted: bool = False) -> ndarray:
    """Calculate the gas fraction.

    For dust/gas mixtures this is from the dust fraction.

    Parameters
    ----------
    snap
        The Snap object.
    ignore_accreted : optional
        Ignore accreted particles. Default is False.

    Returns
    -------
    ndarray
        The gas fraction on the particles.
    """
    if snap.properties['dust_method'] == 'dust/gas mixture':
        _dust_fraction: ndarray = snap['dust_fraction']
        _gas_fraction = 1 - _dust_fraction.sum(axis=1)

    elif snap.properties['dust_method'] == 'dust as separate sets of particles':
        particle_type = snap['type']
        _gas_fraction = np.ones(len(snap))
        _gas_fraction[particle_type != snap.particle_type['gas']] = 0

    else:
        _gas_fraction = np.ones(len(snap))

    if ignore_accreted:
        h: ndarray = snap['smoothing_length']
        return _gas_fraction[h > 0]
    return _gas_fraction


def gas_mass(snap: SnapLike, ignore_accreted: bool = False) -> ndarray:
    """Calculate the gas mass.

    For dust/gas mixtures this is from the dust fraction.

    Parameters
    ----------
    snap
        The Snap object.
    ignore_accreted : optional
        Ignore accreted particles. Default is False.

    Returns
    -------
    ndarray
        The gas mass on the particles.
    """
    if snap.properties['dust_method'] == 'dust/gas mixture':
        _dust_fraction: ndarray = snap['dust_fraction']
        _gas_fraction = 1 - _dust_fraction.sum(axis=1)

    elif snap.properties['dust_method'] == 'dust as separate sets of particles':
        particle_type = snap['type']
        _gas_fraction = np.zeros(len(snap))
        _gas_fraction[particle_type == snap.particle_type['gas']] = 1

    mass: ndarray = snap['mass']
    if ignore_accreted:
        h: ndarray = snap['smoothing_length']
        return (mass * _gas_fraction)[h > 0]
    return mass * _gas_fraction


def gas_density(snap: SnapLike, ignore_accreted: bool = False) -> ndarray:
    """Calculate the gas density.

    For dust/gas mixtures this is from the dust fraction.

    Parameters
    ----------
    snap
        The Snap object.
    ignore_accreted : optional
        Ignore accreted particles. Default is False.

    Returns
    -------
    ndarray
        The gas density on the particles.
    """
    if snap.properties['dust_method'] == 'dust/gas mixture':
        _dust_fraction: ndarray = snap['dust_fraction']
        _gas_fraction = 1 - _dust_fraction.sum(axis=1)

    elif snap.properties['dust_method'] == 'dust as separate sets of particles':
        particle_type = snap['type']
        _gas_fraction = np.zeros(len(snap))
        _gas_fraction[particle_type == snap.particle_type['gas']] = 1

    density: ndarray = snap['density']
    if ignore_accreted:
        h: ndarray = snap['smoothing_length']
        return (density * _gas_fraction)[h > 0]
    return density * _gas_fraction


def dust_fraction(snap: SnapLike, ignore_accreted: bool = False) -> ndarray:
    """Calculate the dust fraction.

    For dust/gas mixtures this exists already.

    Parameters
    ----------
    snap
        The Snap object.
    ignore_accreted : optional
        Ignore accreted particles. Default is False.

    Returns
    -------
    ndarray
        The gas fraction on the particles.
    """
    if snap.properties['dust_method'] == 'dust/gas mixture':
        _dust_fraction: ndarray = snap['dust_fraction']

    elif snap.properties['dust_method'] == 'dust as separate sets of particles':
        n_dust = len(snap.properties.get('grain_size', []))
        sub_type = snap['sub_type']
        type_mask = snap['type'] == snap.particle_type['dust']
        _dust_fraction = np.zeros((len(snap), n_dust))
        for idx in range(1, n_dust + 1):
            _dust_fraction[type_mask & (sub_type == idx), idx] = 1

    else:
        raise ValueError('No dust available')

    if ignore_accreted:
        h: ndarray = snap['smoothing_length']
        return _dust_fraction[h > 0]
    return _dust_fraction


def dust_mass(snap: SnapLike, ignore_accreted: bool = False) -> ndarray:
    """Calculate the dust mass per species.

    For dust/gas mixtures this is from the dust fraction.

    Parameters
    ----------
    snap
        The Snap object.
    ignore_accreted : optional
        Ignore accreted particles. Default is False.

    Returns
    -------
    ndarray
        The dust mass per species on the particles.
    """
    if snap.properties['dust_method'] == 'dust/gas mixture':
        _dust_fraction: ndarray = snap['dust_fraction']

    elif snap.properties['dust_method'] == 'dust as separate sets of particles':
        n_dust = len(snap.properties.get('grain_size', []))
        sub_type = snap['sub_type']
        type_mask = snap['type'] == snap.particle_type['dust']
        _dust_fraction = np.zeros((len(snap), n_dust))
        for idx in range(1, n_dust + 1):
            _dust_fraction[type_mask & (sub_type == idx), idx] = 1

    else:
        raise ValueError('No dust available')

    mass: ndarray = snap['mass']
    if ignore_accreted:
        h: ndarray = snap['smoothing_length']
        _dust_mass = mass[:, np.newaxis] * _dust_fraction
        return _dust_mass[h > 0]
    return mass[:, np.newaxis] * _dust_fraction


def dust_density(snap: SnapLike, ignore_accreted: bool = False) -> ndarray:
    """Calculate the dust density per species.

    For dust/gas mixtures this is from the dust fraction.

    Parameters
    ----------
    snap
        The Snap object.
    ignore_accreted : optional
        Ignore accreted particles. Default is False.

    Returns
    -------
    ndarray
        The dust density per species on the particles.
    """
    if snap.properties['dust_method'] == 'dust/gas mixture':
        _dust_fraction: ndarray = snap['dust_fraction']

    elif snap.properties['dust_method'] == 'dust as separate sets of particles':
        n_dust = len(snap.properties.get('grain_size', []))
        sub_type = snap['sub_type']
        type_mask = snap['type'] == snap.particle_type['dust']
        _dust_fraction = np.zeros((len(snap), n_dust))
        for idx in range(1, n_dust + 1):
            _dust_fraction[type_mask & (sub_type == idx), idx] = 1

    else:
        raise ValueError('No dust available')

    density: ndarray = snap['density']
    if ignore_accreted:
        h: ndarray = snap['smoothing_length']
        _dust_mass = density[:, np.newaxis] * _dust_fraction
        return _dust_mass[h > 0]
    return density[:, np.newaxis] * _dust_fraction


def radial_distance(
    snap: SnapLike, coordinates: str = 'cylindrical', ignore_accreted: bool = False
) -> ndarray:
    """Calculate the radial distance.

    Can compute the radial distance in cylindrical and spherical
    coordinates.

    Parameters
    ----------
    snap
        The Snap object.
    coordinates : optional
        Either 'cylindrical' or 'spherical'. Default is 'cylindrical'.
    ignore_accreted : optional
        Ignore accreted particles. Default is False.

    Returns
    -------
    ndarray
        The radial distance on the particles.
    """
    if ignore_accreted:
        h: ndarray = snap['smoothing_length']
        pos: ndarray = snap['position'][h > 0]
    else:
        pos = snap['position']

    x, y, z = pos[:, 0], pos[:, 1], pos[:, 2]

    if coordinates == 'cylindrical':
        return np.sqrt(x ** 2 + y ** 2)
    if coordinates == 'spherical':
        return np.sqrt(x ** 2 + y ** 2 + z ** 2)
    raise ValueError('Cannot determine coordinates')


def azimuthal_angle(snap: SnapLike, ignore_accreted: bool = False) -> ndarray:
    """Calculate the azimuthal angle.

    Parameters
    ----------
    snap
        The Snap object.
    ignore_accreted : optional
        Ignore accreted particles. Default is False.

    Returns
    -------
    ndarray
        The azimuthal angle on the particles.
    """
    if ignore_accreted:
        h: ndarray = snap['smoothing_length']
        pos: ndarray = snap['position'][h > 0]
    else:
        pos = snap['position']

    x, y = pos[:, 0], pos[:, 1]
    return np.arctan2(y, x)


def polar_angle(snap: SnapLike, ignore_accreted: bool = False) -> ndarray:
    """Calculate the polar angle.

    Parameters
    ----------
    snap
        The Snap object.
    ignore_accreted : optional
        Ignore accreted particles. Default is False.

    Returns
    -------
    ndarray
        The azimuthal angle on the particles.
    """
    if ignore_accreted:
        h: ndarray = snap['smoothing_length']
        pos: ndarray = snap['position'][h > 0]
    else:
        pos = snap['position']

    x, y, z = pos[:, 0], pos[:, 1], pos[:, 2]
    r = np.sqrt(x ** 2 + y ** 2 + z ** 2)
    return np.arccos(z / r)


def radial_velocity(
    snap: SnapLike, coordinates: str = 'cylindrical', ignore_accreted: bool = False
) -> ndarray:
    """Calculate the radial velocity.

    Can compute the radial velocity in cylindrical and spherical
    coordinates.

    Parameters
    ----------
    snap
        The Snap object.
    coordinates : optional
        Either 'cylindrical' or 'spherical'. Default is 'cylindrical'.
    ignore_accreted : optional
        Ignore accreted particles. Default is False.

    Returns
    -------
    ndarray
        The radial velocity on the particles.
    """
    if ignore_accreted:
        h: ndarray = snap['smoothing_length']
        pos: ndarray = snap['position'][h > 0]
        vel: ndarray = snap['velocity'][h > 0]
    else:
        pos = snap['position']
        vel = snap['velocity']

    x, y, z = pos[:, 0], pos[:, 1], pos[:, 2]
    vx, vy, vz = vel[:, 0], vel[:, 1], vel[:, 2]

    if coordinates == 'cylindrical':
        vr = (x * vx + y * vy) / np.sqrt(x ** 2 + y ** 2)
    elif coordinates == 'spherical':
        vr = (x * vx + y * vy + z * vz) / np.sqrt(x ** 2 + y ** 2 + z ** 2)
    else:
        raise ValueError('Cannot determine coordinates')

    return vr


def angular_velocity(snap: SnapLike, ignore_accreted: bool = False) -> ndarray:
    """Calculate the angular velocity.

    Parameters
    ----------
    snap
        The Snap object.
    ignore_accreted : optional
        Ignore accreted particles. Default is False.

    Returns
    -------
    ndarray
        The angular velocity on the particles.
    """
    if ignore_accreted:
        h: ndarray = snap['smoothing_length']
        pos: ndarray = snap['position'][h > 0]
        vel: ndarray = snap['velocity'][h > 0]
    else:
        pos = snap['position']
        vel = snap['velocity']

    x, y = pos[:, 0], pos[:, 1]
    vx, vy = vel[:, 0], vel[:, 1]

    vphi = (x * vy - y * vx) / (x ** 2 + y ** 2)

    return vphi


def temperature(
    snap: SnapLike, molecular_weight: float = 2.381, ignore_accreted: bool = False
) -> ndarray:
    """Calculate the gas temperature.

    Parameters
    ----------
    snap
        The Snap object.
    molecular_weight
        The gas molecular weight in gram / mole. E.g. 2.381 for
        molecular hydrogen.
    ignore_accreted : optional
        Ignore accreted particles. Default is False.

    Returns
    -------
    ndarray
        The gas temperature on the particles.
    """
    if ignore_accreted:
        h: ndarray = snap['smoothing_length']
        cs: ndarray = snap['sound_speed'][h > 0]
    else:
        cs = snap['sound_speed']

    gamma = snap.properties['adiabatic_index']

    molecular_weight = molecular_weight * plonk_units('gram / mole')
    specific_gas_constant = (plonk_units.R / molecular_weight).to_base_units()

    return cs ** 2 / (gamma * specific_gas_constant)
