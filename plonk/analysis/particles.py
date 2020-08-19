"""Calculate extra quantities on the particles."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from .._logging import logger
from .._units import Quantity
from .._units import units as plonk_units
from ..utils.math import cross, norm

if TYPE_CHECKING:
    from ..snap.snap import SnapLike

ORIGIN = (0, 0, 0) * plonk_units.au
MOLECULAR_HYDROGEN_WEIGHT = 2.381

# Derived quantities require some arrays already present
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
    'radius_cylindrical': ('position',),
    'radius_spherical': ('position',),
    'semi_major_axis': ('position', 'velocity'),
    'specific_angular_momentum': ('position', 'velocity'),
    'specific_kinetic_energy': 'velocity',
    'stokes_number': ('position', 'stopping_time'),
    'temperature': ('sound_speed',),
    'velocity_radial_cylindrical': ('position', 'velocity'),
    'velocity_radial_spherical': ('position', 'velocity'),
}

# Arrays which represent quantities with x-, y-, z-components in space
vector_arrays = ['angular_momentum', 'momentum', 'specific_angular_momentum']

# Arrays which represent dust quantities with columns per dust species
dust_arrays = ['dust_density', 'dust_fraction', 'dust_mass', 'stokes_number']


def angular_momentum(
    snap: SnapLike, origin: Quantity = None, ignore_accreted: bool = False
) -> Quantity:
    """Calculate the angular momentum.

    Parameters
    ----------
    snap
        The Snap object.
    origin : optional
        The origin around which to compute the angular momentum as a
        Quantity like (x, y, z) * au. Default is (0, 0, 0).
    ignore_accreted : optional
        Ignore accreted particles. Default is False.

    Returns
    -------
    Quantity
        The angular momentum on the particles.
    """
    if ignore_accreted:
        h: Quantity = snap['smoothing_length']
        mass: Quantity = snap['mass'][h > 0]
        pos: Quantity = snap['position'][h > 0]
        vel: Quantity = snap['velocity'][h > 0]
    else:
        mass = snap['mass']
        pos = snap['position']
        vel = snap['velocity']

    origin = snap.translation if snap.translation is not None else ORIGIN
    pos = pos - origin

    return mass[:, np.newaxis] * cross(pos, vel)


def angular_velocity(snap: SnapLike, ignore_accreted: bool = False) -> Quantity:
    """Calculate the angular velocity.

    Parameters
    ----------
    snap
        The Snap object.
    ignore_accreted : optional
        Ignore accreted particles. Default is False.

    Returns
    -------
    Quantity
        The angular velocity on the particles.
    """
    if ignore_accreted:
        h: Quantity = snap['smoothing_length']
        pos: Quantity = snap['position'][h > 0]
        vel: Quantity = snap['velocity'][h > 0]
    else:
        pos = snap['position']
        vel = snap['velocity']

    x, y = pos[:, 0], pos[:, 1]
    vx, vy = vel[:, 0], vel[:, 1]

    vphi = (x * vy - y * vx) / (x ** 2 + y ** 2)

    return vphi


def azimuthal_angle(snap: SnapLike, ignore_accreted: bool = False) -> Quantity:
    """Calculate the azimuthal angle.

    Parameters
    ----------
    snap
        The Snap object.
    ignore_accreted : optional
        Ignore accreted particles. Default is False.

    Returns
    -------
    Quantity
        The azimuthal angle on the particles.
    """
    if ignore_accreted:
        h: Quantity = snap['smoothing_length']
        pos: Quantity = snap['position'][h > 0]
    else:
        pos = snap['position']

    x, y = pos[:, 0], pos[:, 1]
    return np.arctan2(y, x)


def dust_density(snap: SnapLike, ignore_accreted: bool = False) -> Quantity:
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
    Quantity
        The dust density per species on the particles.
    """
    if snap.properties['dust_method'] == 'dust/gas mixture':
        _dust_fraction: Quantity = snap['dust_fraction']

    elif snap.properties['dust_method'] == 'dust as separate sets of particles':
        n_dust = len(snap.properties.get('grain_size', []))
        sub_type = snap['sub_type']
        type_mask = snap['type'] == snap.particle_type['dust']
        _dust_fraction = np.zeros((len(snap), n_dust))
        for idx in range(n_dust):
            _dust_fraction[type_mask & (sub_type == idx + 1), idx] = 1

    else:
        raise ValueError('No dust available')

    density: Quantity = snap['density']
    if ignore_accreted:
        h: Quantity = snap['smoothing_length']
        _dust_mass = density[:, np.newaxis] * _dust_fraction
        return _dust_mass[h > 0]
    return density[:, np.newaxis] * _dust_fraction


def dust_fraction(snap: SnapLike, ignore_accreted: bool = False) -> Quantity:
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
    Quantity
        The gas fraction on the particles.
    """
    if snap.properties['dust_method'] == 'dust/gas mixture':
        _dust_fraction: Quantity = snap['dust_fraction']

    elif snap.properties['dust_method'] == 'dust as separate sets of particles':
        n_dust = len(snap.properties.get('grain_size', []))
        sub_type = snap['sub_type']
        type_mask = snap['type'] == snap.particle_type['dust']
        _dust_fraction = np.zeros((len(snap), n_dust))
        for idx in range(n_dust):
            _dust_fraction[type_mask & (sub_type == idx + 1), idx] = 1

    else:
        raise ValueError('No dust available')

    if ignore_accreted:
        h: Quantity = snap['smoothing_length']
        return _dust_fraction[h > 0]
    return _dust_fraction


def dust_mass(snap: SnapLike, ignore_accreted: bool = False) -> Quantity:
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
    Quantity
        The dust mass per species on the particles.
    """
    if snap.properties['dust_method'] == 'dust/gas mixture':
        _dust_fraction: Quantity = snap['dust_fraction']

    elif snap.properties['dust_method'] == 'dust as separate sets of particles':
        n_dust = len(snap.properties.get('grain_size', []))
        sub_type = snap['sub_type']
        type_mask = snap['type'] == snap.particle_type['dust']
        _dust_fraction = np.zeros((len(snap), n_dust))
        for idx in range(n_dust):
            _dust_fraction[type_mask & (sub_type == idx + 1), idx] = 1

    else:
        raise ValueError('No dust available')

    mass: Quantity = snap['mass']
    if ignore_accreted:
        h: Quantity = snap['smoothing_length']
        _dust_mass = mass[:, np.newaxis] * _dust_fraction
        return _dust_mass[h > 0]
    return mass[:, np.newaxis] * _dust_fraction


def eccentricity(
    snap: SnapLike,
    gravitational_parameter: Quantity = None,
    origin: Quantity = None,
    ignore_accreted: bool = False,
) -> Quantity:
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
        The origin around which to compute the eccentricity as a
        Quantity like (x, y, z) * au. Default is (0, 0, 0).
    ignore_accreted : optional
        Ignore accreted particles. Default is False.

    Returns
    -------
    Quantity
        The eccentricity on the particles.
    """
    if ignore_accreted:
        h: Quantity = snap['smoothing_length']
        pos: Quantity = snap['position'][h > 0]
        vel: Quantity = snap['velocity'][h > 0]
    else:
        pos = snap['position']
        vel = snap['velocity']

    origin = snap.translation if snap.translation is not None else ORIGIN
    pos = pos - origin

    mu = snap.properties.get('gravitational_parameter')
    if mu is None:
        raise ValueError(
            'must pass in gravitational_parameter or '
            'set on Snap with set_gravitational_parameter'
        )

    radius = norm(pos, axis=1)

    _specific_angular_momentum = cross(pos, vel)
    specific_angular_momentum_magnitude = norm(_specific_angular_momentum, axis=1)

    _specific_kinetic_energy = 1 / 2 * norm(vel, axis=1) ** 2
    specific_potential_energy = -mu / radius
    specific_energy = _specific_kinetic_energy + specific_potential_energy

    term = specific_energy * (specific_angular_momentum_magnitude / mu) ** 2
    return np.sqrt(1 + 2 * term)


def gas_density(snap: SnapLike, ignore_accreted: bool = False) -> Quantity:
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
    Quantity
        The gas density on the particles.
    """
    if snap.properties['dust_method'] == 'dust/gas mixture':
        _dust_fraction: Quantity = snap['dust_fraction']
        _gas_fraction = 1 - _dust_fraction.sum(axis=1)

    elif snap.properties['dust_method'] == 'dust as separate sets of particles':
        particle_type = snap['type']
        _gas_fraction = np.zeros(len(snap))
        _gas_fraction[particle_type == snap.particle_type['gas']] = 1

    density: Quantity = snap['density']
    if ignore_accreted:
        h: Quantity = snap['smoothing_length']
        return (density * _gas_fraction)[h > 0]
    return density * _gas_fraction


def gas_fraction(snap: SnapLike, ignore_accreted: bool = False) -> Quantity:
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
    Quantity
        The gas fraction on the particles.
    """
    if snap.properties['dust_method'] == 'dust/gas mixture':
        _dust_fraction: Quantity = snap['dust_fraction']
        _gas_fraction = 1 - _dust_fraction.sum(axis=1)

    elif snap.properties['dust_method'] == 'dust as separate sets of particles':
        particle_type = snap['type']
        _gas_fraction = np.ones(len(snap))
        _gas_fraction[particle_type != snap.particle_type['gas']] = 0

    else:
        _gas_fraction = np.ones(len(snap))

    if ignore_accreted:
        h: Quantity = snap['smoothing_length']
        return _gas_fraction[h > 0]
    return _gas_fraction


def gas_mass(snap: SnapLike, ignore_accreted: bool = False) -> Quantity:
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
    Quantity
        The gas mass on the particles.
    """
    if snap.properties['dust_method'] == 'dust/gas mixture':
        _dust_fraction: Quantity = snap['dust_fraction']
        _gas_fraction = 1 - _dust_fraction.sum(axis=1)

    elif snap.properties['dust_method'] == 'dust as separate sets of particles':
        particle_type = snap['type']
        _gas_fraction = np.zeros(len(snap))
        _gas_fraction[particle_type == snap.particle_type['gas']] = 1

    mass: Quantity = snap['mass']
    if ignore_accreted:
        h: Quantity = snap['smoothing_length']
        return (mass * _gas_fraction)[h > 0]
    return mass * _gas_fraction


def inclination(snap: SnapLike, ignore_accreted: bool = False) -> Quantity:
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
    Quantity
        The inclination on the particles.
    """
    if ignore_accreted:
        h: Quantity = snap['smoothing_length']
        mass: Quantity = snap['mass'][h > 0]
        pos: Quantity = snap['position'][h > 0]
        vel: Quantity = snap['velocity'][h > 0]
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


def keplerian_frequency(
    snap: SnapLike,
    gravitational_parameter: Quantity = None,
    origin: Quantity = None,
    ignore_accreted: bool = False,
) -> Quantity:
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
        The origin around which to compute the Keplerian frequency as a
        Quantity like (x, y, z) * au. Default is (0, 0, 0).
    ignore_accreted : optional
        Ignore accreted particles. Default is False.

    Returns
    -------
    Quantity
        The Keplerian frequency on the particles.
    """
    if ignore_accreted:
        h: Quantity = snap['smoothing_length']
        pos: Quantity = snap['position'][h > 0]
    else:
        pos = snap['position']

    origin = snap.translation if snap.translation is not None else ORIGIN
    pos = pos - origin

    mu = snap.properties.get('gravitational_parameter')
    if mu is None:
        raise ValueError(
            'must pass in gravitational_parameter or '
            'set on Snap with set_gravitational_parameter'
        )

    radius = norm(pos, axis=1)
    return np.sqrt(mu / radius ** 3)


def kinetic_energy(snap: SnapLike, ignore_accreted: bool = False) -> Quantity:
    """Calculate the kinetic energy.

    Parameters
    ----------
    snap
        The Snap object.
    ignore_accreted : optional
        Ignore accreted particles. Default is False.

    Returns
    -------
    Quantity
        The kinetic energy on the particles.
    """
    if ignore_accreted:
        h: Quantity = snap['smoothing_length']
        mass: Quantity = snap['mass'][h > 0]
        vel: Quantity = snap['velocity'][h > 0]
    else:
        mass = snap['mass']
        vel = snap['velocity']

    return 1 / 2 * mass * norm(vel, axis=1) ** 2


def momentum(snap: SnapLike, ignore_accreted: bool = False) -> Quantity:
    """Calculate the momentum.

    Parameters
    ----------
    snap
        The Snap object.
    ignore_accreted : optional
        Ignore accreted particles. Default is False.

    Returns
    -------
    Quantity
        The linear momentum on the particles.
    """
    if ignore_accreted:
        h: Quantity = snap['smoothing_length']
        mass: Quantity = snap['mass'][h > 0]
        vel: Quantity = snap['velocity'][h > 0]
    else:
        mass = snap['mass']
        vel = snap['velocity']

    return mass[:, np.newaxis] * vel


def polar_angle(snap: SnapLike, ignore_accreted: bool = False) -> Quantity:
    """Calculate the polar angle.

    Parameters
    ----------
    snap
        The Snap object.
    ignore_accreted : optional
        Ignore accreted particles. Default is False.

    Returns
    -------
    Quantity
        The azimuthal angle on the particles.
    """
    if ignore_accreted:
        h: Quantity = snap['smoothing_length']
        pos: Quantity = snap['position'][h > 0]
    else:
        pos = snap['position']

    x, y, z = pos[:, 0], pos[:, 1], pos[:, 2]
    r = np.sqrt(x ** 2 + y ** 2 + z ** 2)
    return np.arccos(z / r)


def radius_cylindrical(
    snap: SnapLike, origin: Quantity = None, ignore_accreted: bool = False,
) -> Quantity:
    """Calculate the cylindrical radial distance.

    Parameters
    ----------
    snap
        The Snap object.
    ignore_accreted : optional
        Ignore accreted particles. Default is False.

    Returns
    -------
    Quantity
        The radial distance on the particles.
    """
    if ignore_accreted:
        h: Quantity = snap['smoothing_length']
        pos: Quantity = snap['position'][h > 0]
    else:
        pos = snap['position']

    origin = snap.translation if snap.translation is not None else ORIGIN
    pos = pos - origin
    x, y = pos[:, 0], pos[:, 1]

    return np.sqrt(x ** 2 + y ** 2)


def radius_spherical(
    snap: SnapLike, origin: Quantity = None, ignore_accreted: bool = False,
) -> Quantity:
    """Calculate the spherical radial distance.

    Parameters
    ----------
    snap
        The Snap object.
    ignore_accreted : optional
        Ignore accreted particles. Default is False.

    Returns
    -------
    Quantity
        The radial distance on the particles.
    """
    if ignore_accreted:
        h: Quantity = snap['smoothing_length']
        pos: Quantity = snap['position'][h > 0]
    else:
        pos = snap['position']

    origin = snap.translation if snap.translation is not None else ORIGIN
    pos = pos - origin
    x, y, z = pos[:, 0], pos[:, 1], pos[:, 2]

    return np.sqrt(x ** 2 + y ** 2 + z ** 2)


def semi_major_axis(
    snap: SnapLike,
    gravitational_parameter: Quantity = None,
    origin: Quantity = None,
    ignore_accreted: bool = False,
) -> Quantity:
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
        The origin around which to compute the semi-major axis as a
        Quantity like (x, y, z) * au. Default is (0, 0, 0).
    ignore_accreted : optional
        Ignore accreted particles. Default is False.

    Returns
    -------
    Quantity
        The semi-major axis on the particles.
    """
    if ignore_accreted:
        h: Quantity = snap['smoothing_length']
        pos: Quantity = snap['position'][h > 0]
        vel: Quantity = snap['velocity'][h > 0]
    else:
        pos = snap['position']
        vel = snap['velocity']

    origin = snap.translation if snap.translation is not None else ORIGIN
    pos = pos - origin

    mu = snap.properties.get('gravitational_parameter')
    if mu is None:
        raise ValueError(
            'must pass in gravitational_parameter or '
            'set on Snap with set_gravitational_parameter'
        )

    radius = norm(pos, axis=1)

    _specific_angular_momentum = cross(pos, vel)
    specific_angular_momentum_magnitude = norm(_specific_angular_momentum, axis=1)

    _specific_kinetic_energy = 1 / 2 * norm(vel, axis=1) ** 2
    specific_potential_energy = -mu / radius
    specific_energy = _specific_kinetic_energy + specific_potential_energy

    term = specific_energy * (specific_angular_momentum_magnitude / mu) ** 2

    _eccentricity = np.sqrt(1 + 2 * term)

    return specific_angular_momentum_magnitude ** 2 / (mu * (1 - _eccentricity ** 2))


def specific_angular_momentum(
    snap: SnapLike, origin: Quantity = None, ignore_accreted: bool = False
) -> Quantity:
    """Calculate the specific angular momentum.

    Parameters
    ----------
    snap
        The Snap object.
    origin : optional
        The origin around which to compute the angular momentum as a
        Quantity like (x, y, z) * au. Default is (0, 0, 0).
    ignore_accreted : optional
        Ignore accreted particles. Default is False.

    Returns
    -------
    Quantity
        The specific angular momentum on the particles.
    """
    if ignore_accreted:
        h: Quantity = snap['smoothing_length']
        pos: Quantity = snap['position'][h > 0]
        vel: Quantity = snap['velocity'][h > 0]
    else:
        pos = snap['position']
        vel = snap['velocity']

    origin = snap.translation if snap.translation is not None else ORIGIN
    pos = pos - origin

    return cross(pos, vel)


def specific_kinetic_energy(snap: SnapLike, ignore_accreted: bool = False) -> Quantity:
    """Calculate the specific kinetic energy.

    Parameters
    ----------
    snap
        The Snap object.
    ignore_accreted : optional
        Ignore accreted particles. Default is False.

    Returns
    -------
    Quantity
        The specific kinetic energy on the particles.
    """
    if ignore_accreted:
        h: Quantity = snap['smoothing_length']
        vel: Quantity = snap['velocity'][h > 0]
    else:
        vel = snap['velocity']

    return 1 / 2 * norm(vel, axis=1) ** 2


def stokes_number(
    snap: SnapLike,
    gravitational_parameter: Quantity = None,
    origin: Quantity = None,
    ignore_accreted: bool = False,
) -> Quantity:
    """Calculate the Stokes number.

    Parameters
    ----------
    snap
        The Snap object.
    gravitational_parameter
        The gravitational parameter (mu = G M) as a Pint quantity.
    origin : optional
        The origin around which to compute the Stokes number as a
        Quantity like (x, y, z) * au. Default is (0, 0, 0).
    ignore_accreted : optional
        Ignore accreted particles. Default is False.

    Returns
    -------
    Quantity
        The Stokes number on the particles.
    """
    if ignore_accreted:
        h: Quantity = snap['smoothing_length']
        pos: Quantity = snap['position'][h > 0]
        t_s: Quantity = snap['stopping_time'][h > 0]
    else:
        pos = snap['position']
        t_s = snap['stopping_time']

    origin = snap.translation if snap.translation is not None else ORIGIN
    pos = pos - origin

    mu = snap.properties.get('gravitational_parameter')
    if mu is None:
        raise ValueError(
            'must pass in gravitational_parameter or '
            'set on Snap with set_gravitational_parameter'
        )

    radius = norm(pos, axis=1)
    Omega_k = np.sqrt(mu / radius ** 3)

    Stokes = t_s * Omega_k[:, np.newaxis]
    return Stokes


def temperature(
    snap: SnapLike, molecular_weight: float = None, ignore_accreted: bool = False
) -> Quantity:
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
    Quantity
        The gas temperature on the particles.
    """
    if ignore_accreted:
        h: Quantity = snap['smoothing_length']
        cs: Quantity = snap['sound_speed'][h > 0]
    else:
        cs = snap['sound_speed']

    gamma = snap.properties['adiabatic_index']

    molecular_weight = snap.properties.get('molecular_weight')
    if molecular_weight is None:
        logger.warning(
            'molecular_weight not set, assuming 2.381 for molecular hydrogen\n'
            'use set_molecular_weight to set on Snap'
        )
        molecular_weight = MOLECULAR_HYDROGEN_WEIGHT

    molecular_weight = molecular_weight * plonk_units('gram / mole')
    specific_gas_constant = (plonk_units.R / molecular_weight).to_base_units()

    return cs ** 2 / (gamma * specific_gas_constant)


def velocity_radial_cylindrical(
    snap: SnapLike, ignore_accreted: bool = False
) -> Quantity:
    """Calculate the cylindrical radial velocity.

    Parameters
    ----------
    snap
        The Snap object.
    ignore_accreted : optional
        Ignore accreted particles. Default is False.

    Returns
    -------
    Quantity
        The radial velocity on the particles.
    """
    if ignore_accreted:
        h: Quantity = snap['smoothing_length']
        pos: Quantity = snap['position'][h > 0]
        vel: Quantity = snap['velocity'][h > 0]
    else:
        pos = snap['position']
        vel = snap['velocity']

    x, y = pos[:, 0], pos[:, 1]
    vx, vy = vel[:, 0], vel[:, 1]

    return (x * vx + y * vy) / np.sqrt(x ** 2 + y ** 2)


def velocity_radial_spherical(
    snap: SnapLike, ignore_accreted: bool = False
) -> Quantity:
    """Calculate the spherical radial velocity.

    Parameters
    ----------
    snap
        The Snap object.
    ignore_accreted : optional
        Ignore accreted particles. Default is False.

    Returns
    -------
    Quantity
        The radial velocity on the particles.
    """
    if ignore_accreted:
        h: Quantity = snap['smoothing_length']
        pos: Quantity = snap['position'][h > 0]
        vel: Quantity = snap['velocity'][h > 0]
    else:
        pos = snap['position']
        vel = snap['velocity']

    x, y, z = pos[:, 0], pos[:, 1], pos[:, 2]
    vx, vy, vz = vel[:, 0], vel[:, 1], vel[:, 2]

    return (x * vx + y * vy + z * vz) / np.sqrt(x ** 2 + y ** 2 + z ** 2)
