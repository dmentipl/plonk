"""Calculate extra quantities on the particles."""

from __future__ import annotations

from typing import TYPE_CHECKING, Union

import numpy as np

from .._logging import logger
from .._units import Quantity
from .._units import units as plonk_units
from ..utils.math import cross, norm

if TYPE_CHECKING:
    from ..snap.snap import Sinks, SnapLike

ORIGIN = (0, 0, 0) * plonk_units('meter')
MOLECULAR_HYDROGEN_WEIGHT = 2.381

# Derived quantities require some arrays already present
# Ignoring type, sub_type, (optional) smoothing_length
array_requires = {
    'angular_momentum': ['mass', 'position', 'velocity'],
    'angular_velocity': ['position', 'velocity'],
    'azimuthal_angle': ['position'],
    'dust_density': ['density', 'dust_fraction'],
    'dust_mass': ['dust_fraction', 'mass'],
    'gas_density': ['density', 'dust_fraction'],
    'gas_fraction': ['dust_fraction'],
    'gas_mass': ['dust_fraction', 'mass'],
    'kinetic_energy': ['mass', 'velocity'],
    'momentum': ['mass', 'velocity'],
    'polar_angle': ['position'],
    'radius_cylindrical': ['position'],
    'radius_spherical': ['position'],
    'specific_angular_momentum': ['position', 'velocity'],
    'specific_kinetic_energy': ['velocity'],
    'temperature': ['sound_speed'],
    'velocity_radial_cylindrical': ['position', 'velocity'],
    'velocity_radial_spherical': ['position', 'velocity'],
}

# Arrays which represent quantities with x-, y-, z-components in space
vector_arrays = ['angular_momentum', 'momentum', 'specific_angular_momentum']

# Arrays which represent dust quantities with columns per dust species
dust_arrays = ['dust_density', 'dust_mass', 'stokes_number']


def angular_momentum(
    snap: Union[SnapLike, Sinks], origin: Quantity = None, ignore_accreted: bool = False
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

    if origin is None:
        origin = ORIGIN
    pos = pos - origin

    if pos.ndim == 1:
        return mass * cross(pos, vel)
    return mass[:, np.newaxis] * cross(pos, vel)


def angular_velocity(
    snap: Union[SnapLike, Sinks], origin: Quantity = None, ignore_accreted: bool = False
) -> Quantity:
    """Calculate the angular velocity.

    Parameters
    ----------
    snap
        The Snap object.
    origin : optional
        The origin around which to compute the angular velocity as a
        Quantity like (x, y, z) * au. Default is (0, 0, 0).
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

    if origin is None:
        origin = ORIGIN
    pos = pos - origin

    x, y = pos[:, 0], pos[:, 1]
    vx, vy = vel[:, 0], vel[:, 1]

    vphi = (x * vy - y * vx) / (x ** 2 + y ** 2)

    return vphi


def azimuthal_angle(
    snap: Union[SnapLike, Sinks], origin: Quantity = None, ignore_accreted: bool = False
) -> Quantity:
    """Calculate the azimuthal angle.

    Parameters
    ----------
    snap
        The Snap object.
    origin : optional
        The origin around which to compute the azimuthal angle as a
        Quantity like (x, y, z) * au. Default is (0, 0, 0).
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

    if origin is None:
        origin = ORIGIN
    pos = pos - origin

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
        raise ValueError(
            'Dust method is "dust as separate sets of particles"\n'
            'To get dust_density create a SubSnap with snap.family("dust")\n'
            'and then access "density"'
        )

    else:
        raise ValueError('No dust available')

    density: Quantity = snap['density']
    if ignore_accreted:
        h: Quantity = snap['smoothing_length']
        _dust_mass = density[:, np.newaxis] * _dust_fraction
        return _dust_mass[h > 0]
    return density[:, np.newaxis] * _dust_fraction


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
        raise ValueError(
            'Dust method is "dust as separate sets of particles"\n'
            'To get dust_mass create a SubSnap with snap.family("dust")\n'
            'and then access "mass"'
        )

    else:
        raise ValueError('No dust available')

    mass: Quantity = snap['mass']
    if ignore_accreted:
        h: Quantity = snap['smoothing_length']
        _dust_mass = mass[:, np.newaxis] * _dust_fraction
        return _dust_mass[h > 0]
    return mass[:, np.newaxis] * _dust_fraction


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
        raise ValueError(
            'Dust method is "dust as separate sets of particles"\n'
            'To get gas_density create a SubSnap with snap.family("gas")\n'
            'and then access "density"'
        )

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
        raise ValueError(
            'Dust method is "dust as separate sets of particles"\n'
            'The gas_fraction is 1 on gas particles, and 0 otherwise'
        )

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
        raise ValueError(
            'Dust method is "dust as separate sets of particles"\n'
            'To get gas_density create a SubSnap with snap.family("gas")\n'
            'and then access "mass"'
        )

    mass: Quantity = snap['mass']
    if ignore_accreted:
        h: Quantity = snap['smoothing_length']
        return (mass * _gas_fraction)[h > 0]
    return mass * _gas_fraction


def kinetic_energy(
    snap: Union[SnapLike, Sinks], ignore_accreted: bool = False
) -> Quantity:
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


def momentum(snap: Union[SnapLike, Sinks], ignore_accreted: bool = False) -> Quantity:
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

    if vel.ndim == 1:
        return mass * vel
    return mass[:, np.newaxis] * vel


def polar_angle(
    snap: Union[SnapLike, Sinks], origin: Quantity = None, ignore_accreted: bool = False
) -> Quantity:
    """Calculate the polar angle.

    Parameters
    ----------
    snap
        The Snap object.
    origin : optional
        The origin around which to compute the polar angle as a
        Quantity like (x, y, z) * au. Default is (0, 0, 0).
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

    if origin is None:
        origin = ORIGIN
    pos = pos - origin

    x, y, z = pos[:, 0], pos[:, 1], pos[:, 2]
    r = np.sqrt(x ** 2 + y ** 2 + z ** 2)
    return np.arccos(z / r)


def radius_cylindrical(
    snap: Union[SnapLike, Sinks],
    origin: Quantity = None,
    ignore_accreted: bool = False,
) -> Quantity:
    """Calculate the cylindrical radial distance.

    Parameters
    ----------
    snap
        The Snap object.
    origin : optional
        The origin around which to compute the cylindrical radius as a
        Quantity like (x, y, z) * au. Default is (0, 0, 0).
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

    if origin is None:
        origin = ORIGIN
    pos = pos - origin

    x, y = pos[:, 0], pos[:, 1]

    return np.sqrt(x ** 2 + y ** 2)


def radius_spherical(
    snap: Union[SnapLike, Sinks],
    origin: Quantity = None,
    ignore_accreted: bool = False,
) -> Quantity:
    """Calculate the spherical radial distance.

    Parameters
    ----------
    snap
        The Snap object.
    origin : optional
        The origin around which to compute the spherical radius as a
        Quantity like (x, y, z) * au. Default is (0, 0, 0).
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

    if origin is None:
        origin = ORIGIN
    pos = pos - origin

    x, y, z = pos[:, 0], pos[:, 1], pos[:, 2]

    return np.sqrt(x ** 2 + y ** 2 + z ** 2)


def specific_angular_momentum(
    snap: Union[SnapLike, Sinks], origin: Quantity = None, ignore_accreted: bool = False
) -> Quantity:
    """Calculate the specific angular momentum.

    Parameters
    ----------
    snap
        The Snap object.
    origin : optional
        The origin around which to compute the specific angular momentum
        as a Quantity like (x, y, z) * au. Default is (0, 0, 0).
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

    if origin is None:
        origin = ORIGIN
    pos = pos - origin

    return cross(pos, vel)


def specific_kinetic_energy(
    snap: Union[SnapLike, Sinks], ignore_accreted: bool = False
) -> Quantity:
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
    snap: Union[SnapLike, Sinks], origin: Quantity = None, ignore_accreted: bool = False
) -> Quantity:
    """Calculate the cylindrical radial velocity.

    Parameters
    ----------
    snap
        The Snap object.
    origin : optional
        The origin around which to compute the cylindrical radial
        velocity as a Quantity like (x, y, z) * au. Default is
        (0, 0, 0).
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

    if origin is None:
        origin = ORIGIN
    pos = pos - origin

    x, y = pos[:, 0], pos[:, 1]
    vx, vy = vel[:, 0], vel[:, 1]

    return (x * vx + y * vy) / np.sqrt(x ** 2 + y ** 2)


def velocity_radial_spherical(
    snap: Union[SnapLike, Sinks], origin: Quantity = None, ignore_accreted: bool = False
) -> Quantity:
    """Calculate the spherical radial velocity.

    Parameters
    ----------
    snap
        The Snap object.
    origin : optional
        The origin around which to compute the spherical radial velocity
        as a Quantity like (x, y, z) * au. Default is (0, 0, 0).
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

    if origin is None:
        origin = ORIGIN
    pos = pos - origin

    x, y, z = pos[:, 0], pos[:, 1], pos[:, 2]
    vx, vy, vz = vel[:, 0], vel[:, 1], vel[:, 2]

    return (x * vx + y * vy + z * vz) / np.sqrt(x ** 2 + y ** 2 + z ** 2)
