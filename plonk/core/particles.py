"""
particles.py

Daniel Mentiplay, 2019.
"""

import numpy as np

I_GAS = 1
I_DUST = 7


def calculate_extra_quantity(dump, quantity, **kwargs):
    """
    Calculate extra quantity.

    Computes an extra quantity on the dump specified by a string.

    Parameters
    ----------
    dump : plonk.Dump
        The plonk dump object.

    quantity : str
        A string specifying the extra quantity to calculate.

    **kwargs
        Extra arguments to functions to calculate specific quantites.

    Examples
    --------
    Calculating the angular momentum.

    >>> quantity = 'angular momentum'
    >>> calculate_extra_quantity(dump, quantity)
    """

    quantities = [
        'spherical radius',
        'cylindrical radius',
        'velocity magnitude',
        'momentum',
        'momentum magnitude',
        'angular velocity',
        'angular momentum',
        'angular momentum magnitude',
        'specific angular momentum',
        'specific angular momentum magnitude',
    ]

    if quantity not in quantities:
        print(f'{quantity} not available')
        return None

    if quantity in ['spherical radius']:
        data = (dump.particles['xyz'],)
        func = _spherical_radius
        kwargs = {}

    elif quantity in ['cylindrical radius']:
        data = (dump.particles['xyz'],)
        func = _cylindrical_radius
        kwargs = {}

    elif quantity in ['velocity magnitude']:
        data = (dump.particles['vxyz'],)
        func = _velocity_magnitude
        kwargs = {}

    elif quantity in ['momentum']:
        data = dump.particles['vxyz'], dump.particle_mass
        func = _momentum
        kwargs = {}

    elif quantity in ['momentum magnitude']:
        data = dump.particles['vxyz'], dump.particle_mass
        func = _momentum_magnitude
        kwargs = {}

    elif quantity in ['angular velocity']:
        data = dump.particles['xyz'], dump.particles['vxyz']
        func = _angular_velocity
        kwargs = {}

    elif quantity in ['angular momentum']:
        data = dump.particles['xyz'], dump.particles['vxyz'], dump.particle_mass
        func = _angular_momentum
        kwargs = {}

    elif quantity in ['angular momentum magnitude']:
        data = dump.particles['xyz'], dump.particles['vxyz'], dump.particle_mass
        func = _angular_momentum
        kwargs = {}

    elif quantity in ['specific angular momentum']:
        data = dump.particles['xyz'], dump.particles['vxyz']
        func = _specific_angular_momentum
        kwargs = {}

    elif quantity in ['specific angular momentum magnitude']:
        data = dump.particles['xyz'], dump.particles['vxyz']
        func = _specific_angular_momentum_magnitude
        kwargs = {}

    return _call_function_on_data(*data, func=func, **kwargs)


def _call_function_on_data(*data, func, **kwargs):
    return func(*data, **kwargs)


def _spherical_radius(position):
    return np.linalg.norm(position, axis=1)


def _cylindrical_radius(position):
    return np.linalg.norm(position[:, 0:2], axis=1)


def _cylindrical_azimuthal_angle(position):
    return np.arctan2(position[:, 1], position[:, 0])


def _velocity_magnitude(velocity):
    return np.linalg.norm(velocity, axis=1)


def _momentum(velocity, mass):
    if isinstance(mass, float) or isinstance(mass, int):
        return mass * velocity
    if isinstance(mass, np.ndarray):
        if mass.ndim == 0:
            return mass * velocity
        elif mass.ndim == 1:
            return mass[:, np.newaxis] * velocity
    raise ValueError('Check inputs, probably mass')


def _momentum_magnitude(velocity, mass):
    return mass * _velocity_magnitude(velocity)


def _specific_angular_momentum(position, velocity):
    if position.ndim > 1:
        if position.shape[1] != 3:
            raise ValueError('Wrong shape array')
    if position.shape != velocity.shape:
        raise ValueError('Position and velocity array shapes must match')
    return np.cross(position, velocity)


def _specific_angular_momentum_magnitude(position, velocity):
    return np.linalg.norm(
        _specific_angular_momentum(position, velocity), axis=1
    )


def _angular_velocity(position, velocity):
    return _specific_angular_momentum_magnitude(
        position, velocity
    ) / _spherical_radius(position)


def _angular_momentum(position, velocity, mass):
    if isinstance(mass, float) or isinstance(mass, int):
        return mass * _specific_angular_momentum(position, velocity)
    if isinstance(mass, np.ndarray):
        if mass.ndim == 0:
            return mass * _specific_angular_momentum(position, velocity)
        elif mass.ndim == 1:
            return mass[:, np.newaxis] * _specific_angular_momentum(
                position, velocity
            )
    raise ValueError('Check inputs, probably mass')


def _angular_momentum_magnitude(position, velocity, mass):
    return np.linalg.norm(_angular_momentum(position, velocity, mass), axis=1)


def _eccentricity(position, velocity, gravitational_parameter):
    kinetic_energy = 1 / 2 * _velocity_magnitude(velocity) ** 2
    potential_energy = -gravitational_parameter / _spherical_radius(position)
    energy = kinetic_energy + potential_energy
    term = (
        2
        * energy
        * _specific_angular_momentum_magnitude(position, velocity) ** 2
        / gravitational_parameter ** 2
    )
    return np.sqrt(1 + term)
