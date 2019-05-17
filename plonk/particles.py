"""
particles.py

Daniel Mentiplay, 2019.
"""

import numpy as np

I_GAS = 1
I_DUST = 7


def calculate_extra_quantity(particles, quantity):
    if quantity in ['momentum']:
        return
    pass


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
