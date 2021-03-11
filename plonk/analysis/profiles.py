"""Extra pre-defined profiles."""

import numpy as np

from .._units import Quantity
from .._units import units as plonk_units

G = (1 * plonk_units.newtonian_constant_of_gravitation).to_base_units()

# Some profiles are only appropriate in some geometries (Profile.ndim)
appropriate_ndim = {
    'alpha_shakura_sunyaev': [2],
    'angular_momentum_phi': [2, 3],
    'angular_momentum_theta': [2, 3],
    'aspect_ratio': [2],
    'disc_viscosity': [2],
    'dust_mass': [1, 2, 3],
    'dust_surface_density': [2],
    'epicyclic_frequency': [2],
    'gas_mass': [1, 2, 3],
    'gas_surface_density': [2],
    'linear_density': [1],
    'mass': [1, 2, 3],
    'midplane_stokes_number': [2],
    'scale_height': [2],
    'surface_density': [2],
    'toomre_q': [2],
}

# Dust profiles: can be 'mixture', 'mixture (gas)', 'separate', or 'both'
dust_profiles = {
    'dust_mass': 'mixture',
    'dust_surface_density': 'mixture',
    'gas_mass': 'mixture (gas)',
    'gas_surface_density': 'mixture (gas)',
    'midplane_stokes_number': 'both',
}


def alpha_shakura_sunyaev(prof) -> Quantity:
    """Shakura-Sunyaev alpha disc viscosity."""
    try:
        alpha = prof.snap._file_pointer['header/alpha'][()]
    except KeyError:
        raise ValueError('Cannot determine artificial viscosity alpha')
    return alpha / 10 * prof['smoothing_length'] / prof['scale_height']


def angular_momentum_phi(prof) -> Quantity:
    """Angle between specific angular momentum and x-axis in xy-plane."""
    angular_momentum_x = prof['angular_momentum_x']
    angular_momentum_y = prof['angular_momentum_y']
    return np.arctan2(angular_momentum_y, angular_momentum_x)


def angular_momentum_theta(prof) -> Quantity:
    """Angle between specific angular momentum and xy-plane."""
    angular_momentum_z = prof['angular_momentum_z']
    angular_momentum_magnitude = prof['angular_momentum_mag']
    return np.arccos(angular_momentum_z / angular_momentum_magnitude)


def aspect_ratio(prof) -> Quantity:
    """Aspect ratio profile."""
    H = prof['scale_height']
    R = prof['radius']
    return H / R


def disc_viscosity(prof) -> Quantity:
    """Disc viscosity. I.e. nu = alpha_SS c_s H."""
    return prof['alpha_shakura_sunyaev'] * prof['sound_speed'] * prof['scale_height']


def dust_mass(idx, prof) -> Quantity:
    """Dust mass profile."""
    M = prof.snap[f'dust_mass_{idx+1:03}']
    return prof.particles_to_binned_quantity('sum', M)


def dust_surface_density(idx, prof) -> Quantity:
    """Dust surface density profile.

    Units are [mass / length ** ndim], which depends on ndim of profile.
    """
    return prof[f'dust_mass_{idx+1:03}'] / prof['size']


def epicyclic_frequency(prof) -> Quantity:
    """Epicyclic frequency."""
    Omega = prof['keplerian_frequency']
    R = prof['radius']
    return np.sqrt(2 * Omega / R * np.gradient(R ** 2 * Omega, R))


def gas_mass(prof) -> Quantity:
    """Gas mass profile."""
    M = prof.snap['gas_mass']
    return prof.particles_to_binned_quantity('sum', M)


def gas_surface_density(prof) -> Quantity:
    """Gas surface density profile.

    Units are [mass / length ** ndim], which depends on ndim of profile.
    """
    return prof['gas_mass'] / prof['size']


def linear_density(prof) -> Quantity:
    """Linear density profile.

    Units are [mass / length].
    """
    if prof.ndim != 1:
        raise ValueError('linear_density is appropriate for linear profiles')
    return prof['mass'] / prof['size']


def mass(prof) -> Quantity:
    """Mass profile."""
    M = prof.snap['mass']
    return prof.particles_to_binned_quantity('sum', M)


def midplane_stokes_number(idx, prof) -> Quantity:
    """Midplane Stokes number profile."""
    gamma = prof.snap.properties['adiabatic_index']
    grain_density = prof.snap.properties['grain_density'][idx]
    grain_size = prof.snap.properties['grain_size'][idx]
    return (
        np.pi
        * np.sqrt(gamma)
        * grain_density
        * grain_size
        / 2
        / prof['surface_density']
    )


def scale_height(prof) -> Quantity:
    """Scale height profile."""
    z = prof.snap['z']
    return prof.particles_to_binned_quantity('std', z)


def surface_density(prof) -> Quantity:
    """Surface density profile.

    Units are [mass / length ** 2].
    """
    if prof.ndim != 2:
        raise ValueError('line_density is appropriate for linear profiles')
    return prof['mass'] / prof['size']


def toomre_q(prof) -> Quantity:
    """Toomre Q parameter."""
    return (
        prof['sound_speed']
        * prof['keplerian_frequency']
        / (np.pi * G * prof['surface_density'])
    )
