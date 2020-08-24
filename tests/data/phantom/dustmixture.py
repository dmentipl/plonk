"""Test data for reading Phantom snapshots."""

import numpy as np

filename = 'dustmixture_00000.h5'

array_name_map = {
    'deltavxyz': 'differential_velocity',
    'divv': 'velocity_divergence',
    'dt': 'timestep',
    'dustfrac': 'dust_fraction',
    'divv': 'velocity_divergence',
    'dt': 'timestep',
    'h': 'smoothing_length',
    'tstop': 'stopping_time',
    'vxyz': 'velocity',
    'xyz': 'position',
}

mean_array_values = {
    'dt': 27.806782,
    'dustfrac': 0.0006161144306891705,
    'h': 14.530762,
    'tstop': 139.21310126271644,
}

std_array_values = {
    'deltavxyz': 0.00472864287869192,
    'vxyz': 0.08860771690319386,
    'xyz': 52.7968040666829,
}

properties = {
    'adiabatic_index': 1.0,
    'dust_method': 'dust/gas mixture',
    'equation_of_state': 'locally isothermal disc',
    'grain_density': np.array([3000.0, 3000.0, 3000.0, 3000.0, 3000.0, 3000.0]),
    'grain_size': np.array(
        [
            2.15443469e-06,
            1.00000000e-05,
            4.64158883e-05,
            2.15443469e-04,
            1.00000000e-03,
            4.64158883e-03,
        ]
    ),
    'smoothing_length_factor': 1.0,
    'time': 0.0,
}

available_arrays = [
    'angular_momentum',
    'angular_velocity',
    'azimuthal_angle',
    'density',
    'differential_velocity',
    'dust_density',
    'dust_fraction',
    'dust_mass',
    'eccentricity',
    'gas_density',
    'gas_fraction',
    'gas_mass',
    'id',
    'inclination',
    'keplerian_frequency',
    'kinetic_energy',
    'mass',
    'momentum',
    'polar_angle',
    'position',
    'pressure',
    'radius_cylindrical',
    'radius_spherical',
    'semi_major_axis',
    'smoothing_length',
    'sound_speed',
    'specific_angular_momentum',
    'specific_kinetic_energy',
    'stokes_number',
    'stopping_time',
    'sub_type',
    'temperature',
    'timestep',
    'type',
    'velocity',
    'velocity_divergence',
    'velocity_radial_cylindrical',
    'velocity_radial_spherical',
]

loaded_arrays = ['angular_momentum', 'mass', 'position', 'smoothing_length', 'velocity']

position_shape = (1000, 3)
