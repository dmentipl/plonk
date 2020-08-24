"""Test data: Phantom adiabatic gas shock in a periodic box."""

filename = 'adiabatic_00000.h5'

profile_file = 'adiabatic_00000_profile.csv'

array_name_map = {
    'alpha': 'alpha_viscosity_numerical',
    'divv': 'velocity_divergence',
    'h': 'smoothing_length',
    'u': 'specific_internal_energy',
    'vxyz': 'velocity',
    'xyz': 'position',
}

mean_array_values = {
    'alpha': 0.01663276,
    'divv': 0.0,
    'h': 0.015515282,
    'u': 1.466666659333334,
    'vxyz': 0.0,
    'xyz': -0.0648464781342126,
}


std_array_values = {
    'alpha': 0.10590645,
    'divv': 0.0,
    'h': 0.004334426,
    'u': 0.09428090368680181,
    'vxyz': 0.0,
    'xyz': 0.17061091469222273,
}

properties = {
    'adiabatic_index': 1.66666667,
    'equation_of_state': 'adiabatic',
    'smoothing_length_factor': 1.0,
    'time': 0.0,
}

available_arrays = [
    'alpha_viscosity_numerical',
    'angular_momentum',
    'angular_velocity',
    'azimuthal_angle',
    'density',
    'eccentricity',
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
    'specific_internal_energy',
    'specific_kinetic_energy',
    'sub_type',
    'temperature',
    'type',
    'velocity',
    'velocity_divergence',
    'velocity_radial_cylindrical',
    'velocity_radial_spherical',
]

loaded_arrays = ['angular_momentum', 'mass', 'position', 'smoothing_length', 'velocity']

position_shape = (20736, 3)

num_sinks = 0

len_gas = 16416

length_unit = 0.01
