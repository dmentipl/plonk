"""Test data: Phantom MHD jet."""

filename = 'mhd_00000.h5'

profile_file = 'mhd_00000_profile.csv'

array_name_map = {
    'Bxyz': 'magnetic_field',
    'alpha': 'alpha_viscosity_numerical',
    'curlBxyz': 'magnetic_field_curl',
    'divB': 'magnetic_field_divergence',
    'divv': 'velocity_divergence',
    'dt': 'timestep',
    'h': 'smoothing_length',
    'poten': 'gravitational_potential',
    'psi': 'magnetic_field_psi',
    'vxyz': 'velocity',
    'xyz': 'position',
}

mean_array_values = {
    'Bxyz': -0.0029753947863998204,
    'alpha': 1.0,
    'divB': -2.7119378e-09,
    'divv': -8.750603e-06,
    'dt': 0.08885766,
    'h': 0.85343355,
    'poten': -4.8283255e-05,
    'psi': 0.0,
    'vxyz': -2.309158043633835e-06,
}

std_array_values = {
    'Bxyz': 0.004207843660662612,
    'alpha': 0.0,
    'curlBxyz': 1.0034642e-07,
    'divB': 1.3255743e-07,
    'divv': 0.0018527617,
    'dt': 0.0,
    'h': 0.49173638,
    'poten': 1.2438333e-05,
    'psi': 0.0,
    'vxyz': 0.018668088990800993,
    'xyz': 3.0172118440124573,
}

properties = {
    'adiabatic_index': 1.0,
    'equation_of_state': 'isothermal',
    'smoothing_length_factor': 1.2,
    'time': 0.0,
}

available_arrays = [
    'alpha_viscosity_numerical',
    'angular_momentum',
    'angular_velocity',
    'azimuthal_angle',
    'density',
    'eccentricity',
    'gravitational_potential',
    'id',
    'inclination',
    'keplerian_frequency',
    'kinetic_energy',
    'magnetic_field',
    'magnetic_field_curl',
    'magnetic_field_divergence',
    'magnetic_field_psi',
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

position_shape = (4843, 3)

num_sinks = 0

len_gas = 4843

length_unit = 100000000000000.0
