"""Test data: Phantom dust as separate particles."""

import numpy as np

filename = 'dustseparate_00000.h5'

array_name_map = {
    'divv': 'velocity_divergence',
    'dt': 'timestep',
    'dustfrac': 'dust_to_gas_ratio',
    'h': 'smoothing_length',
    'tstop': 'stopping_time',
    'vxyz': 'velocity',
    'xyz': 'position',
}

mean_array_values = {
    'dt': 27.693,
    'dustfrac': 0.004310665694734535,
    'h': 14.473696,
    'tstop': 426.2011872672413,
}

std_array_values = {
    'vxyz': 0.09247972116933928,
    'xyz': 52.5885856475364,
}

properties = {
    'adiabatic_index': 1.0,
    'dust_method': 'dust as separate sets of particles',
    'equation_of_state': 'locally isothermal disc',
    'grain_density': np.array([3000.0]),
    'grain_size': np.array([0.01]),
    'smoothing_length_factor': 1.0,
    'time': 0.0,
}

available_arrays = [
    'angular_momentum',
    'angular_velocity',
    'azimuthal_angle',
    'density',
    'dust_to_gas_ratio',
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

position_shape = (2000, 3)

num_sinks = 1

len_gas = 1000

length_unit = 149600000000.0

mean_ev_values = {
    'time': 5.918015573791112,
    'energy_kinetic': 0.0005606857456011111,
    'energy_thermal': 0.0037500000000000003,
    'energy_magnetic': 0.0,
    'energy_potential': -0.001125769931522222,
    'energy_total': 0.0031849158140888883,
    'momentum': 1.9756579854757646e-07,
    'angular_momentum': 0.41664812648222216,
    'density_max': 2.5104522733888893e-06,
    'density_average': 4.625551587055556e-08,
    'timestep': 0.31312251713111106,
    'entropy': 1.6873464356e-05,
    'mach_number_rms': 7.227724436144445,
    'velocity_rms': 0.1519984266444444,
    'center_of_mass_x': 1.17666825593873e-06,
    'center_of_mass_y': -8.014586524520609e-07,
    'center_of_mass_z': -8.36395109508012e-08,
    'gas_density_max': 2.5104522733888893e-06,
    'gas_density_average': 9.119164037311111e-08,
    'dust_density_max': 3.164142650588889e-08,
    'dust_density_average': 1.3193913680666666e-09,
}
