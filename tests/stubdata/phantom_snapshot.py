"""Stub data for reading Phantom snapshots."""

import numpy as np

array_name_map = {
    'divv': 'velocity_divergence',
    'dt': 'timestep',
    'dustfrac': 'dust_fraction',
    'h': 'smoothing_length',
    'tstop': 'stopping_time',
    'vxyz': 'velocity',
    'xyz': 'position',
}

mean_array_values = {
    'divv': -0.00014556303,
    'dt': 27.693,
    'dustfrac': 0.004310665694734535,
    'h': 14.473696,
    'tstop': 426.2011872672413,
    'vxyz': -1.5469107476443848e-17,
    'xyz': 2.622376390111943e-14,
}

properties = {
    'adiabatic_index': 1.0,
    'grain_density': np.array([3.0]),
    'grain_size': np.array([1.0]),
    'smoothing_length_factor': 1.0,
    'polytropic_constant': 0.0025,
    'sound_speed_index': 0.25,
    'time': 0.0,
}
