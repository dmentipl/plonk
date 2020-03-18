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
    'gamma': 1.0,
    'grain density': np.array([5049628.37866372]),
    'grain size': np.array([6.68449198e-14]),
    'hfact': 1.0,
    'ieos': 3,
    'polyk': 0.0025,
    'qfacdisc': 0.25,
    'time': 0.0,
}
