"""
Stub data for Phantom Arrays object.
"""

import numpy as np

dimensions = {
    'divv': 'T^-1',
    'dt': 'T',
    'dustfrac': None,
    'h': 'L',
    'itype': None,
    'tstop': 'T',
    'vxyz': 'L T^-1',
    'xyz': 'L',
}

dtype = {
    'divv': np.dtype('<f4'),
    'dt': np.dtype('<f4'),
    'dustfrac': np.dtype('<f8'),
    'h': np.dtype('<f4'),
    'itype': np.dtype('int8'),
    'tstop': np.dtype('<f8'),
    'vxyz': np.dtype('<f8'),
    'xyz': np.dtype('<f8'),
}

fields = ('divv', 'dt', 'dustfrac', 'h', 'itype', 'tstop', 'vxyz', 'xyz')

number = 2000

shape = {
    'divv': (2000,),
    'dt': (2000,),
    'dustfrac': (2000, 1),
    'h': (2000,),
    'itype': (2000,),
    'tstop': (2000, 1),
    'vxyz': (2000, 3),
    'xyz': (2000, 3),
}

structured_array_dtype = np.dtype(
    [
        ('divv', '<f4'),
        ('dt', '<f4'),
        ('dustfrac', '<f8', (1,)),
        ('h', '<f4'),
        ('itype', 'i1'),
        ('tstop', '<f8', (1,)),
        ('vxyz', '<f8', (3,)),
        ('xyz', '<f8', (3,)),
    ]
)

extra_quantities = {
    '|l|': 8.444529243418579,
    '|v|': 0.13784272262033528,
    'l': 2.762380515521682,
    'R': 79.44477126261133,
    'p': -1.5469107476443848e-17,
    '|L|': 8.444529243418579,
    '|p|': 0.13784272262033528,
    'L': 2.762380515521682,
    'r': 80.40303456661074,
}
