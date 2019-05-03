"""
disc_00000.py

Test data for reading in Phantom HDF dump file.
"""

import numpy as np


class DumpTest:
    """
    Test Phantom dump.
    """

    def __init__(self):

        self.parameters = {
            'Bextx':
            0.0,
            'Bexty':
            0.0,
            'Bextz':
            0.0,
            'C_cour':
            0.3,
            'C_force':
            0.25,
            'RK2':
            0.0037499999999999994,
            'alpha':
            0.23580036600499324,
            'alphaB':
            1.0,
            'alphau':
            1.0,
            'angtot_in':
            0.4156965032125065,
            'dtmax':
            1154.294847149428,
            'dum':
            0.0,
            'etot_in':
            0.0031034525718330545,
            'fileident':
            b'fulldump: Phantom 1.2.1 7b32a1c (hydro+dust): 26/02/2019 '
            b'09:42:31.2',
            'gamma':
            1.0,
            'get_conserv':
            -1.0,
            'graindens':
            np.array([
                5049628.37866372, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                0., 0., 0., 0., 0., 0., 0., 0., 0., 0.
            ]),
            'grainsize':
            np.array([
                6.68449198e-14, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
                0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
                0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
                0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
                0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
                0.00000000e+00, 0.00000000e+00
            ]),
            'hfact':
            1.0,
            'idust':
            7,
            'ieos':
            3,
            'iexternalforce':
            0,
            'isink':
            0,
            'majorv':
            1,
            'massoftype':
            np.array([
                1.e-07, 0.e+00, 0.e+00, 0.e+00, 0.e+00, 0.e+00, 5.e-09, 0.e+00,
                0.e+00, 0.e+00, 0.e+00, 0.e+00, 0.e+00, 0.e+00, 0.e+00, 0.e+00,
                0.e+00
            ]),
            'mdust_in':
            np.array([
                0.0005, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                0., 0., 0., 0., 0., 0., 0.
            ]),
            'microv':
            1,
            'minorv':
            2,
            'nblocks':
            1,
            'ndustlarge':
            1,
            'ndustsmall':
            0,
            'npartoftype':
            np.array(
                [500000, 0, 0, 0, 0, 0, 100000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                dtype=np.int32),
            'nparttot':
            600000,
            'nptmass':
            1,
            'ntypes':
            17,
            'polyk2':
            0.0,
            'qfacdisc':
            0.25,
            'rhozero':
            0.0,
            'time':
            0.0,
            'tolh':
            0.0001,
            'totmom_in':
            6.812299010933773e-20,
            'udist':
            14960000000000.0,
            'umagfd':
            8138.187949213826,
            'umass':
            1.9891e+33,
            'utime':
            5022728.790082334,
            'xmax':
            0.5,
            'xmin':
            -0.5,
            'ymax':
            0.5,
            'ymin':
            -0.5,
            'zmax':
            0.5,
            'zmin':
            -0.5
        }
