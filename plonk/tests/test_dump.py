'''
test_dump.py
'''

import unittest

import numpy as np

from plonk.dump import Dump
from .test_data.disc_00000 import DumpTest

class TestDump(unittest.TestCase):
    """
    Test reading Phantom dump files.
    """

    def test_read_dump(self):
        """
        Testing reading Phantom HDF5 dump file.
        """

        test_dump = 'tests/test_data/disc_00000.h5'
        dump = Dump(test_dump)
        test_data = DumpTest()

        for para in dump.parameters:
            if isinstance(dump.parameters[para], np.ndarray):
                np.testing.assert_allclose(dump.parameters[para],
                                           test_data.parameters[para])
            else:
                self.assertEqual(dump.parameters[para],
                                 test_data.parameters[para])

if __name__ == '__main__':
    unittest.main(verbosity=2)
