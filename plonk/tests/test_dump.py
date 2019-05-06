"""
test_dump.py
"""

import os
import unittest

import numpy as np

from plonk.dump import Dump
from .test_data.disc_00000 import DumpTest


class TestDump(unittest.TestCase):
    """Test reading Phantom HDF dump files."""

    def test_read_dump(self):
        """Testing reading Phantom HDF dump file."""

        file_dir = os.path.dirname(os.path.abspath(__file__))
        test_dump = os.path.join(file_dir, 'test_data', 'disc_00000.h5')
        dump = Dump(test_dump)
        test_data = DumpTest()

        for para in dump.parameters:
            if isinstance(dump.parameters[para], np.ndarray):
                np.testing.assert_allclose(
                    dump.parameters[para], test_data.parameters[para]
                )
            else:
                self.assertEqual(
                    dump.parameters[para], test_data.parameters[para]
                )


if __name__ == '__main__':
    unittest.main(verbosity=2)
