"""
Testing reading and writing dumps.
"""

import pathlib
import unittest

import numpy as np
import plonk


class TestPhantomDump(unittest.TestCase):
    """Test reading Phantom HDF dump files."""

    def test_read_dump_parameters(self):
        """Testing reading Phantom HDF dump file parameters."""

        test_file = pathlib.Path(__file__).parent / 'test_data' / 'disc_00000.h5'

        test_header = plonk.Dump(test_file).header

        for para in test_header:
            if isinstance(test_header[para], np.ndarray):
                np.testing.assert_allclose(test_header[para], reference_header[para])
            else:
                self.assertEqual(test_header[para], reference_header[para])


if __name__ == '__main__':
    unittest.main(verbosity=2)
