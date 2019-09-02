"""
Testing Dump.
"""

import pathlib
import unittest

import numpy as np
import plonk

from .stubdata.phantom_dump import header, array_keys, mean_array_values


class TestReadPhantomDump(unittest.TestCase):
    """Test reading Phantom HDF5 dump files."""

    def test_read_dump(self):
        """Testing reading Phantom HDF5 dumps."""

        # Read from pathlib.Path
        test_file = pathlib.Path(__file__).parent / 'stubdata/phantom_00000.h5'
        dump = plonk.Dump(test_file)
        dump._close_file()

        # Read from str
        test_file = str(pathlib.Path(__file__).parent) + '/stubdata/phantom_00000.h5'
        dump = plonk.Dump(test_file)
        dump._close_file()

        # Not exists
        test_file = 'does_not_exist.h5'
        self.assertRaises(FileNotFoundError, plonk.Dump, test_file)
        dump._close_file()

    def test_read_header(self):
        """Testing reading Phantom HDF5 dump header."""

        test_file = pathlib.Path(__file__).parent / 'stubdata/phantom_00000.h5'
        dump = plonk.Dump(test_file)

        for para in dump.header:
            if isinstance(dump.header[para], np.ndarray):
                np.testing.assert_allclose(dump.header[para], header[para])
            else:
                self.assertEqual(dump.header[para], header[para])

        dump._close_file()

    def test_read_particle_arrays(self):
        """Testing reading Phantom HDF5 dump particle arrays."""

        test_file = pathlib.Path(__file__).parent / 'stubdata/phantom_00000.h5'
        dump = plonk.Dump(test_file)

        self.assertEqual(set(dump.particles.arrays), array_keys)

        for key, val in dump.particles.arrays.items():
            np.testing.assert_allclose(val[()].mean(), mean_array_values[key])

        dump._close_file()


if __name__ == '__main__':
    unittest.main(verbosity=2)
