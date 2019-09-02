"""
Testing Arrays.
"""

import pathlib
import unittest

import h5py
import numpy as np

import plonk

from .stubdata.phantom_arrays import (
    dimensions,
    dtype,
    extra_quantities,
    fields,
    number,
    shape,
    structured_array_dtype,
)


class TestArrays(unittest.TestCase):
    """Test Arrays class."""

    def test_init_arrays(self):
        """Testing initialising Arrays object."""

        test_file = pathlib.Path(__file__).parent / 'stubdata/phantom_00000.h5'
        file_handle = h5py.File(test_file)
        arrays = plonk.core.particles.Arrays(
            file_handle=file_handle, arrays_label='particles'
        )

        self.assertEqual(arrays.dimensions, dimensions)
        self.assertEqual(arrays.dtype, dtype)
        self.assertEqual(arrays.number, number)
        self.assertEqual(arrays.fields, fields)
        self.assertEqual(arrays.shape, shape)

        file_handle.close()

    def test_to_structured_array(self):
        """Testing Arrays object to numpy structured array."""

        test_file = pathlib.Path(__file__).parent / 'stubdata/phantom_00000.h5'
        file_handle = h5py.File(test_file)
        arrays = plonk.core.particles.Arrays(
            file_handle=file_handle, arrays_label='particles'
        )

        self.assertEqual(arrays.to_structured_array().dtype, structured_array_dtype)

        file_handle.close()

    def test_extra_quantities(self):
        """Testing calculating extra quantities."""

        test_file = pathlib.Path(__file__).parent / 'stubdata/phantom_00000.h5'
        file_handle = h5py.File(test_file)
        arrays = plonk.core.particles.Arrays(
            file_handle=file_handle, arrays_label='particles'
        )

        for quantity in {'L', 'R', 'l', 'p', 'r', '|L|', '|l|', '|p|', '|v|'}:
            self.assertEqual(
                arrays.extra_quantity(quantity, mass=np.array([1.0]))[0].mean(),
                extra_quantities[quantity],
            )

        file_handle.close()


if __name__ == '__main__':
    unittest.main(verbosity=2)
