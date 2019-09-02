"""
Testing utils.
"""

import unittest

import numpy as np

import plonk


class TestRotateVector(unittest.TestCase):
    """Test rotating a vector."""

    def test_rotate_vector(self):

        u = np.array([1, 0, 0])
        v = np.array([0, 1, 0])
        theta = np.pi / 2

        w = plonk.core.utils.rotate_vector_arbitrary_axis(u, v, theta)
        np.testing.assert_array_almost_equal(w, np.array([0, 0, -1]))

        u = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
        w = plonk.core.utils.rotate_vector_arbitrary_axis(u, v, theta)
        np.testing.assert_array_almost_equal(
            w, np.array([[0, 0, -1], [0, 1, 0], [1, 0, 0]])
        )

        u = np.array([1, 0, 0])
        v = np.array([0, 0, 0])
        theta = np.pi / 2

        w = plonk.core.utils.rotate_vector_arbitrary_axis(u, v, theta)
        np.testing.assert_array_almost_equal(w, np.array([0, 0, 0]))


class TestNormalizeVector(unittest.TestCase):
    """Test normalizing a vector."""

    def test_normalize_vector(self):

        u = np.random.rand(99)
        v = plonk.core.utils.normalize_vector(u)

        np.testing.assert_array_almost_equal(v, v / np.linalg.norm(v))


if __name__ == '__main__':
    unittest.main(verbosity=2)
