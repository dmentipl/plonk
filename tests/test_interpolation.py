"""
Testing visualization.interpolation.
"""

import unittest

import numpy as np

import plonk

from .stubdata.interpolation_arrays import (
    scalar,
    scalar_cross_section,
    scalar_perspective,
    vector,
    vector_cross_section,
)

N = 10

XYZ = np.ones((N, 3))
H = np.ones((N,))
W = np.ones_like(H)
M = np.ones_like(H)
DX = DY = (0, 1)
PIX = (10, 10)

SDAT = np.ones_like(H)
VDAT = np.ones_like(XYZ)


class TestScalarInterpolation(unittest.TestCase):
    """Test interpolation of scalar quantities."""

    def test_scalar_interpolation(self):

        im = plonk.visualization.interpolation.scalar_interpolation(
            XYZ, H, W, SDAT, M, DX, DY, number_pixels=PIX
        )

        np.testing.assert_allclose(im, scalar, rtol=1e-5)

    def test_scalar_cross_section(self):

        cross_section = True
        slice_position = 0.5

        im = plonk.visualization.interpolation.scalar_interpolation(
            XYZ,
            H,
            W,
            SDAT,
            M,
            DX,
            DY,
            number_pixels=PIX,
            cross_section=cross_section,
            slice_position=slice_position,
        )

        np.testing.assert_allclose(im, scalar_cross_section, rtol=1e-5)

    def test_scalar_perspective(self):

        perspective = True
        observer_distance = 10.0

        im = plonk.visualization.interpolation.scalar_interpolation(
            XYZ,
            H,
            W,
            SDAT,
            M,
            DX,
            DY,
            number_pixels=PIX,
            perspective=perspective,
            observer_distance=observer_distance,
        )

        np.testing.assert_allclose(im, scalar_perspective, rtol=1e-5)


@unittest.skip('Test fails due to bugs in Splash')
class TestVectorInterpolation(unittest.TestCase):
    """Test interpolation of vector quantities."""

    def test_vector_interpolation(self):

        vec = plonk.visualization.interpolation.vector_interpolation(
            XYZ, H, W, VDAT, DX, DY, number_pixels=PIX
        )

        np.testing.assert_allclose(vec, vector, rtol=1e-5)

    def test_vector_cross_section(self):

        cross_section = True
        slice_position = 0.5

        vec = plonk.visualization.interpolation.vector_interpolation(
            XYZ,
            H,
            W,
            VDAT,
            DX,
            DY,
            number_pixels=PIX,
            cross_section=cross_section,
            slice_position=slice_position,
        )

        np.testing.assert_allclose(vec, vector_cross_section, rtol=1e-5)


if __name__ == '__main__':
    unittest.main(verbosity=2)
