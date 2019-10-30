"""Testing visualization.interpolation."""

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

XX = np.ones(N)
YY = np.ones(N)
ZZ = np.ones(N)
HH = np.ones(N)
WW = np.ones(N)
MM = np.ones(N)
EXTENT = (0, 1, 0, 1)
PIX = (10, 10)

S_DATA = np.ones(N)
X_DATA = np.ones(N)
Y_DATA = np.ones(N)

ZSLICE = 0.5


class TestScalarInterpolation(unittest.TestCase):
    """Test interpolation of scalar quantities."""

    def test_scalar_interpolation(self):

        im = plonk.visualization.interpolation.scalar_interpolation(
            data=S_DATA,
            x_position=XX,
            y_position=YY,
            z_position=ZZ,
            extent=EXTENT,
            smoothing_length=HH,
            particle_mass=MM,
            number_of_pixels=PIX,
        )

        np.testing.assert_allclose(im, scalar, rtol=1e-5)

    def test_scalar_cross_section(self):

        im = plonk.visualization.interpolation.scalar_interpolation(
            data=S_DATA,
            x_position=XX,
            y_position=YY,
            z_position=ZZ,
            extent=EXTENT,
            smoothing_length=HH,
            particle_mass=MM,
            number_of_pixels=PIX,
            cross_section=ZSLICE,
        )

        np.testing.assert_allclose(im, scalar_cross_section, rtol=1e-5)


class TestVectorInterpolation(unittest.TestCase):
    """Test interpolation of vector quantities."""

    def test_vector_interpolation(self):

        vec = plonk.visualization.interpolation.vector_interpolation(
            x_data=X_DATA,
            y_data=Y_DATA,
            x_position=XX,
            y_position=YY,
            z_position=ZZ,
            extent=EXTENT,
            smoothing_length=HH,
            particle_mass=MM,
            number_of_pixels=PIX,
        )

        np.testing.assert_allclose(vec, vector, rtol=1e-5)

    def test_vector_cross_section(self):

        vec = plonk.visualization.interpolation.vector_interpolation(
            x_data=X_DATA,
            y_data=Y_DATA,
            x_position=XX,
            y_position=YY,
            z_position=ZZ,
            extent=EXTENT,
            smoothing_length=HH,
            particle_mass=MM,
            number_of_pixels=PIX,
            cross_section=ZSLICE,
        )

        np.testing.assert_allclose(vec, vector_cross_section, rtol=1e-5)


if __name__ == '__main__':
    unittest.main(verbosity=2)
