"""Testing visualization.interpolation."""

import numpy as np

import plonk

from .stubdata.interpolation_arrays import (
    scalar_cross_section,
    scalar_projection,
    vector_cross_section,
    vector_projection,
)

N = 10

XX = np.ones(N)
YY = np.ones(N)
ZZ = np.ones(N)
HH = np.ones(N)
WW = np.ones(N)
MM = np.ones(N)

S_DATA = np.ones(N)
X_DATA = np.ones(N)
Y_DATA = np.ones(N)

EXTENT = (0, 1, 0, 1)
PIX = (10, 10)
ZSLICE = 0.5
HFACT = 1.0


def test_scalar_interpolation_projection():
    """Test projection interpolation."""
    im = plonk.visualize.interpolation.scalar_interpolation(
        data=S_DATA,
        x_coordinate=XX,
        y_coordinate=YY,
        z_coordinate=ZZ,
        extent=EXTENT,
        smoothing_length=HH,
        particle_mass=MM,
        hfact=HFACT,
        number_of_pixels=PIX,
    )

    np.testing.assert_allclose(im, scalar_projection, rtol=1e-5)


def test_scalar_interpolation_cross_section():
    """Test cross section interpolation."""
    im = plonk.visualize.interpolation.scalar_interpolation(
        data=S_DATA,
        x_coordinate=XX,
        y_coordinate=YY,
        z_coordinate=ZZ,
        extent=EXTENT,
        smoothing_length=HH,
        particle_mass=MM,
        hfact=HFACT,
        number_of_pixels=PIX,
        cross_section=ZSLICE,
    )

    np.testing.assert_allclose(im, scalar_cross_section, rtol=1e-5)


def test_vector_interpolation_projection():
    """Test projection interpolation."""
    vec = plonk.visualize.interpolation.vector_interpolation(
        x_data=X_DATA,
        y_data=Y_DATA,
        x_coordinate=XX,
        y_coordinate=YY,
        z_coordinate=ZZ,
        extent=EXTENT,
        smoothing_length=HH,
        particle_mass=MM,
        hfact=HFACT,
        number_of_pixels=PIX,
    )

    np.testing.assert_allclose(vec, vector_projection, rtol=1e-5)


def test_vector_interpolation_cross_section():
    """Test cross section interpolation."""
    vec = plonk.visualize.interpolation.vector_interpolation(
        x_data=X_DATA,
        y_data=Y_DATA,
        x_coordinate=XX,
        y_coordinate=YY,
        z_coordinate=ZZ,
        extent=EXTENT,
        smoothing_length=HH,
        particle_mass=MM,
        hfact=HFACT,
        number_of_pixels=PIX,
        cross_section=ZSLICE,
    )

    np.testing.assert_allclose(vec, vector_cross_section, rtol=1e-5)
