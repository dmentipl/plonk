"""Testing visualization.interpolation."""

import numpy as np

from plonk.visualize.interpolation import scalar_interpolation, vector_interpolation

from .stubdata.interpolation_arrays import (
    scalar_slice,
    scalar_projection,
    vector_slice,
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
    im = scalar_interpolation(
        quantity=S_DATA,
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


def test_scalar_interpolation_slice():
    """Test cross section interpolation."""
    im = scalar_interpolation(
        quantity=S_DATA,
        x_coordinate=XX,
        y_coordinate=YY,
        z_coordinate=ZZ,
        extent=EXTENT,
        smoothing_length=HH,
        particle_mass=MM,
        hfact=HFACT,
        number_of_pixels=PIX,
        z_slice=ZSLICE,
    )

    np.testing.assert_allclose(im, scalar_slice, rtol=1e-5)


def test_vector_interpolation_projection():
    """Test projection interpolation."""
    vec = vector_interpolation(
        quantity_x=X_DATA,
        quantity_y=Y_DATA,
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


def test_vector_interpolation_slice():
    """Test cross section interpolation."""
    vec = vector_interpolation(
        quantity_x=X_DATA,
        quantity_y=Y_DATA,
        x_coordinate=XX,
        y_coordinate=YY,
        z_coordinate=ZZ,
        extent=EXTENT,
        smoothing_length=HH,
        particle_mass=MM,
        hfact=HFACT,
        number_of_pixels=PIX,
        z_slice=ZSLICE,
    )

    np.testing.assert_allclose(vec, vector_slice, rtol=1e-5)
