"""Data to test interpolation functions."""

import pathlib

import numpy as np

DIR = pathlib.Path(__file__).parent

scalar_projection = np.loadtxt(DIR / 'scalar_projection.csv')
scalar_slice = np.loadtxt(DIR / 'scalar_slice.csv')
vector_projection = np.array(
    (
        np.loadtxt(DIR / 'vector_x_projection.csv'),
        np.loadtxt(DIR / 'vector_y_projection.csv'),
    )
)
vector_slice = np.array(
    (np.loadtxt(DIR / 'vector_x_slice.csv'), np.loadtxt(DIR / 'vector_y_slice.csv'),)
)
