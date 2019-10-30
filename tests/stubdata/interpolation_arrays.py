"""Data to test interpolation functions."""

import pathlib

import numpy as np

DIR = pathlib.Path(__file__).parent

scalar_projection = np.loadtxt(DIR / 'scalar_projection.csv')
scalar_cross_section = np.loadtxt(DIR / 'scalar_cross_section.csv')
vector_projection = np.array(
    (
        np.loadtxt(DIR / 'vector_x_projection.csv'),
        np.loadtxt(DIR / 'vector_y_projection.csv'),
    )
)
vector_cross_section = np.array(
    (
        np.loadtxt(DIR / 'vector_x_cross_section.csv'),
        np.loadtxt(DIR / 'vector_y_cross_section.csv'),
    )
)
