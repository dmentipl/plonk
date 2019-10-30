"""Produce a density rendering of an SPH dump.

Daniel Mentiplay, 2019.
"""

import pathlib

import matplotlib.pyplot as plt

import plonk

# ---------------------------------------------------------------------------- #

# --- File

DATA_FILE = pathlib.Path('~/runs/directory/file_12345')

# --- Read dump file

dump = plonk.Dump(DATA_FILE)

position = dump.particles.arrays['xyz'][:]
smoothing_length = dump.particles.arrays['h'][:]
density = dump.density
particle_mass = dump.mass

x_position = position[:, 0]
y_position = position[:, 1]

# --- Options

scalar_options = {'norm': 'linear', 'cmap': 'gist_heat'}
interpolation_options = {'number_of_pixels': (512, 512)}
extent = (-100.0, 100.0, -100.0, 100.0)

# --- Plot image

viz_density = plonk.Visualization(
    scalar_data=dump.density,
    x_coordinate=x_position,
    y_coordinate=y_position,
    extent=extent,
    particle_mass=particle_mass,
    smoothing_length=smoothing_length,
    scalar_options=scalar_options,
    interpolation_options=interpolation_options,
)

plt.show()
