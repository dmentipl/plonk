"""Radial binning of disc quantities.

Produces similar output to the Phantom disc analysis module.

Daniel Mentiplay, 2019.
"""

import pathlib

import matplotlib.pyplot as plt
import plonk

# ---------------------------------------------------------------------------- #

# --- Files

DIRECTORY = pathlib.Path('~/runs/directory')
PREFIX = 'prefix'

# --- Options

NUMBER_RADIAL_BINS = 100
RADIUS_IN = 5
RADIUS_OUT = 250

# ---------------------------------------------------------------------------- #

# --- Read dump files

simulation = plonk.Simulation(prefix=PREFIX, directory=DIRECTORY)

# --- Perform analysis

print('\nPerforming disc analysis...\n')
radial_averages = list()
for dump in simulation.dumps:
    print(f'{dump.file_name}')
    radial_averages.append(
        plonk.analysis.disc(
            dump=dump,
            radius_in=RADIUS_IN,
            radius_out=RADIUS_OUT,
            number_radial_bins=NUMBER_RADIAL_BINS,
        )
    )

# --- Plot data

print('\nPlotting surface density and scale height...')
figure, axis = plt.subplots(1, 2, figsize=(10, 4))

for df in radial_averages:
    axis[0].plot(df['R'], df['sigma'])
    axis[1].plot(df['R'], df['H'])

axis[0].set_xlabel('radius')
axis[0].set_ylabel('surface density')
axis[1].set_xlabel('radius')
axis[1].set_ylabel('scale height')

plt.show()
