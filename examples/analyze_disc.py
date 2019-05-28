"""
analyze_disc.py

Produces similar output to the phantomanalysis module in
analysis_disc.f90.

Daniel Mentiplay, 2019.
"""

import pathlib

import matplotlib.pyplot as plt
import plonk

# ---------------------------------------------------------------------------- #

# --- Files

DIRECTORY = pathlib.Path('~/runs/prefix/dir')
PREFIX = 'prefix'

# --- Options

NUMBER_RADIAL_BINS = 100
RADIUS_IN = 5
RADIUS_OUT = 250

# ---------------------------------------------------------------------------- #

# --- Read dump files

SIMULATION = plonk.Simulation(prefix=PREFIX, directory=DIRECTORY)

# --- Perform analysis

print('\nPerforming disc analysis...\n')
RADIAL_AVERAGES = list()
for dump in SIMULATION.dumps:
    print(f'{dump.file_name}')
    RADIAL_AVERAGES.append(
        plonk.analysis.disc(
            dump=dump,
            radius_in=RADIUS_IN,
            radius_out=RADIUS_OUT,
            number_radial_bins=NUMBER_RADIAL_BINS,
        )
    )

# --- Plot data

print('\nPlotting surface density and scale height...')
FIG, AX = plt.subplots(1, 2, figsize=(10, 4))

for df in RADIAL_AVERAGES:
    AX[0].plot(df['R'], df['sigma'])
    AX[1].plot(df['R'], df['H'])

AX[0].set_xlabel('radius')
AX[0].set_ylabel('surface density')
AX[1].set_xlabel('radius')
AX[1].set_ylabel('scale height')

plt.show()
