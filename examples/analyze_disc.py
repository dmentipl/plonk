"""
analyze_disc.py

Produces similar output to the phantomanalysis module in analysis_disc.f90.

Daniel Mentiplay, 2019.
"""

import os

import matplotlib.pyplot as plt

import plonk

# ---------------------------------------------------------------------------- #

# --- Files

PREFIX = 'disc'
FILE_NUMBERS = range(1)

# --- Options

NUMBER_RADIAL_BINS = 200
RADIUS_IN = 5
RADIUS_OUT = 150

# ---------------------------------------------------------------------------- #

# --- Read dump files

DATA_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data')
DUMPFILES = [PREFIX + f'_{i:05}.h5' for i in FILE_NUMBERS]

DUMPS = list()
for dumpfile in DUMPFILES:
    file = os.path.join(DATA_PATH, dumpfile)
    print('\nReading in data from dumpfile: ' + file + '...')
    dump = plonk.Dump(file)
    DUMPS.append(dump)

# --- Perform analysis

RADIAL_AVERAGES = list()
for dump in DUMPS:
    print('\nPerforming disc analysis...\n')
    RADIAL_AVERAGES.append(
        plonk.analysis.disc(
            radius_in=RADIUS_IN,
            radius_out=RADIUS_OUT,
            number_radial_bins=NUMBER_RADIAL_BINS,
            dump=dump,
        )
    )

# --- Plot data

print('\nPlotting surface density and scale height...\n')
FIG, AX = plt.subplots(1, 2, figsize=(10, 4))

for df in RADIAL_AVERAGES:
    AX[0].plot(df['R'], df['sigma'])
    AX[1].plot(df['R'], df['H'])

AX[0].set_xlabel('radius')
AX[0].set_ylabel('surface density')
AX[1].set_xlabel('radius')
AX[1].set_ylabel('scale height')

plt.show()
