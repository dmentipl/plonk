"""
visualize_disc.py

Produces similar rendering to Splash.

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

RENDER = 'rho'  # Render density
RENDER_FRACTION_MAX = 0.05  # Colour bar scaled to % of max.
IMAGE_RANGE = 150  # Figure width in code units

# ---------------------------------------------------------------------------- #

# --- Read dump files

DATA_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data')
DUMPFILES = [PREFIX + f'_{i:05}.h5' for i in FILE_NUMBERS]

DUMPS = list()
for dumpfile in DUMPFILES:
    file = os.path.join(DATA_PATH, dumpfile)
    print('\nReading in data from dumpfile: ' + file + '...')
    DUMPS.append(plonk.Dump(file))

# --- Plot images

print(f'\nPlotting dump files...\n')
for dump in DUMPS:
    plt.figure()
    plonk.plot(
        dump,
        render=RENDER,
        render_fraction_max=RENDER_FRACTION_MAX,
        image_range=IMAGE_RANGE,
    )

plt.show()
