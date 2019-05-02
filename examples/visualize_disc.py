'''
visualize_disc.py

Produces similar rendering to Splash.

Daniel Mentiplay, 2019.
'''

import os

import matplotlib.pyplot as plt

from plonk.dump import Dump
from plonk.visualization.image import plot

# ---------------------------------------------------------------------------- #

#--- Options for plotting

RENDER              = 'rho'           # Render density
RENDER_FRACTION_MAX = 0.05            # Colour bar scaled to % of max.
IMAGE_RANGE         = 150             # Figure width in code units

#--- Dump file names: e.g. disc_00000.h5, ...

DATA_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data')

PREFIX = 'disc'
FILE_NUMBERS = range(1)

dumpfiles = [PREFIX + f'_{i:05}.h5' for i in FILE_NUMBERS]

# ---------------------------------------------------------------------------- #

#--- Read dump files

dumps = list()
for dumpfile in dumpfiles:
    file = os.path.join(DATA_PATH, dumpfile)
    print('\nReading in data from dumpfile: ' + file + '...')
    dumps.append(Dump(file))

# ---------------------------------------------------------------------------- #

#--- Plot images

print(f'\nPlotting dump files...\n')
for dump in dumps:
    plt.figure()
    plot(dump, render=RENDER, render_fraction_max=RENDER_FRACTION_MAX,
         image_range=IMAGE_RANGE)

plt.show()
