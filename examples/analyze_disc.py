'''
analyze_disc.py

Produces similar output to the phantomanalysis module in analysis_disc.f90.

Daniel Mentiplay, 2019.
'''

import os

import matplotlib.pyplot as plt

from plonk.analysis.disc import disc_analysis
from plonk.dump import Dump

# ---------------------------------------------------------------------------- #

#--- Options for disc_analysis

NUMBER_RADIAL_BINS = 200
RADIUS_IN          = 5
RADIUS_OUT         = 150

#--- Dump file names: e.g. disc_00000.h5, ...

DATA_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data')

PREFIX = 'disc'
FILE_NUMBERS = range(1)

dumpfiles = [PREFIX + f'_{i:05}.h5' for i in FILE_NUMBERS]

# ---------------------------------------------------------------------------- #

dumps = list()

for dumpfile in dumpfiles:

#--- Read dump files

    file = os.path.join(DATA_PATH, dumpfile)

    print('\nReading in data from dumpfile: ' + file + '...')

    dump = Dump(file)

    dumps.append(dump)

# ---------------------------------------------------------------------------- #

#--- Perform analysis

radial_averages = list()

for dump in dumps:

    print('\nPerforming disc analysis...\n')
    radial_averages.append(
        disc_analysis( radius_in          = RADIUS_IN,
                       radius_out         = RADIUS_OUT,
                       number_radial_bins = NUMBER_RADIAL_BINS,
                       dump               = dump ) )

# ---------------------------------------------------------------------------- #

#--- Plot data

print('\nPlotting surface density and scale height...\n')
fig, ax = plt.subplots(1, 2, figsize=(10, 4))

for df in radial_averages:
    ax[0].plot(df['R'], df['sigma'])
    ax[1].plot(df['R'], df['H'])

ax[0].set_xlabel('radius')
ax[0].set_ylabel('surface density')
ax[1].set_xlabel('radius')
ax[1].set_ylabel('scale height')

plt.show()
