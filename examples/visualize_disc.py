"""
visualize_disc.py

Produces similar rendering to Splash.

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

RENDER = 'density'
EXTENT = [-150, 150, -150, 150]

# ---------------------------------------------------------------------------- #

# --- Read dump files

simulation = plonk.Simulation(prefix=PREFIX, directory=DIRECTORY)

# --- Plot images

print(f'\nPlotting dump files...\n')
for dump in simulation.dumps:
    print(f'{dump.file_name}')
    plt.figure()
    plonk.Visualization(dump, render=RENDER, extent=EXTENT)

plt.show()
