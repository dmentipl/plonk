"""
Plonk
=====

**Plonk** is a Python package for analyzing and visualizing smoothed particle
hydrodynamics simulation data.

Features
--------

  - Read in Phantom HDF dump files.
  - Access particle arrays as pandas DataFrame objects.
  - Access simulation parameters and units.
  - Compute extra quantities on particles.
  - Visualize data using interpolation provided by Splash.

See https://plonk.readthedocs.io/ for documentation. Code is available at
https://github.com/dmentipl/plonk.
"""

from . import analysis
from .dump import Dump
from .simulation import Simulation
from .visualization.image import plot

__all__ = ['Dump', 'Simulation', 'analysis', 'plot']
