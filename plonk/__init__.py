"""
Plonk
=====

**Plonk** is a Python package for analyzing and visualizing smoothed
particle hydrodynamics simulation data.

Features
--------

- Read in Phantom HDF snapshot files.
- Read in global evolution files.
- Encapsulate entire simulation data as Simulation object.
- Access particle and sink arrays.
- Access simulation parameters and units.
- Compute extra quantities on particles.
- Visualize data using kernel interpolation.

Classes
-------

- Snap
    Represents a smoothed particle hydrodynamics snapshot file,
    containing particles, sinks, and file header information.

- Evolution
    Represents globally averaged quantities as time series data.

- Simulation
    Represents an entire smoothed particle hydrodynamics simulation.
    It contains instances of Snap and Evolution objects.

- Visualization
    Represents a visualization of a Snap object.

Subpackages
-----------

- analysis
    Contains classes and functions for performing analysis on snapshot
    files.

- visualization
    Contains classes and functions for visualization of snapshot files.

Documentation
-------------

See https://plonk.readthedocs.io/ for documentation. The source code is
available at https://github.com/dmentipl/plonk.
"""

from . import analysis, simulation, snap, utils, visualization
from .snap import Snap, load_snap
from .simulation import Evolution, Simulation, load_ev, load_sim
from .visualization import Visualization

__all__ = [
    'Evolution',
    'Simulation',
    'Snap',
    'analysis',
    'snap',
    'load_snap',
    'load_ev',
    'load_sim',
    'simulation',
    'utils',
    'visualization',
    'Visualization',
]

# Canonical version number
__version__ = '0.1.0'
