"""
Plonk
=====

**Plonk** is a Python package for analyzing and visualizing smoothed
particle hydrodynamics simulation data.

Features
--------

- Read in Phantom HDF dump files.
- Read in global evolution files.
- Encapsulate entire simulation data as Simulation object.
- Access particle and sink arrays.
- Access simulation parameters and units.
- Compute extra quantities on particles.
- Visualize data using kernel interpolation.

Classes
-------

- Dump
    Represents a smoothed particle hydrodynamics dump file,
    containing particles, sinks, and file header information.

- Evolution
    Represents globally averaged quantities as time series data.

- Simulation
    Represents an entire smoothed particle hydrodynamics simulation.
    It contains instances of Dump and Evolution objects.

- Visualization
    Represents a visualization of a Dump object.

Subpackages
-----------

- analysis
    Contains classes and functions for performing analysis on dump
    files.

- visualization
    Contains classes and functions for visualization of dump files.

Documentation
-------------

See https://plonk.readthedocs.io/ for documentation. The source code is
available at https://github.com/dmentipl/plonk.
"""

from . import analysis, dump, simulation, utils, visualization
from .dump import load_dump
from .simulation import load_ev, load_sim

__all__ = [
    'analysis',
    'dump',
    'load_dump',
    'load_ev',
    'load_sim',
    'simulation',
    'utils',
    'visualization',
]

# Canonical version number
__version__ = '0.1.0'
