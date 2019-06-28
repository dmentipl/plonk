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
- Visualize data using interpolation provided by Splash.

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

- Units
    Represents the units of a Dump or Simulation.

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

from . import analysis, visualization
from .core.constants import constants
from .core.dump import Dump
from .core.evolution import Evolution
from .core.simulation import Simulation
from .core.units import Units
from .visualization import Visualization

__all__ = [
    'Dump',
    'Evolution',
    'Simulation',
    'Units',
    'Visualization',
    'analysis',
    'constants',
    'visualization',
]

# Canonical version number
__version__ = '0.1.0'
