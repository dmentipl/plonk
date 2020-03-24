"""
Plonk
=====

**Plonk** is a Python package for analyzing and visualizing smoothed
particle hydrodynamics simulation data.

Features
--------

- Read in Phantom HDF snapshot files.
- Read in global quantity and sink evolution files.
- Encapsulate entire simulation data as Simulation object.
- Access particle and sink arrays.
- Access simulation parameters and units.
- Compute extra quantities on particles.
- Visualize data using kernel interpolation.
- Generate radial profiles.

Classes
-------

- Profile
    Represents a radial profile through the snapshot in either
    cylindrical or spherical coordinates.

- Snap
    Represents a smoothed particle hydrodynamics snapshot file,
    containing particles, sinks, and file header information.

- Simulation
    Represents an entire smoothed particle hydrodynamics simulation.
    It contains a list of Snap objects, and time series data as pandas
    dataframes.

- Visualization
    Represents a visualization of a Snap object.

Subpackages
-----------

- analysis
    Contains classes and functions for performing analysis on snapshot
    files.

- simulation
    Contains classes and functions for accessing multiple simulation
    files as a coherent object.

- snap
    Contains classes and functions for reading and accessing snapshot
    files.

- utils
    Contains utility classes and functions.

- visualize
    Contains classes and functions for visualization of snapshot files.

Documentation
-------------

See https://plonk.readthedocs.io/ for documentation. The source code is
available at https://github.com/dmentipl/plonk.
"""

from typing import Any

import pint

units: Any = pint.UnitRegistry(system='cgs')
Quantity: Any = units.Quantity

from . import analysis, simulation, snap, utils, visualize
from .analysis import Profile, load_profile
from .simulation import Simulation, load_ev, load_sim
from .snap import Snap, load_snap
from .visualize import Visualization


__all__ = (
    ['Profile', 'Simulation', 'Snap', 'Visualization']  # Classes
    + ['analysis', 'simulation', 'snap', 'utils', 'visualize']  # Packages
    + ['load_ev', 'load_profile', 'load_snap', 'load_sim']  # User functions
    + ['units', 'Quantity']
)

# Canonical version number
__version__ = '0.4.1'
