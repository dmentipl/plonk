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

import logging
import platform
from typing import Any

import pint

# Canonical version number
__version__ = '0.5.1'

# Units
units: Any = pint.UnitRegistry(system='cgs')
Quantity: Any = units.Quantity

# Logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

console_handler = logging.StreamHandler()
file_handler = logging.FileHandler('.plonk.log')

console_handler.setLevel(logging.INFO)
file_handler.setLevel(logging.DEBUG)

console_format = logging.Formatter('%(levelname)s - %(message)s')
file_format = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
console_handler.setFormatter(console_format)
file_handler.setFormatter(file_format)

logger.addHandler(console_handler)
logger.addHandler(file_handler)


def get_os_info():
    """Get the operating system version for logging."""
    system = platform.system()
    if system == 'Darwin':
        system = 'macOS'
    release = platform.release()
    return f'{system} version: {release}'


logger.debug(f'Plonk v{__version__} on Python {platform.python_version()}')
logger.debug(f'{get_os_info()}, {platform.machine()}')

from . import analysis, simulation, snap, utils, visualize
from .analysis import Profile, load_profile
from .simulation import Simulation, load_ev, load_sim
from .snap import Snap, load_snap

__all__ = (
    ['Profile', 'Simulation', 'Snap']  # Classes
    + ['analysis', 'simulation', 'snap', 'utils', 'visualize']  # Packages
    + ['load_ev', 'load_profile', 'load_snap', 'load_sim']  # User functions
    + ['units', 'Quantity']
)
