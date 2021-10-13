"""
Plonk
=====

**Plonk** is a Python package for analyzing and visualizing smoothed
particle hydrodynamics simulation data.

Features
--------
Here are some of the features of Plonk:

- Read in Phantom HDF snapshot files.
- Read in global quantity and sink time series data files.
- Encapsulate entire simulation data as Simulation object.
- Access particle and sink arrays.
- Access simulation parameters and units.
- Compute extra quantities on particles.
- Visualize data using kernel interpolation.
- Generate profiles.

Classes
-------
Profile
    Represents a profile through the snapshot in Cartesian, cylindrical
    or spherical coordinates.
Sinks
    Represents sinks in a snapshot file.
Snap
    Represents a smoothed particle hydrodynamics snapshot file,
    containing particles, sinks, and file header information.
SubSnap
    Represents a subset of particles in a snapshot file.
Simulation
    Represents an entire smoothed particle hydrodynamics simulation.
    It contains a list of Snap objects, and time series data as pandas
    dataframes.

Subpackages
-----------
analysis
    Perform analysis on SPH snapshot data.
simulation
    Access multiple simulation files as a coherent object.
snap
    Read and accessing snapshot files.
utils
    Utility functions.
visualize
    Visualize of snapshot files.

Documentation
-------------
See https://plonk.readthedocs.io/ for documentation. The source code is
available at https://github.com/dmentipl/plonk.
"""

from ._config import read_config, write_config
from ._logging import logger_init as _logger_init
from ._units import Quantity, add_units, array_units, units
from .analysis.profile import Profile, load_profile
from .simulation.simulation import Simulation, load_sim, load_simulation
from .simulation.time_series import load_ev, load_time_series
from .snap import load_snap
from .snap.snap import Sinks, Snap, SnapLike, SubSnap
from .visualize.animation import animate
from .visualize.interpolation import interpolate
from .visualize.simulation import visualize_sim
from .visualize.visualization import image, plot, vector

__version__ = '0.7.4'

_logger_init(__version__)

add_units()

__all__ = [
    'Profile',
    'Quantity',
    'Simulation',
    'Sinks',
    'Snap',
    'SnapLike',
    'SubSnap',
    'add_units',
    'animate',
    'array_units',
    'image',
    'interpolate',
    'load_ev',
    'load_profile',
    'load_sim',
    'load_simulation',
    'load_snap',
    'load_time_series',
    'plot',
    'read_config',
    'units',
    'vector',
    'visualize_sim',
    'write_config',
]
