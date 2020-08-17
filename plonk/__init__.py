"""
Plonk
=====

**Plonk** is a Python package for analyzing and visualizing smoothed
particle hydrodynamics simulation data.

Features
--------
Here are some of the features of Plonk:

- Read in Phantom HDF snapshot files.
- Read in global quantity and sink evolution files.
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

import importlib_metadata

from ._logging import logger_init as _logger_init
from ._units import units, units_defaults
from .analysis.profile import Profile, load_profile
from .simulation.evolution import load_ev
from .simulation.simulation import Simulation, load_sim
from .snap.readers import load_snap
from .snap.snap import Sinks, Snap, SnapLike, SubSnap
from .visualize.animation import animation, animation_particles, animation_profiles
from .visualize.interpolation import interpolate
from .visualize.simulation import visualize_sim
from .visualize.visualization import image, plot, vector

__version__ = importlib_metadata.version('plonk')

_logger_init(__version__)

__all__ = [
    'Profile',
    'Simulation',
    'Sinks',
    'Snap',
    'SnapLike',
    'SubSnap',
    'animation',
    'animation_particles',
    'animation_profiles',
    'image',
    'interpolate',
    'load_ev',
    'load_profile',
    'load_sim',
    'load_snap',
    'plot',
    'units',
    'units_defaults',
    'vector',
    'visualize_sim',
]
