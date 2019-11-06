"""The snap sub-package.

It contains Plonk the implementation of smoothed particle hydrodynamics
snapshot file.
"""

from .readers import load_snap
from .snap import Snap

__all__ = ['Snap', 'load_snap']
