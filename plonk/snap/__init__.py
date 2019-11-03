"""The snap sub-package.

It contains Plonk the implementation of smoothed particle hydrodynamics
snapshot file.
"""

from .snap import Snap
from .readers import load_snap

__all__ = ['Snap', 'load_snap']
