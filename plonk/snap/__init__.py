"""Snapshot data access and manipulation.

It contains Plonk the implementation of smoothed particle hydrodynamics
snapshot file.
"""

from pathlib import Path
from typing import Union

from ..utils.strings import is_documented_by
from .snap import Sinks, Snap, SnapLike, SubSnap


@is_documented_by(Snap.load_snap)
def load_snap(
    filename: Union[str, Path],
    data_source: str = 'phantom',
    config: Union[str, Path] = None,
):
    return Snap().load_snap(filename=filename, data_source=data_source, config=config)


__all__ = [
    'Sinks',
    'Snap',
    'SnapLike',
    'SubSnap',
    'load_snap',
]
