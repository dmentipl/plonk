"""Reader for snapshot files."""

from pathlib import Path
from typing import Union

from ..snap import Snap
from .phantom_hdf5 import PhantomHDF5Snap

_available_formats = ('phantom',)


def load_snap(filename: Union[str, Path], filetype: str = 'phantom') -> Snap:
    """Load a snapshot from file.

    Parameters
    ----------
    filename
        Path to snapshot file.
    filetype
        The SPH file type.
    """
    if filetype.lower() in ('phantom', 'phantom_h5', 'phantom_hdf', 'phantom_hdf5'):
        filetype = 'phantom'
    if filetype not in _available_formats:
        raise ValueError('Unknown file type')
    if filetype == 'phantom':
        return PhantomHDF5Snap().generate_snap_from_file(filename)
    else:
        raise ValueError('Cannot load snapshot')
