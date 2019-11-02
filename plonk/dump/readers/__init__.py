"""Reader for dump files."""

from pathlib import Path
from typing import Union

from ..dump import Dump
from .phantom_hdf5 import PhantomHDF5Dump

_available_formats = ('phantom',)


def load_dump(filename: Union[str, Path], dump_type: str = 'phantom') -> Dump:
    """Load a dump from file.

    Parameters
    ----------
    filename
        Path to dump file.
    dump_type
        The type of SPH dump.
    """
    if dump_type.lower() in ('phantom', 'phantom_h5', 'phantom_hdf', 'phantom_hdf5'):
        dump_type = 'phantom'
    if dump_type not in _available_formats:
        raise ValueError('Unknown dump type')
    if dump_type == 'phantom':
        return PhantomHDF5Dump().generate_dump_from_file(filename)
    else:
        raise ValueError('Cannot load dump')
