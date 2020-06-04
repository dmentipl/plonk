"""Reader for snapshot files."""

from pathlib import Path
from typing import Union

from ... import logger
from ..snap import Snap
from .phantom import generate_snap_from_file as read_phantom

_data_sources = ('Phantom',)


def load_snap(filename: Union[str, Path], data_source: str = 'Phantom') -> Snap:
    """Load a snapshot from file.

    Parameters
    ----------
    filename
        Path to snapshot file.
    data_source : optional
        The SPH software that produced the data. Default is 'Phantom'.

    Returns
    -------
    Snap
        The Snap object.
    """
    if data_source not in _data_sources:
        raise ValueError(
            f'Unknown data source. Available data sources:\n{_data_sources}'
        )

    if data_source == 'Phantom':
        try:
            return read_phantom(filename)
        except OSError:
            # Catch errors raised by h5py due to file corruption
            logger.error(f'Cannot read file: {filename}')
    raise ValueError('Cannot load snapshot')
