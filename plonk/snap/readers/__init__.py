"""Reader for snapshot files."""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING, Union

from ..._logging import logger
from .phantom import generate_snap_from_file as read_phantom

if TYPE_CHECKING:
    from ..snap import Snap

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
        except FileNotFoundError as e:
            logger.error(f'File not found: {filename}')
            raise e
        except OSError as e:
            # Catch errors raised by h5py due to file corruption
            logger.error(f'File likely corrupted: {filename}')
            raise e
    raise RuntimeError('Cannot load snap')
