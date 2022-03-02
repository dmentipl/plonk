"""Time series data for global or sink quantites.

This module contains a function for loading global quantities and sink
particle time series data typical of Phantom simulations. These files
track averaged quantities that are more frequently output than snapshot
files.
"""

from __future__ import annotations

import warnings
from pathlib import Path
from typing import TYPE_CHECKING, List, Tuple, Union

from pandas import DataFrame

from .._logging import logger
from ._phantom_ev import load_data_from_file as load_data_from_file_phantom
from ._phantom_ev import time_series_units as time_series_units_phantom

if TYPE_CHECKING:
    from .simulation import Simulation


def load_ev(
    filenames: Union[str, Path, Tuple[str], Tuple[Path], List[str], List[Path]],
    data_source: str = 'Phantom',
    config: Union[str, Path] = None,
) -> DataFrame:
    """Load time series data from file(s).

    Time series files track global quantities, such as energy, momentum,
    and density, over time. The time increments in these files is
    smaller than the snapshot file output time. These files are
    typically stored as text files.

    The data is stored as a pandas DataFrame.

    Parameters
    ----------
    filename(s)
        Collection of paths to time series file(s) in chronological
        order. These should all contain the same columns.
    data_source : optional
        The code used to produce the data. Default is 'Phantom'.
    config : optional
        The path to a Plonk config.toml file.

    Returns
    -------
    dataframe
        A pandas DataFrame with the time series data.

    Examples
    --------
    Reading a single time series file into a pandas DataFrame.

    >>> file_name = 'simulation.ev'
    >>> ts = plonk.load_time_series(file_name)

    Reading a collection of time series files into a pandas DataFrame.

    >>> file_names = ('sim01.ev', 'sim02.ev', 'sim03.ev')
    >>> ts = plonk.load_time_series(file_names)
    """
    msg = (
        'load_ev is deprecated and will be removed in v0.7.4, '
        'please use load_time_series instead'
    )
    logger.warning(msg)
    warnings.warn(msg, DeprecationWarning)
    return load_time_series(filenames=filenames, data_source=data_source, config=config)


def load_time_series(
    filenames: Union[str, Path, Tuple[str], Tuple[Path], List[str], List[Path]],
    data_source: str = 'Phantom',
    config: Union[str, Path] = None,
) -> DataFrame:
    """Load time series data from file(s).

    Time series files track global quantities, such as energy, momentum,
    and density, over time. The time increments in these files is
    smaller than the snapshot file output time. These files are
    typically stored as text files.

    The data is stored as a pandas DataFrame.

    Parameters
    ----------
    filename(s)
        Collection of paths to time series file(s) in chronological
        order. These should all contain the same columns.
    data_source : optional
        The code used to produce the data. Default is 'Phantom'.
    config : optional
        The path to a Plonk config.toml file.

    Returns
    -------
    dataframe
        A pandas DataFrame with the time series data.

    Examples
    --------
    Reading a single time series file into a pandas DataFrame.

    >>> file_name = 'simulation.ev'
    >>> ts = plonk.load_time_series(file_name)

    Reading a collection of time series files into a pandas DataFrame.

    >>> file_names = ('sim01.ev', 'sim02.ev', 'sim03.ev')
    >>> ts = plonk.load_time_series(file_names)
    """
    if data_source == 'Phantom':
        return load_data_from_file_phantom(filenames=filenames, config=config)
    raise ValueError('Cannot determine code used to produce time series data')


def time_series_units(
    sim: Simulation, data_source: str = 'Phantom', config: Union[str, Path] = None
):
    """Get units of time series data from Simulation object.

    Parameters
    ----------
    sim
        The Simulation object.
    data_source : optional
        The code used to produce the data. Default is 'Phantom'.
    config : optional
        The path to a Plonk config.toml file.

    Returns
    -------
    Dict
    """
    if data_source == 'Phantom':
        return time_series_units_phantom(sim=sim, config=config)
    raise ValueError('Cannot determine code used to produce time series data')
