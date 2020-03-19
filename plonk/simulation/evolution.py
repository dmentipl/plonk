"""Time series data for global or sink quantites.

This module contains a function for loading global quantities and sink
particle time series data. These files track averaged quantities that
are more frequently output than snapshot files.
"""

from __future__ import annotations

from pathlib import Path
from typing import List, Tuple, Union

import numpy as np
import pandas as pd
from pandas import DataFrame


class _Evolution:
    """Smoothed particle hydrodynamics simulation time series object."""

    def __init__(
        self,
        filenames: Union[str, Path, Tuple[str], Tuple[Path], List[str], List[Path]],
    ):
        if isinstance(filenames, (str, Path)):
            _filenames = [filenames]
        elif isinstance(filenames, (list, tuple)):
            _filenames = list(filenames)
        else:
            raise ValueError('filenames is not a known type')

        _file_paths = list()
        for filename in _filenames:
            path = Path(filename)
            _file_paths.append(path.resolve())
        self.file_paths = tuple(_file_paths)

        _check_file_consistency(self.file_paths)

        columns = _get_columns(self.file_paths[0])
        self.dataframe = self._get_data(columns)

    def _get_data(self, columns) -> DataFrame:

        times = list()
        for filename in self.file_paths:
            times.append(np.loadtxt(filename, usecols=0))

        _skiprows = [0]
        if len(times) > 1:
            for t1, t2 in zip(times, times[1:]):
                _skiprows.append(np.where(t2 < t1[-1])[0][-1] + 2)

        df = pd.concat(
            (
                pd.read_csv(
                    f,
                    names=columns,
                    skiprows=skiprows,
                    skipinitialspace=True,
                    delim_whitespace=True,
                    comment='#',
                )
                for f, skiprows in zip(self.file_paths, _skiprows)
            )
        )

        df.reset_index(inplace=True, drop=True)
        return df


def _get_columns(filename: Path) -> Tuple[str, ...]:

    with open(filename) as f:
        column_line = f.readline().strip('\n')

    _column_line = [item.strip('] ')[2:].strip(' ') for item in column_line.split('[')]

    return tuple(_column_line[1:])


def _check_file_consistency(filenames: Tuple[Path, ...]) -> None:

    columns = _get_columns(filenames[0])
    for filename in filenames:
        columns_previous = columns
        columns = _get_columns(filename)
        if columns != columns_previous:
            raise ValueError('files have different columns')


def load_ev(
    filenames: Union[str, Path, Tuple[str], Tuple[Path], List[str], List[Path]],
) -> DataFrame:
    """Load time evolution data from file(s).

    Evolution files track global quantities, such as energy, momentum,
    and density, over time. The time increments in these files is
    smaller than the snapshot file output time. These files are
    typically stored as text files.

    The data is stored as a pandas DataFrame.

    Parameters
    ----------
    filename(s)
        Collection of paths to evolution file(s) in chronological order.
        These should all contain the same columns.

    Returns
    -------
    dataframe
        A pandas DataFrame with the time evolution data.

    Examples
    --------
    Reading a single evolution file into a pandas DataFrame.

    >>> file_name = 'simulation.ev'
    >>> ev = plonk.load_ev(file_name)

    Reading a collection of evolution files into a pandas DataFrame.

    >>> file_names = ('sim01.ev', 'sim02.ev', 'sim03.ev')
    >>> ev = plonk.load_ev(file_names)
    """
    return _Evolution(filenames).dataframe
