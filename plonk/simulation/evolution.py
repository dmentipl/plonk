"""Evolution class for global quantites.

This module contains the Evolution class for tracking global quantities
and sink particle time series data. These files track averaged
quantities that are more frequently output than snapshot files.
"""

from __future__ import annotations

from pathlib import Path
from typing import List, Tuple, Union

import numpy as np
import pandas as pd
from numpy import ndarray
from pandas import DataFrame, Series


class Evolution:
    """Smoothed particle hydrodynamics simulation time series object.

    Evolution files track global quantities, such as energy, momentum,
    and density, over time. The time increments in these files is
    smaller than the snapshot file output time. These files are
    typically stored as text files.

    The data is stored as a pandas DataFrame.

    Examples
    --------
    Reading a single evolution file into an Evolution object.

    >>> file_name = 'simulation.ev'
    >>> ev = plonk.load_ev(file_name)

    Reading a collection of evolution files into an Evolution object.

    >>> file_names = ('sim01.ev', 'sim02.ev', 'sim03.ev')
    >>> ev = plonk.load_ev(file_names)

    See what columns and what data are available.

    >>> ev.columns
    >>> ev.data

    Accessing part of the data as a pandas DataFrame or Series.

    >>> time = ev['time']
    >>> first_100 = ev[:100]

    Plot columns using pandas plotting methods.

    >>> ev.plot('time', ['x', 'y'])
    """

    def __init__(self):

        self.file_paths: Tuple[Path, ...]
        self.file_names: Tuple[str, ...]

    def load_from_file(
        self,
        filenames: Union[str, Path, Tuple[str], Tuple[Path], List[str], List[Path]],
    ) -> Evolution:
        """Load from file(s).

        Parameters
        ----------
        filename(s)
            Collection of paths to evolution file(s) in chronological order.
            These should all contain the same columns.
        """
        if isinstance(filenames, (str, Path)):
            _filenames = [filenames]
        elif isinstance(filenames, (list, tuple)):
            _filenames = list(filenames)
        else:
            raise ValueError('filenames is not a known type')

        _file_paths = list()
        _file_names = list()
        for filename in _filenames:
            path = Path(filename)
            _file_paths.append(path.resolve())
            _file_names.append(path.name)

        self.file_paths = tuple(_file_paths)
        self.file_names = tuple(_file_names)
        _check_file_consistency(self.file_paths)

        self._columns = _get_columns(self.file_paths[0])
        self._data = self._get_data()

        return self

    @property
    def columns(self) -> Tuple[str, ...]:
        """List of available time evolution data."""
        return self._columns

    @property
    def data(self) -> DataFrame:
        """Time evolution data as a pandas DataFrame."""
        return self._data

    def plot(self, *args, **kwargs):
        """Plot using pandas."""
        return self.data.plot(*args, **kwargs)

    def __getitem__(
        self, inp: Union[str, List[str], ndarray, int, slice]
    ) -> Union[DataFrame, Series]:
        """Return a pandas object."""
        if isinstance(inp, str):
            return self.data[inp]
        elif isinstance(inp, (ndarray, int, slice)):
            return self.data.iloc[inp]
        elif isinstance(inp, list):
            if isinstance(inp[0], str):
                return self.data[inp]
            elif isinstance(inp[0], int):
                return self.data.iloc[inp]
        raise ValueError('Cannot determine item to return')

    def __setitem__(self, name: str, item: Union[ndarray, Series]):
        """Set a column with ndarray or pandas Series."""
        if not isinstance(item, (ndarray, Series)):
            raise ValueError('"item" must be NumPy ndarray or pandas Series')
        if item.shape[0] != len(self):
            raise ValueError('Length of array does not match length')
        if name in self.columns:
            raise ValueError('Attempting overwrite existing column')
        else:
            self.data[name] = item

    def _get_data(self) -> DataFrame:

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
                    names=self._columns,
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

    def __repr__(self):
        """Dunder repr method."""
        return self.__str__()

    def __str__(self):
        """Dunder str method."""
        return f'<plonk.Evolution: "{self.file_names}">'

    def __len__(self):
        """Dunder len method."""
        return len(self.data)


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
) -> Evolution:
    """Load Evolution from file(s).

    Parameters
    ----------
    filename(s)
        Collection of paths to evolution file(s) in chronological order.
        These should all contain the same columns.
    """
    return Evolution().load_from_file(filenames)
