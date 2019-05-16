"""
evolution.py

Daniel Mentiplay, 2019.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

FLOAT_TYPE = '<f8'


class Evolution:
    """
    Smoothed particle hydrodynamics simulation time evolution object.

    Evolution files track global quantities, such as energy, momentum,
    and density, over time. The time increments in these files is
    smaller than the dump file output time. These files are typically
    stored as text files.

    Parameters
    ----------
    filename(s) : list of str
        List of paths to evolution file(s) in chronological order.
        These should all contain the same columns.

    Examples
    --------
    Reading a single evolution file into an Evolution object.

    >>> file_name = 'simulation.ev'
    >>> evol = plonk.Evolution(file_name)

    Reading a list of evolution files into an Evolution object.

    >>> file_names = ['sim01.ev', 'sim02.ev', 'sim03.ev']
    >>> evol = plonk.Evolution(file_names)

    Plotting kinetic and thermal energy against time.

    >>> evol.plot('ekin', 'etherm')
    """

    def __init__(self, filenames):

        if isinstance(filenames, str):
            filenames = [filenames]

        self._file_paths = list()
        self._file_names = list()
        for filename in filenames:
            if not isinstance(filename, str) and not isinstance(filename, Path):
                raise TypeError(
                    'filenames must be a list of str or pathlib.Path'
                )
            path = Path(filename)
            self._file_paths.append(path.resolve())
            self._file_names.append(path.name)

        _check_file_consistency(filenames)
        self._columns = _get_columns(filenames[0])

        self._data = self._get_data()

    @property
    def data(self):
        """Evolution data, e.g. time, energy, momentum."""
        return self._data

    def plot(self, *args, **kwargs):
        """
        Plot evolution data as lines or markers against time.

        Parameters
        ----------
        *args : str
            The evolution data to plot against time.
        **kwargs
            Keyword arguments to pass to matplotlib.pyplot.plot.
        """

        ydata = list()
        for arg in args:
            if isinstance(arg, str):
                if arg in self._columns:
                    ydata.append(self.data[arg])
                else:
                    raise ValueError('data not available')
            else:
                raise TypeError('must specifiy column name with a string')

        xdat = self.data['time']
        for ydat in ydata:
            plt.plot(xdat, ydat, **kwargs)

    def _get_data(self):

        times = list()
        for filename in self._file_names:
            times.append(np.loadtxt(filename, usecols=0))

        final_row_index = [
            np.where(t1 < t2[0])[0][-1] for t1, t2 in zip(times, times[1:])
        ]

        dtype = [(column, FLOAT_TYPE) for column in self._columns]

        arr = [
            np.loadtxt(filename)[: final_row_index[idx]]
            for idx, filename in enumerate(self._file_names[:-1])
        ]
        arr.append(np.loadtxt(self._file_names[-1]))
        arr = np.concatenate(arr)

        return np.core.records.fromarrays(arr.T, dtype=dtype)

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        return (
            f'<plonk.Evolution: "{self._file_names}", '
            f'path="{self._file_paths}">'
        )


def _get_columns(filename):

    with open(filename) as f:
        column_line = f.readline().strip('\n')

    column_line = [
        item.strip('] ')[2:].strip(' ') for item in column_line.split('[')
    ]

    return column_line[1:]


def _check_file_consistency(filenames):

    columns = _get_columns(filenames[0])
    for filename in filenames:
        columns_previous = columns
        columns = _get_columns(filename)
        if columns != columns_previous:
            raise ValueError('files have different columns')
