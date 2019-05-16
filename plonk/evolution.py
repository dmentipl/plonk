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
    Smoothed particle hydrodynamics evolution (time series) file object.

    Parameters
    ----------
    filename : str
        Path to evolution file.

    Examples
    --------
    Reading an evolution file into an Evolution object.

    >>> file_name = 'simulation01.ev'
    >>> evol = plonk.Evolution(file_name)

    Plotting kinetic and thermal energy against time.

    >>> evol.plot('ekin', 'etherm')
    """

    def __init__(self, filename):

        if not isinstance(filename, str) and not isinstance(filename, Path):
            raise TypeError('filename must be str or pathlib.Path')

        path = Path(filename)
        self._file_path = path.resolve()
        self._file_name = path.name

        self._columns = self._get_columns()

    @property
    def data(self):
        """Evolution data, e.g. time, energy, momentum."""
        return self._get_data()

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

    def _get_columns(self):

        with open(self._file_name) as f:
            column_line = f.readline().strip('\n')

        column_line = [
            item.strip('] ')[2:].strip(' ') for item in column_line.split('[')
        ]

        return column_line[1:]

    def _get_data(self):
        dtype = [(column, FLOAT_TYPE) for column in self._columns]
        return np.loadtxt(self._file_name, dtype=dtype)

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        return (
            f'<plonk.Evolution: "{self._file_name}", '
            f'path="{self._file_path}">'
        )
