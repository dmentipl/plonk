"""
Multiple plots in one figure.
"""

import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1 import AxesGrid

from .visualization import Visualization

SCALE = 3


class MultiPlot:
    """
    Multiple plots of the same quantity over different dump files.

    Parameters
    ----------
    dumps : numpy.ndarray of Dump
        A 1d or 2d numpy array of Dump objects arranged as required for
        the multiplot.
    **kwargs
        The required arguments to plonk.Visualization.

    Examples
    --------
    Density rendering multiple dumps from plonk.Simuation.

    >>> dumps = np.array([dump for dump in simulation.dumps])
    >>> options = {
            'render': 'density',
            'extent': [-100, 100, -100, 100],
        }
    >>> multiplot = MultiPlot(dumps, **options)
    """

    def __init__(self, dumps, **kwargs):

        if not isinstance(dumps, np.ndarray):
            raise TypeError('dumps argument must be numpy.ndarray')
        if dumps.ndim != 2:
            if dumps.ndim == 1:
                dumps = dumps[np.newaxis, :]

        self.plots = np.empty(dumps.shape, dtype=object)

        nrows = dumps.shape[0]
        ncols = dumps.shape[1]

        self.figure = plt.figure(figsize=(ncols * SCALE, nrows * SCALE))

        self.grid = AxesGrid(
            self.figure,
            111,
            nrows_ncols=(nrows, ncols),
            axes_pad=0.05,
            cbar_mode='single',
            cbar_location='right',
            cbar_pad=0.05,
        )

        vmin = np.zeros_like(dumps)
        vmax = np.zeros_like(dumps)

        for idxi in range(nrows):
            for idxj in range(ncols):

                ax = self.grid.axes_row[idxi][idxj]
                dump = dumps[idxi, idxj]
                plot = Visualization(
                    dump=dump, **kwargs, colorbar=False, axis=ax
                )
                vmin[idxi, idxj] = plot._vmin
                vmax[idxi, idxj] = plot._vmax
                self.plots[idxi, idxj] = plot

                if not idxi == nrows - 1:
                    ax.set_xlabel('')
                    ax.tick_params(labelbottom=False)
                if not idxj == 0:
                    ax.set_ylabel('')
                    ax.tick_params(labelleft=False)

        for idxi in range(nrows):
            for idxj in range(ncols):
                self.plots[idxi, idxj].set_render_range(
                    vmin=vmin.mean(), vmax=vmax.mean()
                )

        self.colorbar = self.grid.cbar_axes[0].colorbar(self.plots[0, 0].image)
