"""
Multiple plots in one figure.
"""

import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1 import AxesGrid

from ..core.dump import Dump
from .visualization import Visualization

SCALE = 3.0


class MultiPlot:
    """
    Multiple plots in the same figure.

    Multiple plots of the same quantity over different dump files, or
    multiple plots of different quantities over the same dump file.

    Parameters
    ----------
    dumps : Dump or list of Dump
        A single Dump or a list of Dump objects in order for the
        multiplot.
    options : dict or list of dict
        A single dict or a list of dicts with keyword arguments to
        plonk.Visualization.

    Other parameters
    ----------------
    shape : tuple, optional
        Specify the number of rows and columns to tile the multiplot.
        The default is one row with multiple columns.
    scale : float, optional (default 3.0)
        Figure size scaling factor.
    xy_factor : float, optional (default 1.0)
        Figure width/height ratio. Greater than unity makes width wider,
        and height shorter. Useful for tweaking plots.

    Examples
    --------
    Rendering density on multiple dumps.

    >>> dumps = [dump for dump in simulation.dumps]
    >>> options = {
            'render': 'density',
            'extent': [-100, 100, -100, 100],
        }
    >>> multiplot = MultiPlot(dumps, options)

    Rendering multiple quantites on a single dump.

    >>> dump = simulation.dumps[-1]
    >>> options = list()
    >>> options.append({'render': 'density'})
    >>> options.append({'render': 'divv'})
    >>> options.append({'render': 'pressure'})
    >>> multiplot = MultiPlot(dumps, options)
    """

    def __init__(self, dumps, options, shape=None, scale=None, xy_factor=None):

        if scale is None:
            scale = SCALE

        if xy_factor is None:
            xy_factor = 1.0

        if isinstance(dumps, Dump):
            number_dumps = 1
            dumps = [dumps]
        elif isinstance(dumps, (list, tuple)):
            number_dumps = len(dumps)
        if isinstance(options, dict):
            number_options = 1
            options = [options]
        elif isinstance(options, (list, tuple)):
            number_options = len(options)

        if number_dumps > 1 and number_options > 1:
            raise ValueError(
                'Must have multiple dumps and one set of options OR '
                'one dump and multiple sets of options.'
            )

        if number_options > 1:
            cbar_mode = 'each'
            cbar_location = 'top'
            axes_pad = 0.20
            cbar_pad = 0.20
        else:
            cbar_mode = 'single'
            cbar_location = 'right'
            axes_pad = 0.05
            cbar_pad = 0.05

        if number_dumps == 1 and number_options > 1:
            dumps = number_options * dumps
        if number_options == 1 and number_dumps > 1:
            options = number_dumps * options

        if shape is None:
            shape = (1, max(number_dumps, number_options))
        else:
            if shape[0] * shape[1] != len(dumps):
                raise ValueError('shape is incompatible with inputs')

        dumps = np.array(dumps).reshape(shape)
        options = np.array(options).reshape(shape)

        self.plots = np.empty(dumps.shape, dtype=object)

        nrows = dumps.shape[0]
        ncols = dumps.shape[1]

        x_scale = scale * xy_factor
        y_scale = scale / xy_factor
        self.figure = plt.figure(figsize=(ncols * x_scale, nrows * y_scale))

        self.grid = AxesGrid(
            self.figure,
            111,
            nrows_ncols=(nrows, ncols),
            axes_pad=axes_pad,
            cbar_mode=cbar_mode,
            cbar_location=cbar_location,
            cbar_pad=cbar_pad,
            cbar_size='5%',
            label_mode='L',
        )

        vmin = np.zeros_like(dumps)
        vmax = np.zeros_like(dumps)

        for idxi in range(nrows):
            for idxj in range(ncols):

                axis = self.grid.axes_row[idxi][idxj]
                dump = dumps[idxi, idxj]
                option = options[idxi, idxj]

                cbar_axis = None
                colorbar = False
                if cbar_mode == 'each':
                    cbar_axis = self.grid.cbar_axes[idxi + idxj]
                    colorbar = True

                plot = Visualization(
                    dump=dump,
                    **option,
                    colorbar=colorbar,
                    cbar_axis=cbar_axis,
                    axis=axis
                )

                vmin[idxi, idxj] = plot._vmin
                vmax[idxi, idxj] = plot._vmax
                self.plots[idxi, idxj] = plot

        if cbar_mode == 'single':

            for idxi in range(nrows):
                for idxj in range(ncols):
                    self.plots[idxi, idxj].set_render_range(
                        vmin=vmin.mean(), vmax=vmax.mean()
                    )

            self.colorbar = self.grid.cbar_axes[0].colorbar(
                self.plots[0, 0].image
            )
