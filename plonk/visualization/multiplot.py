"""
Multiple plots in one figure.
"""

import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1 import AxesGrid

from ..core.dump import Dump
from .visualization import Visualization


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
    axes_pad : float
        Padding between axes.
    cbar_mode : str
        Colorbar mode: 'each' for one colorbar per plot or 'single' for
        one common colorbar for all plots.
    cbar_location : str
        Colorbar position: either 'top, 'bottom', 'left' or 'right'.
    cbar_pad : float
        Padding between axes and colorbar.
    cbar_label : str
        The label for the colorbar.

    Examples
    --------
    Rendering density on multiple dumps.

    >>> dumps = [dump for dump in simulation.dumps]
    >>> options = {
    ...     'render': 'density',
    ...     'extent': [-100, 100, -100, 100],
    ... }
    >>> multiplot = MultiPlot(dumps, options)

    Rendering multiple quantites on a single dump.

    >>> dump = simulation.dumps[-1]
    >>> options = list()
    >>> options.append({'render': 'density'})
    >>> options.append({'render': 'divv'})
    >>> options.append({'render': 'pressure'})
    >>> multiplot = MultiPlot(dump, options)
    """

    def __init__(
        self,
        dumps,
        options,
        shape=None,
        scale=None,
        xy_factor=None,
        axes_pad=None,
        cbar_mode=None,
        cbar_location=None,
        cbar_pad=None,
        cbar_label=None,
    ):

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

        _scale = 3.0
        _xy_factor = 1.0
        if number_options > 1:
            _cbar_mode = 'each'
            _cbar_location = 'top'
            _axes_pad = 0.20
            _cbar_pad = 0.20
        else:
            _cbar_mode = 'single'
            _cbar_location = 'right'
            _axes_pad = 0.05
            _cbar_pad = 0.05

        if scale is None:
            scale = _scale
        if xy_factor is None:
            xy_factor = _xy_factor
        if cbar_mode is None:
            cbar_mode = _cbar_mode
        if cbar_location is None:
            cbar_location = _cbar_location
        if axes_pad is None:
            axes_pad = _axes_pad
        if cbar_pad is None:
            cbar_pad = _cbar_pad

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

                vmin[idxi, idxj] = plot._options.render.render_min
                vmax[idxi, idxj] = plot._options.render.render_max
                self.plots[idxi, idxj] = plot

        if cbar_mode == 'single':

            for idxi in range(nrows):
                for idxj in range(ncols):
                    self.plots[idxi, idxj].set_render_range(
                        vmin=vmin.mean(), vmax=vmax.mean()
                    )

            self.colorbar = self.grid.cbar_axes[0].colorbar(self.plots[0, 0].image)
            if cbar_label is not None:
                self.colorbar.set_label_text(cbar_label)

    def set_render_scale(self, scale):
        """
        Set render scale.

        Parameters
        ----------
        scale : str
            A string representing the render color scale, e.g.
            'linear' or 'log'.
        """

        [viz.set_render_scale(scale) for viz in self.plots.flatten()]

    def set_colorbar_label(self, label):
        """
        Set the colorbar label.

        Parameters
        ----------
        label : str
            The label for the global colorbar.
        """

        self.colorbar.set_label_text(label)

    def set_colormap(self, cmap):
        """
        Set colormap.

        Parameters
        ----------
        cmap : str
            Colormap from Matplotlib.
        """

        [viz.set_colormap(cmap) for viz in self.plots.flatten()]

    def set_render_range(self, vmin=None, vmax=None):
        """
        Set render range for colorbar.

        Parameters
        ----------
        vmin : float
            Minimum for the render colorbar.
        vmax : float
            Maximum for the render colorbar.
        """

        if vmin is not None and vmax is not None:
            [viz.set_render_range(vmin=vmin, vmax=vmax) for viz in self.plots.flatten()]
        if vmin is not None:
            [viz.set_render_range(vmin=vmin) for viz in self.plots.flatten()]
        if vmax is not None:
            [viz.set_render_range(vmax=vmax) for viz in self.plots.flatten()]

    def set_image_size(self, extent=None, size=None):
        """
        Set image size.

        Parameters
        ----------
        extent : list or numpy.ndarray
            Extent is the image size: [xmin, xmax, ymin, ymax].
        size : float
            Extent specified by a single value:
            [-size, size, -size, size].
        """

        if extent is None and size is None:
            raise ValueError('Must set one of extent or size')

        if size is not None:
            extent = [-size, size, -size, size]

        [viz.set_image_size(extent=extent) for viz in self.plots.flatten()]
