"""
Multiple plots in one figure.
"""

import matplotlib.pyplot as plt
import numpy as np

from .visualization import Visualization


class MultiPlot:
    """
    Multiple plots in one figure.

    Parameters
    ----------
    data : numpy.ndarray of dict
        A 2d numpy array of dictionaries that wrap all required
        arguments to plonk.Visualization.

    Examples
    --------
    Density rendering multiple dumps from plonk.Simuation.

    >>> viz = list()
    >>> for dump in simulation.dumps:
    >>>     viz.append(
    >>>         {
    >>>             'dump': dump,
    >>>             'render': 'density',
    >>>             'extent': [-100, 100, -100, 100],
    >>>         }
    >>>     )
    >>> viz = np.array(viz)[np.newaxis, :]
    >>> multiplot = MultiPlot(viz)
    """

    def __init__(self, data):

        self.plots = np.empty(data.shape, dtype=object)

        nrows = data.shape[0]
        ncols = data.shape[1]

        self.figure, self.axes = plt.subplots(
            nrows=nrows, ncols=ncols, squeeze=False
        )

        for idxi in range(nrows):
            for idxj in range(ncols):
                ax = self.axes[idxi, idxj]
                dat = data[idxi, idxj]
                self.plots[idxi, idxj] = Visualization(**dat, axis=ax)
