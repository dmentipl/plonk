"""Visualize multiple snaps together."""

from typing import List, Union

import numpy as np
from numpy import ndarray

from ..snap import SnapLike
from .visualization import Visualization


class MultiVisualization:
    """Visualize multiple snaps.

    Attributes
    ----------
    snaps
        A list of snaps.
    options
        A dictionary of arguments passed to Visualization plot method.
    visualization
        The Visualization object.
    ax
        The matplotlib Axes object.

    Parameters
    ----------
    snaps
        A list of Snap objects.
    quantity
        The quantity to visualize.
    **kwargs
        Keyword arguments to pass to Visualization plot method.
    """

    def __init__(self, snaps: List[SnapLike], quantity: Union[str, ndarray], **kwargs):
        self.snaps = snaps
        self.options = kwargs
        self.quantity = quantity
        viz = Visualization(snap=snaps[0])
        self.visualization = viz.plot(quantity=quantity, **kwargs)
        self.ax = self.visualization.ax

        self._len = -1
        self._where = 0

    def _fn(self, idx: int):
        self.ax.clear()
        cbar = self.visualization.objects['colorbar']
        if cbar is not None:
            cbar.remove()
        viz = Visualization(snap=self.snaps[idx])
        viz.plot(quantity=self.quantity, ax=self.ax, **self.options)
        self.visualization = viz
        return viz

    def next(self, number: int = 1):
        """Visualize next snap."""
        idx = self._where + number
        if idx < len(self):
            self.visualization = self._fn(idx)
            self._where += number
        else:
            print('Too far forward. Going to last snap.')
            self.visualization = self._fn(len(self) - 1)
            self._where = len(self) - 1

    def prev(self, number: int = 1):
        """Visualize previous snap."""
        idx = self._where - number
        if idx > 0:
            self.visualization = self._fn(idx)
            self._where -= number
        else:
            print('Too far back. Going to first snap.')
            self.visualization = self._fn(0)
            self._where = 0

    def goto(self, idx: int):
        """Visualize particular snap by index."""
        if -len(self) < idx < len(self) - 1:
            self.visualization = self._fn(idx)
            self._where = np.mod(idx, len(self))
        else:
            raise ValueError('out of range')

    @property
    def index(self):
        """Current snap index."""
        return self._where

    def __len__(self):
        """Length as number of snaps."""
        if self._len == -1:
            self._len = len(self.snaps)
        return self._len
