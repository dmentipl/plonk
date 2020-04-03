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
        self.quantity = quantity
        self.options = kwargs
        self.ax = Visualization(snap=snaps[0]).plot(quantity=quantity, **kwargs)

        self._len = -1
        self._where = 0

    def _fn(self, idx: int):
        cbar = self.ax.images[0].colorbar
        if cbar is not None:
            cbar.remove()
        self.ax.clear()
        self.ax = Visualization(snap=self.snaps[idx]).plot(
            quantity=self.quantity, ax=self.ax, **self.options
        )

    def next(self, number: int = 1):
        """Visualize next snap."""
        idx = self._where + number
        if idx < len(self):
            self._fn(idx)
            self._where += number
        else:
            print('Too far forward. Going to last snap.')
            self._fn(len(self) - 1)
            self._where = len(self) - 1

    def prev(self, number: int = 1):
        """Visualize previous snap."""
        idx = self._where - number
        if idx > 0:
            self._fn(idx)
            self._where -= number
        else:
            print('Too far back. Going to first snap.')
            self._fn(0)
            self._where = 0

    def goto(self, idx: int):
        """Visualize particular snap by index."""
        if -len(self) < idx < len(self) - 1:
            self._fn(idx)
            self._where = np.mod(idx, len(self))
        else:
            raise ValueError('out of range')

    @property
    def index(self):
        """Get current snap index."""
        return self._where

    def __len__(self):
        """Length as number of snaps."""
        if self._len == -1:
            self._len = len(self.snaps)
        return self._len
