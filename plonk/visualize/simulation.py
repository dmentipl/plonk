"""Visualize a simulation."""

from __future__ import annotations

from typing import TYPE_CHECKING, Any, Dict

import matplotlib.pyplot as plt
import numpy as np
from numpy import ndarray

from .._logging import logger

if TYPE_CHECKING:
    from ..simulation.simulation import Simulation


class VisualizeSimulation:
    """Visualize a simulation.

    Parameters
    ----------
    sim
        The Simulation object.
    """

    def __init__(self, sim: Simulation):
        self.sim = sim
        self.options: Dict[str, Any] = {}
        self.ax: Any = None
        self.snaps = sim.snaps

        self._particle_ids: ndarray = None
        self._kind = ''
        self._len = -1
        self._where = 0
        self._no_new_arrays = True

    @property
    def kind(self):
        """Plot kind."""
        return self._kind

    @kind.setter
    def kind(self, value: str):
        if value not in ('image', 'particle', 'vector'):
            raise ValueError('kind must be one of ("image", "particle", "vector")')
        self._kind = value

    @property
    def particle_ids(self):
        """Particle ids of snaps."""
        return self._particle_ids

    @particle_ids.setter
    def particle_ids(self, value: ndarray):
        subsnaps = [snap[value] for snap in self.snaps]
        self.snaps = subsnaps
        self._particle_ids = value

    @property
    def no_new_arrays(self):
        """Do not allow new arrays on Snap objects."""
        return self._no_new_arrays

    @no_new_arrays.setter
    def no_new_arrays(self, value):
        for snap in self.snaps:
            snap._no_new_arrays = value
        self._no_new_arrays = value

    def _plotting_function(self, kind: str, idx: int):
        if self.ax is None:
            _, self.ax = plt.subplots()
        try:
            im = self.ax.images[0]
            cbar = im.colorbar
            cbar.remove()
        except (IndexError, AttributeError):
            pass
        self.ax.clear()

        snap = self.snaps[idx]
        loaded = set(snap.loaded_arrays())

        if kind == 'image':
            snap.image(ax=self.ax, **self.options)
        elif kind == 'particle':
            snap.plot(ax=self.ax, **self.options)
        elif kind == 'vector':
            snap.vector(ax=self.ax, **self.options)

        new_arrays = set(snap.loaded_arrays()).symmetric_difference(loaded)
        if self._no_new_arrays:
            for array in new_arrays:
                del snap[array]

    def start(self):
        """Start visualization."""
        self._plotting_function(kind=self.kind, idx=0)

    def next(self, number: int = 1):
        """Visualize next snap."""
        idx = self._where + number
        if idx < len(self):
            self._plotting_function(kind=self.kind, idx=idx)
            self._where += number
        else:
            logger.info('Too far forward. Going to last snap.')
            self._plotting_function(kind=self.kind, idx=len(self) - 1)
            self._where = len(self) - 1

    def prev(self, number: int = 1):
        """Visualize previous snap."""
        idx = self._where - number
        if idx > 0:
            self._plotting_function(kind=self.kind, idx=idx)
            self._where -= number
        else:
            logger.info('Too far back. Going to first snap.')
            self._plotting_function(kind=self.kind, idx=0)
            self._where = 0

    def goto(self, idx: int):
        """Visualize particular snap by index."""
        if -len(self) < idx < len(self) - 1:
            self._plotting_function(kind=self.kind, idx=idx)
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


def visualize_sim(sim: Simulation, kind: str, **kwargs) -> VisualizeSimulation:
    """Visualize a simulation.

    Parameters
    ----------
    sim
        The Simulation object.
    kind
        The kind of plot: 'particle', 'image', 'vector'.
    **kwargs
        Keyword arguments to pass to plotting methods such as
        plonk.image, plonk.plot, and plonk.vector.

    Returns
    -------
    VisualizeSimulation

    Examples
    --------
    Initialize object passing in plotting parameters.

    >>> viz = visualize_sim(
    ...     sim=sim,
    ...     kind='image',
    ...     quantity='density',
    ... )

    Go forwards and backwards through snaps.

    >>> viz.next()
    >>> viz.prev()

    Go to a particular snap, or skip ahead.

    >>> viz.goto(10)
    >>> viz.next(5)
    """
    viz = VisualizeSimulation(sim=sim)
    viz.kind = kind
    viz.options = kwargs
    viz.start()

    return viz
