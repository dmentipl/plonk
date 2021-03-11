"""Visualize a simulation."""

from __future__ import annotations

from typing import TYPE_CHECKING, Any, Dict

import matplotlib.pyplot as plt
import numpy as np
from numpy import ndarray

from .._logging import logger
from .visualization import image, plot, vector

if TYPE_CHECKING:
    from ..simulation.simulation import Simulation

KINDS = {'image': image, 'particle': plot, 'vector': vector}


class VisualizeSimulation:
    """Visualize a simulation.

    Parameters
    ----------
    sim
        The Simulation object.
    kind
        A string in 'image', 'particle', or 'vector' which indicates the
        plot kind.
    **kwargs
        Keyword arguments to pass to the plot method.

    Examples
    --------
    Visualize a simulation by density projection images.

    >>> viz = visualize_sim(sim=sim, kind='image', quantity='density')

    Alternatively.

    >>> sim.visualize(kind='image', quantity='density')

    Go forwards and backwards through snaps.

    >>> viz.next()
    >>> viz.prev()

    Go to a particular snap, or skip ahead.

    >>> viz.goto(10)
    >>> viz.next(5)
    """

    def __init__(self, sim: Simulation, kind: str, **kwargs):
        if kind not in KINDS:
            raise ValueError(f'kind must be one of {KINDS.keys()}')

        self.sim = sim
        self.kwargs: Dict[str, Any] = kwargs
        self.ax: Any = None
        self.snaps = sim.snaps

        self._particle_ids: ndarray = None
        self._kind = kind
        self._len = -1
        self._where = 0
        self._no_new_arrays = True

        self._plotting_function(kind=self.kind, idx=0)

    @property
    def kind(self):
        """Plot kind."""
        return self._kind

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
        if isinstance(value, bool):
            self._no_new_arrays = value
        raise ValueError('no_new_arrays must be True or False')

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

        KINDS[kind](snap=snap, ax=self.ax, **self.kwargs)  # type: ignore

        new_arrays = set(snap.loaded_arrays()).symmetric_difference(loaded)
        if self._no_new_arrays:
            for array in new_arrays:
                del snap[array]

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
    Visualize a simulation by density projection images.

    >>> viz = visualize_sim(sim=sim, kind='image', quantity='density')

    Alternatively.

    >>> sim.visualize(kind='image', quantity='density')

    Go forwards and backwards through snaps.

    >>> viz.next()
    >>> viz.prev()

    Go to a particular snap, or skip ahead.

    >>> viz.goto(10)
    >>> viz.next(5)
    """
    return VisualizeSimulation(sim=sim, kind=kind, **kwargs)
