"""
This is the visualization module.

It contains the Plonk implementation for visualizing smoothed particle
hydrodynamics simulations using Splash for interpolation. This is
accessed by the plot function.

Examples
--------
Rendering density on a Dump.

>>> viz = plonk.Visualization(dump, render='density')

Go forwards and backwards through visualizations in a Simulation.

>>> viz_iter = VisualizationIterator(
...     dumps=simulation.dumps,
...     render=render
... )
>>> viz_iter.next()
>>> viz_iter.previous()

Density rendering multiple dumps from plonk.Simuation.

>>> dumps = np.array([dump for dump in simulation.dumps])
>>> options = {
...     'render': 'density',
...     'extent': [-100, 100, -100, 100],
... }
>>> multiplot = MultiPlot(dumps, **options)

Notes
-----
For the Splash source code see https://github.com/danieljprice/splash.
The user guide is at http://users.monash.edu.au/~dprice/splash. When
using any software derived from Splash for academic purposes you should
cite: Price, 2007, Publ. Astron. Soc. Aust., 24, 159-173.
"""

from .multiplot import MultiPlot
from .visualization import Visualization

__all__ = ['MultiPlot', 'Visualization']
