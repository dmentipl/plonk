"""
This is the analysis module.

It contains Plonk implementations of typical smoothed particle
hydrodynamics post-simulation analysis tasks.

Examples
--------
Analyzing a disc by radially binning from the central object.

>>> radial_averages = plonk.analysis.disc(
...     dump=dump,
...     radius_in=radius_in,
...     radius_out=radius_out,
...     number_radial_bins=number_radial_bins,
... )
"""

from .disc import disc_analysis as disc

__all__ = ['disc']
