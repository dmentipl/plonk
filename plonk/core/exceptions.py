"""
Plonk custom exceptions

This module contains Plonk custom exceptions.
"""


class ParticlesError(Exception):
    """Basic exception for errors raised by Particles class."""


class DumpError(Exception):
    """Basic exception for errors raised by Dump class."""


class EvolutionError(Exception):
    """Basic exception for errors raised by Evolution class."""


class SimulationError(Exception):
    """Basic exception for errors raised by Simulation class."""


class UnitsError(Exception):
    """Basic exception for errors raised by Units class."""


class VisualizationError(Exception):
    """Basic exception for errors raised by visualization."""


class InterpolationError(Exception):
    """Basic exception for errors raised by Interpolation class."""


class MultiPlotError(Exception):
    """Basic exception for errors raised by MultiPlot class."""


class AnalysisError(Exception):
    """Basic exception for errors raised by analysis package."""
