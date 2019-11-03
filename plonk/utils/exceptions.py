"""Plonk custom exceptions.

This module contains Plonk custom exceptions.
"""


class AnalysisError(Exception):
    """Basic exception for errors raised by analysis package."""


class EvolutionError(Exception):
    """Basic exception for errors raised by Evolution class."""


class InterpolationError(Exception):
    """Basic exception for errors raised by Interpolation class."""


class SimulationError(Exception):
    """Basic exception for errors raised by Simulation class."""


class SnapError(Exception):
    """Basic exception for errors raised by Snap class."""


class VisualizationError(Exception):
    """Basic exception for errors raised by visualization."""
