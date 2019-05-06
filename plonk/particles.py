"""
particles.py

Daniel Mentiplay, 2019.
"""

import numpy as np

I_GAS = 1
I_DUST = 7


def density_from_smoothing_length(smoothing_length, particle_mass, hfact=1.2):
    """Calculate density from particle mass and smoothing length."""

    return particle_mass * (hfact / np.abs(smoothing_length)) ** 3
