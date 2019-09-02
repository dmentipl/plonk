"""
Testing Simulation.
"""

import pathlib
import unittest

import plonk


class TestReadPhantomSimulation(unittest.TestCase):
    """Test reading Phantom simultion data."""

    def test_init_simulation(self):
        """Testing initialising simulation."""

        dir_path = pathlib.Path(__file__).parent / 'stubdata'
        plonk.Simulation(prefix='phantom', directory=dir_path)
