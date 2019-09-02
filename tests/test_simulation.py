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
        sim = plonk.Simulation(prefix='phantom', directory=dir_path)

        dumps = sim.dumps
        self.assertEqual(len(dumps), 1)
        self.assertEqual(dumps[0].file_name, 'phantom_00000.h5')

        ev = sim.evolution
        self.assertEqual(ev.file_names[0], 'phantom01.ev')


if __name__ == '__main__':
    unittest.main(verbosity=2)
