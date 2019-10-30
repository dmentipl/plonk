"""
Testing Visualization.
"""

import pathlib
import unittest

import plonk

I_GAS = 1
I_DUST = 7

TEST_FILE = pathlib.Path(__file__).parent / 'stubdata/phantom_00000.h5'


class TestInitializeVisualization(unittest.TestCase):
    """Test initialization of Visualization object."""

    def test_initialization(self):

        dump = plonk.Dump(TEST_FILE)

        density = dump.density
        position = dump.particles.arrays['xyz'][:]
        smoothing_length = dump.particles.arrays['h'][:]
        particle_mass = dump.mass
        x_position = position[:, 0]
        y_position = position[:, 1]

        extent = (-150, 150, -150, 150)
        scalar_options = {'norm': 'lin', 'cmap': 'gist_heat'}
        interpolation_options = {
            'number_of_pixels': (512, 512),
            'cross_section': None,
            'density_weighted': False,
        }

        plonk.Visualization(
            scalar_data=density,
            x_coordinate=x_position,
            y_coordinate=y_position,
            extent=extent,
            particle_mass=particle_mass,
            smoothing_length=smoothing_length,
            scalar_options=scalar_options,
            interpolation_options=interpolation_options,
        )


@unittest.skip('Not yet converted to new interpolation method')
class TestInitializeMultiPlot(unittest.TestCase):
    """Test initialization of MultiPlot object."""

    def test_initialization(self):

        dump = plonk.Dump(TEST_FILE)
        options = list()
        options.append({'render': 'density'})
        options.append({'render': 'divv'})

        plonk.visualization.MultiPlot(dump, options)

        dumps = [dump, dump]
        options = {'render': 'density', 'extent': [-100, 100, -100, 100]}

        plonk.visualization.MultiPlot(dumps, options)


if __name__ == '__main__':
    unittest.main(verbosity=2)
