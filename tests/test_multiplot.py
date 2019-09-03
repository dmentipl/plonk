"""
Testing MultiPlot.
"""

import contextlib
import io
import pathlib
import unittest

import plonk

TEST_FILE = pathlib.Path(__file__).parent / 'stubdata/phantom_00000.h5'

TEXT_TRAP = io.StringIO()


class TestInitializeMultiPlot(unittest.TestCase):
    """Test initialization of MultiPlot object."""

    def test_initialization(self):

        dump = plonk.Dump(TEST_FILE)
        options = list()
        options.append({'render': 'density'})
        options.append({'render': 'divv'})

        with contextlib.redirect_stdout(TEXT_TRAP):
            plonk.visualization.MultiPlot(dump, options)

        dumps = [dump, dump]
        options = {'render': 'density', 'extent': [-100, 100, -100, 100]}

        with contextlib.redirect_stdout(TEXT_TRAP):
            plonk.visualization.MultiPlot(dumps, options)


if __name__ == '__main__':
    unittest.main(verbosity=2)
