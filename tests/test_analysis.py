"""
Testing Analysis.
"""

import pathlib
import unittest
import warnings

import pandas as pd

import plonk

I_GAS = 1
I_DUST = 7

TEST_FILE = pathlib.Path(__file__).parent / 'stubdata/phantom_00000.h5'
CSV_FILE = pathlib.Path(__file__).parent / 'stubdata/phantom_00000_analysis.csv'


class TestDiscAnalysis(unittest.TestCase):
    """Test disc analysis."""

    def test_disc_analysis(self):

        dump = plonk.Dump(TEST_FILE)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            av = plonk.analysis.disc(
                dump, radius_in=1, radius_out=150, number_radial_bins=10
            )
        pd.testing.assert_frame_equal(av, pd.read_csv(CSV_FILE, index_col=0))


if __name__ == '__main__':
    unittest.main(verbosity=2)
