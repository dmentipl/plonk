"""Testing evolution files."""

import pathlib

import numpy as np
import pytest

import plonk

from .stubdata.phantom_evolution import mean_values

TEST_FILE = pathlib.Path(__file__).parent / 'stubdata/phantom01.ev'


def test_read_evolution():
    """Test reading Phantom evolution files."""
    # Read from pathlib.Path
    plonk.load_ev(TEST_FILE)

    # Read from str
    plonk.load_ev(str(TEST_FILE))

    # Not exists
    with pytest.raises(FileNotFoundError):
        plonk.load_ev('does_not_exist.ev')


def test_read_evolution_data():
    """Test reading data from Phantom evolution files."""
    ev = plonk.load_ev(TEST_FILE)

    assert set(ev.columns) == mean_values.keys()

    for key in ev.columns:
        np.testing.assert_allclose(ev[key].mean(), mean_values[key])
