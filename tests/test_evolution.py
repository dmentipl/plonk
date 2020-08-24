"""Testing evolution files."""

from pathlib import Path

import numpy as np
import pytest

import plonk

from .data.phantom.dustseparate import mean_ev_values

TEST_FILE = Path(__file__).parent / 'data/phantom/dustseparate01.ev'


def test_read_evolution():
    """Test reading Phantom evolution files."""
    # Read from Path
    plonk.load_ev(TEST_FILE)

    # Read from str
    plonk.load_ev(str(TEST_FILE))

    # Not exists
    with pytest.raises(FileNotFoundError):
        plonk.load_ev('does_not_exist.ev')


def test_read_evolution_data():
    """Test reading data from Phantom evolution files."""
    ev = plonk.load_ev(TEST_FILE)

    assert set(ev.columns) == mean_ev_values.keys()

    for key in ev.columns:
        np.testing.assert_allclose(ev[key].mean(), mean_ev_values[key])
