"""Testing evolution files."""

from pathlib import Path

import numpy as np
import pytest

import plonk

from .data.phantom import dustseparate

DIR = Path(__file__).parent / 'data/phantom'


def test_read_evolution():
    """Test reading Phantom evolution files."""
    # Read from Path
    filename = DIR / dustseparate.ev_file
    plonk.load_ev(filename)

    # Read from str
    plonk.load_ev(str(filename))

    # Not exists
    with pytest.raises(FileNotFoundError):
        plonk.load_ev('does_not_exist.ev')


def test_read_evolution_data():
    """Test reading data from Phantom evolution files."""
    filename = DIR / dustseparate.ev_file
    ev = plonk.load_ev(filename)

    assert set(ev.columns) == dustseparate.mean_ev_values.keys()

    for key in ev.columns:
        np.testing.assert_allclose(ev[key].mean(), dustseparate.mean_ev_values[key])
