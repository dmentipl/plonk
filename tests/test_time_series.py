"""Testing time series data files."""

from pathlib import Path

import numpy as np
import pytest

import plonk

from .data.phantom import dustseparate

DIR = Path(__file__).parent / 'data/phantom'


def test_read_time_series():
    """Test reading Phantom time series data files."""
    # Read from Path
    filename = DIR / dustseparate.ts_file
    plonk.load_time_series(filename)

    # Read from str
    plonk.load_time_series(str(filename))

    # Not exists
    with pytest.raises(FileNotFoundError):
        plonk.load_time_series('does_not_exist.ev')


def test_read_time_series_data():
    """Test reading data from Phantom time series files."""
    filename = DIR / dustseparate.ts_file
    ts = plonk.load_time_series(filename)

    assert set(ts.columns) == dustseparate.mean_ts_values.keys()

    for key in ts.columns:
        np.testing.assert_allclose(ts[key].mean(), dustseparate.mean_ts_values[key])
