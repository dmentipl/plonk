"""Testing SPH neighbours and kernels."""

from pathlib import Path

import numpy as np
import pytest

import plonk
from plonk import analysis

TEST_FILE = Path(__file__).parent / 'stubdata/phantom_00000.h5'
NEIGHBOURS_FILE = Path(__file__).parent / 'stubdata/neighbours.csv'

with open(NEIGHBOURS_FILE) as fp:
    NEIGHBOURS = list()
    lines = fp.readlines()
    for line in lines:
        if line[0] == '#':
            continue
        NEIGHBOURS.append([int(item) for item in line.strip().split(',')])


def test_get_neighbours():
    """Testing getting particle neighbours."""
    snap = plonk.load_snap(TEST_FILE)

    with pytest.raises(ValueError):
        snap.set_kernel('kernel_does_not_exist')

    snap.set_kernel('cubic')

    subsnap = snap['gas']
    indices = np.arange(len(subsnap))
    neigh = subsnap.neighbours(indices)
    for idxi, idxj in enumerate(indices):
        assert set(neigh[idxi]) == set(NEIGHBOURS[idxj])

    snap.close_file()


def test_sph_derivative():
    """Test SPH derivative."""
    snap = plonk.load_snap(TEST_FILE)

    analysis.sph.derivative(
        snap=snap, derivative='div', quantity='velocity', kernel='cubic',
    )

    snap.close_file()
