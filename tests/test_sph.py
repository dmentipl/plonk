"""Testing SPH neighbours and kernels."""

from pathlib import Path

import numpy as np
import pytest

import plonk
from plonk import analysis

from .data.phantom import dustseparate

DIR = Path(__file__).parent / 'data/phantom'


def test_get_neighbours():
    """Testing getting particle neighbours."""
    filename = DIR / dustseparate.filename
    snap = plonk.load_snap(filename)

    filename = DIR / dustseparate.neighbours_file
    neighbours = _neighbours(filename)

    with pytest.raises(ValueError):
        snap.set_kernel('kernel_does_not_exist')

    snap.set_kernel('cubic')

    subsnap = snap['gas']
    indices = np.arange(len(subsnap))
    neigh = subsnap.neighbours(indices)
    for idxi, idxj in enumerate(indices):
        assert set(neigh[idxi]) == set(neighbours[idxj])

    snap.close_file()


def test_sph_derivative():
    """Test SPH derivative."""
    filename = DIR / dustseparate.filename
    snap = plonk.load_snap(filename)

    analysis.sph.derivative(
        snap=snap, derivative='div', quantity='velocity', kernel='cubic',
    )

    snap.close_file()


def _neighbours(neighbours_file):
    with open(neighbours_file) as fp:
        neighbours = list()
        lines = fp.readlines()
        for line in lines:
            if line[0] == '#':
                continue
            neighbours.append([int(item) for item in line.strip().split(',')])
    return neighbours
