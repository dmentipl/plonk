"""Testing Simulation."""

import pathlib

import plonk


def test_init_simulation():
    """Testing initialising simulation."""
    dir_path = pathlib.Path(__file__).parent / 'stubdata'
    sim = plonk.load_sim(prefix='phantom', directory=dir_path)

    snaps = sim.snaps
    assert len(snaps) == 1
    assert len(snaps[0]) == 2000

    assert sim.paths['global_quantities'][0].name == 'phantom01.ev'
