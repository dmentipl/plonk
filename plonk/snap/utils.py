"""Physical units on snapshots."""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING, Union

from .._config import load_config
from .._units import units as plonk_units

if TYPE_CHECKING:
    from .snap import SnapLike


def gravitational_constant_in_code_units(snap: SnapLike) -> float:
    """Gravitational constant in code units.

    Parameters
    ----------
    snap
        The Snap object.

    Returns
    -------
    float
        The gravitational constant in code units.
    """
    G = plonk_units.newtonian_constant_of_gravitation
    G_units = (
        snap.code_units['length'] ** 3
        / snap.code_units['mass']
        / snap.code_units['time'] ** 2
    )
    G = (G / G_units).to_base_units().magnitude
    return G


def add_aliases(snap: SnapLike, filename: Union[str, Path] = None):
    """Add array aliases to a Snap.

    Parameters
    ----------
    snap
        The Snap object.
    config : optional
        The path to a Plonk config.toml file. If None, use the default
        file.
    """
    conf = load_config(filename=filename)
    for key, val in conf['arrays']['aliases'].items():
        snap.add_alias(key, val)
