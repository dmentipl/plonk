"""Physical units on snapshots."""

from __future__ import annotations

from typing import TYPE_CHECKING

from numpy import ndarray

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
    G_units = snap.units['length'] ** 3 / snap.units['mass'] / snap.units['time'] ** 2
    G = (G / G_units).to_base_units().magnitude
    return G


def get_array_in_code_units(snap: SnapLike, name: str) -> ndarray:
    """Get array in code units.

    Parameters
    ----------
    snap
        The Snap or SubSnap.
    name
        The array name.

    Returns
    -------
    ndarray
        The array on the particles in code units.
    """
    return (snap[name] / snap.get_array_code_unit(name)).magnitude
