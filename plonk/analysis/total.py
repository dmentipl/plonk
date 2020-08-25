"""Calculate global (total) quantities on the particles."""

from __future__ import annotations

from typing import TYPE_CHECKING, List

import numpy as np

from .._units import Quantity
from ..utils.math import norm
from . import particles

if TYPE_CHECKING:
    from ..snap.snap import SnapLike, SubSnap


def accreted_mass(snap: SnapLike) -> Quantity:
    """Calculate the accreted mass.

    Parameters
    ----------
    snap
        The Snap object.

    Returns
    -------
    Quantity
        The accreted mass.
    """
    h: Quantity = snap['smoothing_length']
    _mass: Quantity = snap['mass'][~(h > 0)]

    return _mass.sum()


def angular_momentum(
    snap: SnapLike, sinks: bool = True, origin: Quantity = None
) -> Quantity:
    """Calculate the total angular momentum.

    Parameters
    ----------
    snap
        The Snap object.
    sinks : optional
        Include sink particles. Default is True.
    origin : optional
        The origin around which to compute the angular momentum as a
        Quantity like (x, y, z) * au. Default is (0, 0, 0).

    Returns
    -------
    Quantity
        The total angular momentum like (lx, ly, lz).
    """
    angmom = particles.angular_momentum(
        snap=snap, origin=origin, ignore_accreted=True
    ).sum(axis=0)
    if sinks and snap.num_sinks > 0:
        angmom_sinks = particles.angular_momentum(
            snap=snap.sinks, origin=origin, ignore_accreted=False
        )
        if angmom_sinks.ndim == 1:
            angmom = angmom + angmom_sinks
        elif angmom_sinks.ndim == 2:
            angmom = angmom + angmom_sinks.sum(axis=0)
    return angmom


def center_of_mass(snap: SnapLike, sinks: bool = True) -> Quantity:
    """Calculate the center of mass.

    Parameters
    ----------
    snap
        The Snap object.
    sinks : optional
        Include sink particles. Default is True.

    Returns
    -------
    Quantity
        The center of mass as a vector (cx, cy, cz).
    """
    h: Quantity = snap['smoothing_length']
    _mass: Quantity = snap['mass'][h > 0]
    pos: Quantity = snap['position'][h > 0]

    mass_pos = (_mass[:, np.newaxis] * pos).sum(axis=0)
    mass_tot = _mass.sum()

    if sinks and snap.num_sinks > 0:
        mass_sinks = snap.sinks['mass']
        pos_sinks = snap.sinks['position']
        if pos_sinks.ndim == 1:
            mass_pos = mass_pos + mass_sinks * pos_sinks
            mass_tot = mass_tot + mass_sinks
        elif pos_sinks.ndim == 2:
            mass_pos = mass_pos + (mass_sinks[:, np.newaxis] * pos_sinks).sum(axis=0)
            mass_tot = mass_tot + mass_sinks.sum()

    return mass_pos / mass_tot


def dust_mass(snap: SnapLike, squeeze: bool = False) -> Quantity:
    """Calculate the total dust mass per species.

    Parameters
    ----------
    snap
        The Snap object.
    squeeze
        If True return all subtypes in a single array. Default is
        False.

    Returns
    -------
    Quantity
        The total dust mass per species.
    """
    if snap.num_dust_species == 0:
        raise ValueError('No dust available')
    elif snap.num_dust_species > 0:
        if 'dust_fraction' in snap.available_arrays():
            h: Quantity = snap['smoothing_length']
            _mass: Quantity = snap['mass'][h > 0]
            dustfrac: Quantity = snap['dust_fraction'][h > 0]
            dustmass = (_mass[:, np.newaxis] * dustfrac).sum(axis=0)
            if squeeze:
                return dustmass.sum()
            return dustmass
        if snap.num_particles_of_type.get('dust') is not None:
            subsnaps: List[SnapLike] = snap.family('dust')  # type: ignore
            dustmass = [
                subsnap['mass'][subsnap['smoothing_length'] > 0].sum()
                for subsnap in subsnaps
            ]
            if squeeze:
                _dustmass = dustmass[0]
                for val in dustmass[1:]:
                    _dustmass = _dustmass + val
                return _dustmass
            return dustmass
    raise ValueError('Cannot calculate dust mass')


def gas_mass(snap: SnapLike) -> Quantity:
    """Calculate the total gas mass.

    Parameters
    ----------
    snap
        The Snap object.

    Returns
    -------
    Quantity
        The total gas mass.
    """
    gas: SubSnap = snap.family('gas')  # type: ignore
    if snap.num_dust_species == 0:
        return mass(gas, sinks=False)
    elif snap.num_dust_species > 0:
        if 'dust_fraction' in snap.available_arrays():
            h: Quantity = snap['smoothing_length']
            _mass: Quantity = snap['mass'][h > 0]
            dustfrac: Quantity = snap['dust_fraction'][h > 0]
            dustmass = (_mass[:, np.newaxis] * dustfrac).sum()
            mixturemass = mass(gas, sinks=False)
            return mixturemass - dustmass
        else:
            return mass(gas, sinks=False)
    raise ValueError('Cannot calculate gas mass')


def kinetic_energy(snap: SnapLike, sinks: bool = True) -> Quantity:
    """Calculate the total kinetic energy.

    Parameters
    ----------
    snap
        The Snap object.
    sinks : optional
        Include sink particles. Default is True.

    Returns
    -------
    Quantity
        The total kinetic energy.
    """
    ke = particles.kinetic_energy(snap=snap, ignore_accreted=True).sum()
    if sinks and snap.num_sinks > 0:
        if snap.num_sinks == 1:
            m, v = snap.sinks['mass'], snap.sinks['velocity']
            ke = ke + 1 / 2 * m * norm(v) ** 2
        elif snap.num_sinks > 1:
            ke = ke + particles.kinetic_energy(snap=snap.sinks).sum()
    return ke


def mass(snap: SnapLike, sinks: bool = True) -> Quantity:
    """Calculate the total mass.

    Parameters
    ----------
    snap
        The Snap object.
    sinks : optional
        Include sink particles. Default is True.

    Returns
    -------
    Quantity
        The total mass.
    """
    h: Quantity = snap['smoothing_length']
    mass: Quantity = snap['mass'][h > 0]

    mass_sum = mass.sum()
    if sinks and snap.num_sinks > 0:
        mass_sum = mass_sum + np.sum(snap.sinks['mass'])

    return mass_sum


def momentum(snap: SnapLike, sinks: bool = True) -> Quantity:
    """Calculate the total momentum.

    Parameters
    ----------
    snap
        The Snap object.
    sinks : optional
        Include sink particles. Default is True.

    Returns
    -------
    Quantity
        The total linear momentum like (px, py, pz).
    """
    mom = particles.momentum(snap=snap, ignore_accreted=True).sum(axis=0)
    if sinks and snap.num_sinks > 0:
        if snap.num_sinks == 1:
            mom = mom + snap.sinks['mass'] * snap.sinks['velocity']
        elif snap.num_sinks > 1:
            mass = snap.sinks['mass']
            velocity = snap.sinks['velocity']
            mom = mom + (mass[:, np.newaxis] * velocity).sum(axis=0)
    return mom


def specific_angular_momentum(
    snap: SnapLike, sinks: bool = True, origin: Quantity = None
) -> Quantity:
    """Calculate the total specific angular momentum.

    Parameters
    ----------
    snap
        The Snap object.
    sinks : optional
        Include sink particles. Default is True.
    origin : optional
        The origin around which to compute the angular momentum as a
        Quantity like (x, y, z) * au. Default is (0, 0, 0).

    Returns
    -------
    Quantity
        The total specific angular momentum on the particles like
        (hx, hy, hz).
    """
    angmom = particles.specific_angular_momentum(
        snap=snap, origin=origin, ignore_accreted=True
    ).sum(axis=0)
    if sinks and snap.num_sinks > 0:
        angmom_sinks = particles.specific_angular_momentum(
            snap=snap.sinks, origin=origin, ignore_accreted=False
        )
        if angmom_sinks.ndim == 1:
            angmom = angmom + angmom_sinks
        elif angmom_sinks.ndim == 2:
            angmom = angmom + angmom_sinks.sum(axis=0)
    return angmom


def specific_kinetic_energy(snap: SnapLike, sinks: bool = True) -> Quantity:
    """Calculate the total specific kinetic energy.

    Parameters
    ----------
    snap
        The Snap object.
    sinks : optional
        Include sink particles. Default is True.

    Returns
    -------
    Quantity
        The total specific kinetic energy.
    """
    ke = particles.specific_kinetic_energy(snap=snap, ignore_accreted=True).sum()
    if sinks and snap.num_sinks > 0:
        if snap.num_sinks == 1:
            v = snap.sinks['velocity']
            ke = ke + 1 / 2 * norm(v) ** 2
        elif snap.num_sinks > 1:
            ke = ke + particles.specific_kinetic_energy(snap=snap.sinks).sum()
    return ke
