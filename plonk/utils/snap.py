"""Utils for snaps."""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING, List, Union

from .._config import read_config
from .._units import units as plonk_units

if TYPE_CHECKING:
    from ..snap.snap import SnapLike


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
    conf = read_config(filename=filename)
    for key, val in conf['arrays']['aliases'].items():
        snap.add_alias(key, val)


def dust_array_names(
    name: str, num_dust_species: int, add_gas: bool = False
) -> List[str]:
    """List dust array names.

    Parameters
    ----------
    name
        The base array name, e.g. "dust_density" or "stopping_time".
    num_dust_species
        The number of dust species.
    add_gas
        If True add the gas version of the dust name.

    Returns
    -------
    List
        A list of array names with appropriate suffixes.

    Examples
    --------
    Get the dust density strings.

    >>> dust_name_list('dust_density', 5)
    ['dust_density_001',
     'dust_density_002',
     'dust_density_003',
     'dust_density_004',
     'dust_density_005']

    Get the dust density strings with gas.

    >>> dust_name_list(name='dust_density', num_dust_species=5, add_gas=True)
    ['gas_density',
     'dust_density_001',
     'dust_density_002',
     'dust_density_003',
     'dust_density_004',
     'dust_density_005']
    """
    names = list()
    if add_gas:
        names.append(f'{name.replace("dust", "gas")}')
    names += [f'{name}_{n+1:03}' for n in range(num_dust_species)]
    return names


def vector_array_names(name: str, add_mag: bool = False) -> List[str]:
    """List vector array names.

    Parameters
    ----------
    name
        The base array name, e.g. "angular_momentum".
    add_mag
        If True add the magnitude of the array.

    Returns
    -------
    List
        A list of array names with appropriate suffixes.

    Examples
    --------
    Get the angular momentum strings.

    >>> vector_name_list('angular_momentum')
    ['angular_momentum_x',
     'angular_momentum_y',
     'angular_momentum_z']

    Get the angular momentum strings with magnitude.

    >>> vector_name_list(name='angular_momentum', add_mag=True)
    ['angular_momentum_x',
     'angular_momentum_y',
     'angular_momentum_z',
     'angular_momentum_mag']
    """
    names = [f'{name}_{x}' for x in ['x', 'y', 'z']]
    if add_mag:
        names.append(f'{name}_mag')
    return names
