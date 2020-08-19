"""Units."""

from pathlib import Path
from typing import Any, Dict, Union

import pint

from ._config import load_config

units = pint.UnitRegistry()
Quantity = units.Quantity


def add_units(config: Union[str, Path] = None):
    """Add units to the unit registry from a config file.

    Parameters
    ----------
    config : optional
        The path to a Plonk config.toml file. If None, use the default
        file.
    """
    conf = load_config(filename=config)
    for unit, definition in conf['units']['definitions'].items():
        units.define(f'{unit} = {definition}')


def array_units(config: Union[str, Path] = None) -> Dict[str, str]:
    """Return a dictionary of arrays with unit strings.

    Like the following:

        `{'density': 'kg / m ^ 3', ... 'position': 'm', ... }`

    Parameters
    ----------
    config : optional
        The path to a Plonk config.toml file. If None, use the default
        file.

    Returns
    -------
    Dict
    """
    conf = load_config(filename=config)
    d = dict()
    for key, val in conf['arrays']['dimensions'].items():
        dim = _convert_dim_string(val)
        if dim == {'[angle]': 1.0}:
            d[key] = "radian"
        else:
            for key1, val1 in conf['units']['defaults'].items():
                if _dimensionality_comparison(dict(units(val1).dimensionality), dim):
                    d[key] = val1
                    break
    return d


def array_quantities(config: Union[str, Path] = None) -> Dict[str, Any]:
    """Array dimensions dictionary.

    Like the following:

        `{
            'density': {'[mass]': 1, '[length]': -3},
            ...
            'position': {'[length]': 1},
            ...
        }`

    Parameters
    ----------
    config : optional
        The path to a Plonk config.toml file. If None, use the default
        file.

    Returns
    -------
    Dict
    """
    conf = load_config(filename=config)
    arrays = conf['arrays']['dimensions']
    dim = dict()
    for key, val in arrays.items():
        dim[key] = _convert_dim_string(val)
    return dim


def generate_array_code_units(code_units: Dict[str, Any]) -> Dict[str, Any]:
    """Generate array code units dictionary.

    Like the following:

        `{
            'density': 1.0 <Unit('kilogram / meter ** 3')>,
            ...
            'position': 149600000000.0 <Unit('meter')>,
            ...
        }`

    Parameters
    ----------
    code_units

    Returns
    -------
    Dict
        A dictionary of units as Pint quantities.
    """
    _units = dict()
    _array_quantities = array_quantities()
    for arr, unit in _array_quantities.items():
        _units[arr] = _get_code_unit(unit, code_units)
    return _units


def _convert_dim_string(string):
    if string == '':
        return {}
    list_dims = [dim.strip() for dim in string.split(',')]
    keys = list()
    vals = list()
    for dim in list_dims:
        key, val = dim.split(': ')
        val = float(val)
        keys.append(key)
        vals.append(val)
    return {'[' + key + ']': val for key, val in zip(keys, vals)}


def _get_code_unit(dim_dict, code_units):
    unit = 1.0 * units['dimensionless']
    for d in dim_dict:
        if d == '[angle]':
            unit *= units['radian']
        else:
            unit *= code_units[d[1:-1]] ** dim_dict[d]
    return unit


def _dimensionality_comparison(dim1, dim2):
    _dim1 = {key: float(val) for key, val in dim1.items()}
    _dim2 = {key: float(val) for key, val in dim2.items()}
    return _dim1 == _dim2
