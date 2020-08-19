"""Units."""

from pathlib import Path
from typing import Union

import pint

from ._config import load_config

units = pint.UnitRegistry()
Quantity = units.Quantity

# Add useful astronomical units
units.define('solar_mass = 1.9891e30 kg')
units.define('solar_radius = 6.959500e8 m')
units.define('earth_mass = 5.979e24 kg')
units.define('earth_radius = 6.371315e6 m')
units.define('jupiter_mass = 1.89813e27 kg')


def array_units(filename: Union[str, Path] = None):
    """Return a dictionary of arrays with unit strings.

    Like the following:

        `{'density': 'kg / m ^ 3', ... 'position': 'm', ... }`

    This is useful for setting units for plots.
    """
    if filename is None:
        conf = load_config()
    else:
        conf = load_config(filename=filename)
    d = dict()
    for key, val in conf['arrays'].items():
        dim = _convert_dim_string(val)
        if dim == {'[angle]': 1.0}:
            d[key] = "radian"
        else:
            for key1, val1 in conf['units'].items():
                if _dimensionality_comparison(dict(units(val1).dimensionality), dim):
                    d[key] = val1
                    break
    return d


def array_quantities(filename: Union[str, Path] = None):
    """TODO."""
    if filename is None:
        config = load_config()
    else:
        config = load_config(filename=filename)
    arrays = config['arrays']
    dim = dict()
    for key, val in arrays.items():
        dim[key] = _convert_dim_string(val)
    return dim


def generate_array_code_units(code_units):
    """Generate array code units dictionary.

    Parameters
    ----------
    code_units

    Returns
    -------
    units
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
