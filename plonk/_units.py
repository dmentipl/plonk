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

        `{'density': 'kg / m ** 3', ... 'position': 'm', ... }`

    This is useful for setting units for plots.
    """
    if filename is None:
        conf = load_config()
    else:
        conf = load_config(filename=filename)
    d = dict()
    for key, val in conf['arrays'].items():
        d[key] = conf['units'][val]
    return d


def array_quantities(filename: Union[str, Path] = None):
    """TODO."""
    if filename is None:
        config = load_config()
    else:
        config = load_config(filename=filename)
    return config['arrays']


def generate_array_units_dict(units_dictionary):
    """Generate units array dictionary.

    Parameters
    ----------
    units_dictionary
        TODO: See generate_code_units_dict.

    Returns
    -------
    units
        A dictionary of units as Pint quantities.
    """
    _units = dict()
    _array_quantities = array_quantities()
    for arr, unit in _array_quantities.items():
        _units[arr] = units_dictionary[unit]
    return _units


def generate_code_units_dict(length, mass, time, temperature, magnetic_field):
    """Generate units dictionary.

    Parameters
    ----------
    length
        Length unit as a Pint quantity.
    mass
        Mass unit as a Pint quantity.
    time
        Time unit as a Pint quantity.
    temperature
        Temperature unit as a Pint quantity.
    magnetic_field
        Magnetic field unit as a Pint quantity.

    Returns
    -------
    units
        A dictionary of units as Pint quantities.
    """
    _units = dict()

    _units['acceleration'] = length / time ** 2
    _units['angle'] = units('radian')
    _units['angular_momentum'] = mass * length ** 2 / time
    _units['angular_momentum_specific'] = length ** 2 / time
    _units['density'] = mass / length ** 3
    _units['dimensionless'] = units('dimensionless')
    _units['energy'] = (mass * length ** 2 / time ** 2).to('joule')
    _units['energy_specific'] = (length ** 2 / time ** 2).to('joule/kg')
    _units['entropy'] = (mass * length ** 2 / time ** 2).to('joule') / temperature
    _units['force'] = (mass * length / time ** 2).to('newton')
    _units['frequency'] = (1 / time).to('hertz')
    _units['length'] = length
    _units['magnetic_field'] = magnetic_field.to('tesla')
    _units['mass'] = mass
    _units['momentum'] = mass * length / time
    _units['power'] = (mass * length ** 2 / time ** 3).to('watt')
    _units['pressure'] = (mass / time ** 2 / length).to('pascal')
    _units['temperature'] = 1.0 * temperature
    _units['time'] = time
    _units['velocity'] = length / time

    return _units
