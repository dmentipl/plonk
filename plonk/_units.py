"""Units."""

import pint

units = pint.UnitRegistry(system='cgs')
Quantity = units.Quantity

# Add units
units.define('solar_mass = 1.9891e33 g')
units.define('solar_radius = 6.959500e10 cm')
units.define('earth_mass = 5.979e27 g')
units.define('earth_radius = 6.371315e8 cm')
units.define('jupiter_mass = 1.89813e30 g')


def generate_units_dictionary(length, mass, time, magnetic_field):
    """Generate units dictionary.

    Parameters
    ----------
    length
        Length unit as a Pint quantity.
    mass
        Mass unit as a Pint quantity.
    time
        Time unit as a Pint quantity.
    magnetic_field
        Magnetic field unit as a Pint quantity.

    Returns
    -------
    units
        A dictionary of units as Pint quantities.
    """
    _units = {}

    _units['dimensionless'] = length / length
    _units['length'] = length
    _units['time'] = time
    _units['mass'] = mass
    _units['magnetic_field'] = magnetic_field
    _units['frequency'] = 1 / time
    _units['velocity'] = length / time
    _units['momentum'] = mass * length / time
    _units['angular_momentum'] = mass * length ** 2 / time
    _units['specific_angular_momentum'] = length ** 2 / time
    _units['density'] = mass / length ** 3
    _units['acceleration'] = length / time ** 2
    _units['force'] = mass * length / time ** 2
    _units['energy'] = mass * length ** 2 / time ** 2
    _units['specific_energy'] = length ** 2 / time ** 2
    _units['pressure'] = mass / time ** 2 / length
    _units['temperature'] = units.kelvin

    return _units
