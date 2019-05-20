"""
This module contains functions for physical units.
"""

from collections import namedtuple


def units_in_cgs(udist=1.0, umass=1.0, utime=1.0):
    """
    Physical units in cgs.

    Each value is the cgs value of the underlying physical unit.

    Parameters
    ----------
    udist : float
        The distance unit [cgs].
    umass : float
        The mass unit [cgs].
    utime : float
        The time unit [cgs].

    Returns
    -------
    namedtuple
        This contains the units of different physical quantities in cgs.
    """

    quantities = [
        'angular_momentum',
        'distance',
        'energy',
        'force',
        'frequency',
        'mass',
        'mass_density',
        'momentum',
        'pressure',
        'surface_density',
        'time',
        'torque',
        'velocity',
    ]

    Units = namedtuple('Units', quantities)

    units = dict()
    units['distance'] = udist
    units['time'] = utime
    units['mass'] = umass
    units['angular_momentum'] = umass * udist ** 2 / utime
    units['energy'] = umass * udist / utime ** 2
    units['force'] = umass * udist / utime ** 2
    units['frequency'] = 1 / utime
    units['mass_density'] = umass / udist ** 3
    units['momentum'] = umass * udist / utime
    units['pressure'] = umass / (udist * utime ** 2)
    units['surface_density'] = umass / udist ** 2
    units['torque'] = umass * udist ** 2 / utime ** 2
    units['velocity'] = udist / utime

    return Units(**units)


def convert_units(quantity, unit_from, unit_to):
    """
    Convert a quantity from one unit to another:

        new_value = quantity * unit_from / unit_to

    e.g. mass_in_earth_mass =
            mass_in_code_units * units['mass'] / constants.earth_mass

    Parameters
    ----------
    quantity : float
        The quantity to be converted.
    unit_from : float
        The current value of unit in cgs.
    unit_to : float
        The unit to convert quantity to in cgs.

    Returns
    -------
    float
        The original quantity in the new units.
    """

    return quantity * unit_from / unit_to
