"""
This module contains functions for physical units.
"""

from collections import namedtuple

_quantities = [
    'length',
    'L',
    'time',
    'T',
    'mass',
    'M',
    'angular_momentum',
    'J',
    'energy',
    'E',
    'force',
    'F',
    'frequency',
    'f',
    'mass_density',
    'rho',
    'momentum',
    'p',
    'pressure',
    'P',
    'surface_density',
    'sigma',
    'torque',
    'tau',
    'velocity',
    'v',
]

_Units = namedtuple('Units', _quantities)


class Units:
    def __init__(self, ulength=None, umass=None, utime=None):
        """
        Units in cgs.

        Parameters
        ----------
        ulength : float
            The length unit [cgs].
        umass : float
            The mass unit [cgs].
        utime : float
            The time unit [cgs].
        """

        if ulength is None:
            ulength = 1.0
        if utime is None:
            utime = 1.0
        if umass is None:
            umass = 1.0

        self.set_units(ulength, umass, utime)

    def set_units(self, ulength=None, umass=None, utime=None):
        """
        Set units in cgs.

        Parameters
        ----------
        ulength : float
            The length unit [cgs].
        umass : float
            The mass unit [cgs].
        utime : float
            The time unit [cgs].
        """

        if ulength is None:
            ulength = self.units.length
        if utime is None:
            utime = self.units.time
        if umass is None:
            umass = self.units.mass

        _ud = dict()

        _ud['length'] = _ud['L'] = ulength
        _ud['time'] = _ud['T'] = utime
        _ud['mass'] = _ud['M'] = umass

        _ud['angular_momentum'] = _ud['J'] = umass * ulength ** 2 / utime
        _ud['energy'] = _ud['E'] = umass * ulength / utime ** 2
        _ud['force'] = _ud['F'] = umass * ulength / utime ** 2
        _ud['frequency'] = _ud['f'] = 1 / utime
        _ud['mass_density'] = _ud['rho'] = umass / ulength ** 3
        _ud['momentum'] = _ud['p'] = umass * ulength / utime
        _ud['pressure'] = _ud['P'] = umass / (ulength * utime ** 2)
        _ud['surface_density'] = _ud['sigma'] = umass / ulength ** 2
        _ud['torque'] = _ud['tau'] = umass * ulength ** 2 / utime ** 2
        _ud['velocity'] = _ud['v'] = ulength / utime

        self.units = _Units(**_ud)

    def convert_quantity(self, quantity, dimension, new_unit_in_cgs):
        """
        Convert a quantity to new units.

        Parameters
        ----------
        quantity : float
            The quantity to be converted.
        dimension : str
            This can be a string available in units.units, e.g.
            'velocity', or 'energy'. Alternatively it can be a
            combination of 'L', 'M', 'T' with powers separated by
            spaces, e.g. 'L^3 M^-1 T^-2' for the gravitational
            constant.
        new_unit_in_cgs : float
            The unit to convert quantity to expressed in cgs.

        Returns
        -------
        float
            The original quantity in the new units.

        Examples
        --------
        Convert a mass from Units.units to Earth masses
        >>> units.convert_quantity(mass, 'M', plonk.constants.earth_mass)
        """

        return (
            self.convert_quantity_to_cgs(quantity, dimension) / new_unit_in_cgs
        )

    def convert_quantity_to_cgs(self, quantity, dimension):
        """
        Convert quantity from current units to cgs.

        Parameters
        ----------
        quantity : float
            The quantity (in current units) to covert to cgs.
        dimension : str
            This can be a string available in units.units, e.g.
            'velocity', or 'energy'. Alternatively it can be a
            combination of 'L', 'M', 'T' with powers separated by
            spaces, e.g. 'L^3 M^-1 T^-2' for the gravitational
            constant.

        Returns
        -------
        float
            The value of the quantity in cgs.
        """
        if dimension in self.units._fields:
            return quantity * getattr(self.units, dimension)
        return quantity * self._get_cgs_from_dimension(dimension)

    def _get_cgs_from_dimension(self, expression):
        d = _get_dimension_from_string(expression)
        val = 1.0
        for key in d:
            val *= getattr(self.units, key) ** d[key]
        return val


def _get_dimension_from_string(expression):
    units = [unit.split('^') for unit in expression.split()]
    d = {}
    for unit in units:
        if len(unit) > 1:
            d[unit[0]] = int(unit[1])
        elif len(unit) == 1:
            d[unit[0]] = 1
    return d
