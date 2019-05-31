"""
This module contains functions for physical units.
"""

from collections import namedtuple

import numpy as np

from .constants import constants

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

LENGTH_UNITS = (
    ('cm', 1.0),
    ('km', constants.km),
    ('solar_radius', constants.solar_radius),
    ('r_sun', constants.solar_radius),
    ('astronomical_unit', constants.au),
    ('au', constants.au),
    ('parsec', constants.pc),
    ('pc', constants.pc),
)

MASS_UNITS = (
    ('g', 1.0),
    ('earth_mass', constants.earth_mass),
    ('m_earth', constants.earth_mass),
    ('mearth', constants.earth_mass),
    ('earthm', constants.earth_mass),
    ('solar_mass', constants.solar_mass),
    ('m_sun', constants.solar_mass),
    ('msolar', constants.solar_mass),
    ('solarm', constants.solar_mass),
)

TIME_UNITS = (
    ('ms', 1.0e-3),
    ('millisecond', 1.0e-3),
    ('s', 1.0),
    ('second', 1.0),
    ('hr', constants.hour),
    ('hour', constants.hour),
    ('day', constants.day),
    ('yr', constants.year),
    ('year', constants.year),
    ('kyr', 1.0e3 * constants.year),
    ('myr', 1.0e6 * constants.year),
)


class Units:
    """
    Units in cgs.

    Parameters
    ----------
    ulength : float or str
        If float, the length unit in cgs. If str, must be in list of
        available quantities.
    umass : float or str
        If float, the mass unit in cgs. If str, must be in list of
        available quantities.
    utime : float or str
        If float, the time unit in cgs. If str, must be in list of
        available quantities.
    """
    def __init__(self, ulength=None, umass=None, utime=None):

        self.length = None
        self.time = None
        self.mass = None

        self.set_units(ulength, umass, utime)

    def set_units(self, ulength=None, umass=None, utime=None):
        """
        Set units in cgs.

        Parameters
        ----------
        ulength : float or str
            If float, the length unit in cgs. If str, must be in list of
            available quantities.
        umass : float or str
            If float, the mass unit in cgs. If str, must be in list of
            available quantities.
        utime : float or str
            If float, the time unit in cgs. If str, must be in list of
            available quantities.
        """

        if ulength is None:
            if self.length is not None:
                ulength = self.length
            else:
                ulength = 'cm'
        if umass is None:
            if self.mass is not None:
                umass = self.mass
            else:
                umass = 'g'
        if utime is None:
            if self.time is not None:
                utime = self.time
            else:
                utime = 's'

        if isinstance(ulength, str):
            ulength_str = ulength.lower()
            if ulength_str in [unit[0] for unit in LENGTH_UNITS]:
                ulength = [
                    unit[1] for unit in LENGTH_UNITS if unit[0] == ulength_str
                ][0]
            else:
                raise ValueError(f'{ulength} is not available')
        elif isinstance(ulength, (int, float)):
            ulength_str = ''
            for unit in LENGTH_UNITS:
                if np.isclose(ulength, unit[1]):
                    ulength_str = unit[0]
        else:
            raise ValueError('Cannot determine length unit')

        if isinstance(umass, str):
            umass_str = umass.lower()
            if umass_str in [unit[0] for unit in MASS_UNITS]:
                umass = [
                    unit[1] for unit in MASS_UNITS if unit[0] == umass_str
                ][0]
            else:
                raise ValueError(f'{umass} is not available')
        elif isinstance(umass, (int, float)):
            umass_str = ''
            for unit in MASS_UNITS:
                if np.isclose(umass, unit[1]):
                    umass_str = unit[0]
        else:
            raise ValueError('Cannot determine mass unit')

        if isinstance(utime, str):
            utime_str = utime.lower()
            if utime_str in [unit[0] for unit in TIME_UNITS]:
                utime = [
                    unit[1] for unit in TIME_UNITS if unit[0] == utime_str
                ][0]
            else:
                raise ValueError(f'{utime} is not available')
        elif isinstance(utime, (int, float)):
            utime_str = ''
            for unit in TIME_UNITS:
                if np.isclose(utime, unit[1]):
                    utime_str = unit[0]
        else:
            raise ValueError('Cannot determine time unit')

        self.length = (ulength, ulength_str)
        self.mass = (umass, umass_str)
        self.time = (utime, utime_str)

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

    def convert_quantity(self, quantity, dimension, new_units):
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
        new_units : Units
            The units to convert quantity to as a Units object.

        Returns
        -------
        float
            The original quantity in the new units.

        Examples
        --------
        Convert a mass from Units.units to Earth masses
        >>> units.convert_quantity(mass, 'M', plonk.constants.earth_mass)
        """

        return self.convert_quantity_to_cgs(
            quantity, dimension
        ) / new_units.convert_quantity_to_cgs(1.0, dimension)

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
    if len(units) > 3:
        raise ValueError('Cannot interpret string')
    if not set([unit[0] for unit in units]).issubset(set(('L', 'T', 'M'))):
        raise ValueError('Cannot interpret string')
    d = {}
    for unit in units:
        if len(unit) > 1:
            d[unit[0]] = int(unit[1])
        elif len(unit) == 1:
            d[unit[0]] = 1
    return d
