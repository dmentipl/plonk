"""
This module contains the Units class.
"""

from collections import namedtuple

import numpy as np

from .constants import constants

_quantities = [
    ('length', 'L', 'L'),
    ('time', 'T', 'T'),
    ('mass', 'M', 'M'),
    ('angular_momentum', 'J', 'M L^2 T^-1'),
    ('density', 'rho', 'M L^-3'),
    ('energy', 'E', 'M L^2 T^-2'),
    ('force', 'F', 'M L T^-1'),
    ('frequency', 'f', 'T^-1'),
    ('momentum', 'p', 'M L T^-1'),
    ('pressure', 'P', 'M L^-1 T^-1'),
    ('torque', 'tau', 'M L^2 T^-2'),
    ('velocity', 'v', 'L T^-1'),
]

_quantities_long_name = [q[0] for q in _quantities]
_quantities_short_name = [q[1] for q in _quantities]
_quantities_base = [q[2] for q in _quantities]
_quantities_list = _quantities_long_name + _quantities_short_name

_Units = namedtuple(
    'Units', _quantities_list, defaults=(None,) * len(_quantities_list)
)

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
    length : float or str
        If float, the length unit in cgs. If str, must be in list of
        available quantities.
    mass : float or str
        If float, the mass unit in cgs. If str, must be in list of
        available quantities.
    time : float or str
        If float, the time unit in cgs. If str, must be in list of
        available quantities.

    Other parameters
    ----------------
    **kwargs
        Other quantities to set units for.
    """

    def __init__(self, length=None, mass=None, time=None, **kwargs):

        if length is None:
            if self.length is not None:
                length = self.length
            else:
                length = 'cm'
        if mass is None:
            if self.mass is not None:
                mass = self.mass
            else:
                mass = 'g'
        if time is None:
            if self.time is not None:
                time = self.time
            else:
                time = 's'

        if isinstance(length, str):
            length_str = length.lower()
            if length_str in [unit[0] for unit in LENGTH_UNITS]:
                length = [
                    unit[1] for unit in LENGTH_UNITS if unit[0] == length_str
                ][0]
            else:
                raise ValueError(f'{length} is not available')
        elif isinstance(length, (int, float)):
            length_str = ''
            for unit in LENGTH_UNITS:
                if np.isclose(length, unit[1]):
                    length_str = unit[0]
        else:
            raise ValueError('Cannot determine length unit')

        if isinstance(mass, str):
            mass_str = mass.lower()
            if mass_str in [unit[0] for unit in MASS_UNITS]:
                mass = [unit[1] for unit in MASS_UNITS if unit[0] == mass_str][
                    0
                ]
            else:
                raise ValueError(f'{mass} is not available')
        elif isinstance(mass, (int, float)):
            mass_str = ''
            for unit in MASS_UNITS:
                if np.isclose(mass, unit[1]):
                    mass_str = unit[0]
        else:
            raise ValueError('Cannot determine mass unit')

        if isinstance(time, str):
            time_str = time.lower()
            if time_str in [unit[0] for unit in TIME_UNITS]:
                time = [unit[1] for unit in TIME_UNITS if unit[0] == time_str][
                    0
                ]
            else:
                raise ValueError(f'{time} is not available')
        elif isinstance(time, (int, float)):
            time_str = ''
            for unit in TIME_UNITS:
                if np.isclose(time, unit[1]):
                    time_str = unit[0]
        else:
            raise ValueError('Cannot determine time unit')

        self.length = (length, length_str)
        self.mass = (mass, mass_str)
        self.time = (time, time_str)

        _ud = dict()

        _ud['length'] = _ud['L'] = length
        _ud['time'] = _ud['T'] = time
        _ud['mass'] = _ud['M'] = mass

        for quantity, value in kwargs.items():
            if quantity not in _quantities_list:
                raise ValueError(f'{quantity} unavailable')
            if not isinstance(value, float):
                raise TypeError(f'{quantity} must be float')
            if quantity in _quantities_long_name:
                quantity_long = quantity
                quantity_short = _quantities_short_name[
                    _quantities_long_name.index(quantity)
                ]
            elif quantity in _quantities_short_name:
                quantity_short = quantity
                quantity_long = _quantities_long_name[
                    _quantities_short_name.index(quantity)
                ]
            _ud[quantity_short] = value
            _ud[quantity_long] = value

        self.units = _Units(**_ud)

    def convert_quantity_to_new_units(self, quantity, dimension, new_units):
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
        Convert a mass to new units.
        >>> new_units = Units(umass='earth_mass')
        >>> units.convert_quantity_to_new_units(
                mass, 'M', new_units
            )
        """

        return self.convert_quantity_to_cgs(
            quantity, dimension
        ) / new_units.convert_quantity_to_cgs(1.0, dimension)

    def convert_quantity(self, quantity, dimension, unit_factor):
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
        unit_factor : float
            The factor to convert the quantity as cgs.

        Returns
        -------
        float
            The original quantity in the new units.

        Examples
        --------
        Convert a mass from Units.units to Earth masses
        >>> units.convert_quantity_to_new_units(
                mass, 'M', plonk.constants.earth_mass
            )
        """

        return self.convert_quantity_to_cgs(quantity, dimension) * unit_factor

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
            cgs_val = getattr(self.units, dimension)
            if cgs_val is not None:
                return quantity * cgs_val
            else:
                if dimension in _quantities_short_name:
                    dimension = _quantities_base[
                        _quantities_short_name.index(dimension)
                    ]
                elif dimension in _quantities_long_name:
                    dimension = _quantities_base[
                        _quantities_long_name.index(dimension)
                    ]
        return quantity * self._get_cgs_from_dimension(dimension)

    def _get_cgs_from_dimension(self, expression):
        d = dimensions_as_dict(expression)
        val = 1.0
        for key in d:
            val *= getattr(self.units, key) ** d[key]
        return val

    def _units_same(self, a, b):
        if a in _quantities_short_name:
            a = _quantities_base[_quantities_short_name.index(a)]
        elif a in _quantities_long_name:
            a = _quantities_base[_quantities_long_name.index(a)]
        elif a in _quantities_base:
            pass
        elif dimensions_as_dict(a):
            pass
        else:
            raise ValueError(f'{a} unknown dimension/unit')

        if b in _quantities_short_name:
            b = _quantities_base[_quantities_short_name.index(b)]
        elif b in _quantities_long_name:
            b = _quantities_base[_quantities_long_name.index(b)]
        elif b in _quantities_base:
            pass
        elif dimensions_as_dict(b):
            pass
        else:
            raise ValueError(f'{b} unknown dimension/unit')

        return dimensions_as_dict(a) == dimensions_as_dict(b)


def dimensions_as_dict(expression):
    """
    Convert expression like 'L^m M^n T^o' to {'L': m, 'M': n, 'T': o}.
    """

    units = list()
    for unit in expression.split():
        unit = unit.split('^')
        if unit[0] not in ('L', 'M', 'T'):
            raise ValueError(
                'Cannot interpret string: must be like ' '"L^m M^n T^o"'
            )
        if len(unit) == 1:
            units.append([unit[0], 1])
        elif len(unit) == 2:
            try:
                i = int(unit[1])
            except ValueError:
                raise ValueError(
                    'Cannot interpret string: must be like ' '"L^m M^n T^o"'
                )
            units.append([unit[0], i])
        else:
            raise ValueError('Cannot interpret string')

    for dim in ('L', 'M', 'T'):
        if dim not in [unit[0] for unit in units]:
            units.append([dim, 0])

    return dict(sorted(units))
