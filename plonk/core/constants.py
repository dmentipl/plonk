"""
This module contains physical constants, units, and astronomical data
with values in cgs units.
"""

from collections import namedtuple

_constants = [
    'Boltzmann',
    'kB',
    'gas_constant',
    'R',
    'gravitational_constant',
    'G',
    'mass_of_hydrogen',
    'm_H',
    'speed_of_light',
    'c',
    'second',
    'minute',
    'hour',
    'day',
    'year',
    'micron',
    'micrometer',
    'millimeter',
    'meter',
    'kilometer',
    'km',
    'astronomical_unit',
    'au',
    'light_year',
    'ly',
    'parsec',
    'pc',
    'solar_mass',
    'solar_radius',
    'earth_mass',
    'earth_radius',
    'jupiter_mass',
    'saturn_mass',
    'uranus_mass',
    'neptune_mass',
]

_Constants = namedtuple('Constants', _constants)

_cd = dict()

# Physical constants
_cd['Boltzmann'] = _cd['kB'] = 1.380658e-16
_cd['gas_constant'] = _cd['R'] = 8.314e7
_cd['gravitational_constant'] = _cd['G'] = 6.67259e-8
_cd['mass_of_hydrogen'] = _cd['m_H'] = 1.6733e-24
_cd['speed_of_light'] = _cd['c'] = 2.997924e10

# Time units
_cd['second'] = 1.0
_cd['minute'] = 60.0
_cd['hour'] = 3600.0
_cd['day'] = 8.64e4
_cd['year'] = 3.1556926e7

# Length units
_cd['micron'] = _cd['micrometer'] = 1.0e-4
_cd['millimeter'] = 1.0e-1
_cd['meter'] = 100.0
_cd['kilometer'] = _cd['km'] = 1.0e5
_cd['astronomical_unit'] = _cd['au'] = 1.496e13
_cd['light_year'] = _cd['ly'] = 9.4605e17
_cd['parsec'] = _cd['pc'] = 3.086e18

# Solar system
_cd['solar_mass'] = 1.99e33
_cd['solar_radius'] = 6.959500e10
_cd['earth_mass'] = 5.972e27
_cd['earth_radius'] = 6.371315e8
_cd['jupiter_mass'] = 1.898e30
_cd['saturn_mass'] = 5.683e29
_cd['uranus_mass'] = 8.681e29
_cd['neptune_mass'] = 1.024e29

# Create constants namedtuple.
constants = _Constants(**_cd)
