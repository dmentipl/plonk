"""
This module contains physical constants, units, and astronomical data.

The class _Constants is instantiated in the object constants, which
is what should be imported.
"""


class _Constants:
    """This class contains physical constants in cgs units."""

    def __init__(self):

        # Physical constants
        self.Boltzmann = 1.380658e-16
        self.gas_constant = 8.314e7
        self.gravitational_constant = 6.67259e-8
        self.mass_of_hydrogen = 1.6733e-24
        self.speed_of_light = 2.997924e10

        # Time units
        self.seconds = 1.0e0
        self.minutes = 6.0e1
        self.hours = 3.6e3
        self.days = 8.64e4
        self.years = 3.1556926e7

        # Distance units
        self.micron = 1.0e-4
        self.millimeter = 1.0e-1
        self.kilometer = 1.0e5
        self.astronomical_unit = 1.496e13
        self.light_year = 9.4605e17
        self.parsec = 3.086e18

        # Solar system
        self.earth_mass = 5.972e27
        self.neptune_mass = 1.024e29
        self.uranus_mass = 8.681e29
        self.saturn_mass = 5.683e29
        self.jupiter_mass = 1.898e30
        self.solar_mass = 1.99e33
        self.earth_radius = 6.371315e8
        self.solar_radius = 6.959500e10


# Instantiate Constants class to create global constants object.
constants = _Constants()
