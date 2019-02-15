'''
constants.py

D. Mentiplay, 2019.
'''

class _Constants:
    '''
    This class contains physical constants in cgs units.

        Physical constants
        ------------------

        gravitational_constant
        speed_light
        boltzmann
        mass_hydrogen
        gas_constant


        Time scales
        -----------

        seconds
        minutes
        hours
        days
        years


        Distance scales
        ---------------

        astronomical_unit
        light_year
        parsec
        kilometer
        micron
        millimeter


        Solar system
        ------------

        solar_mass
        solar_radius

        jupiter_mass
        saturn_mass
        uranus_mass
        neptune_mass

        earth_mass
        earth_radius
    '''

    # pylint: disable=too-many-instance-attributes

    def __init__(self):

        #--- Physical constants

        self.gravitational_constant = 6.67259e-8
        self.speed_light            = 2.997924e10
        self.boltzmann              = 1.380658e-16
        self.mass_hydrogen          = 1.6733e-24
        self.gas_constant           = 8.314e7

        #--- Time scales

        self.seconds                = 1.e0
        self.minutes                = 6.0e1
        self.hours                  = 3.6e3
        self.days                   = 8.64e4
        self.years                  = 3.1556926e7

        #--- Distance scales

        self.astronomical_unit      = 1.496e13
        self.light_year             = 9.4605e17
        self.parsec                 = 3.086e18
        self.kilometer              = 1.e5
        self.micron                 = 1.e-4
        self.millimeter             = 1.e-1

        #--- Solar system

        self.solar_mass             = 1.99e33
        self.solar_radius           = 6.959500e10

        self.jupiter_mass           = 1.898e30
        self.saturn_mass            = 5.683e29
        self.uranus_mass            = 8.681e29
        self.neptune_mass           = 1.024e29

        self.earth_mass             = 5.972e27
        self.earth_radius           = 6.371315e8

# Instantiate Constants class to create global constants object.
constants = _Constants()
