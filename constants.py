'''
constants.py

D. Mentiplay, 2019.
'''

class Constants:
    '''
    This class contains physical constants in cgs units.

        G                : gravitational constant [cm^3 s^-2 g^-1]
        kB               : Boltzmann constant [erg/K]

        year             : year [s]
        au               : astronomical unit [cm]
        solarMass        : solar mass [g]
        jupiterMass      : Jupiter mass [g]
        earthMass        : Earth mass [g]
    '''

    # pylint: disable=too-many-instance-attributes

    def __init__(self):

        self.G = 6.67259e-8
        self.kB = 1.380658e-16
        self.mH = 1.6733e-24

        self.year = 31557600
        self.au = 1.496e13
        self.solarMass = 1.99e33
        self.jupiterMass = 1.898e30
        self.earthMass = 5.972e27

# Instantiate Constants class to create global constants object.
constants = Constants()
