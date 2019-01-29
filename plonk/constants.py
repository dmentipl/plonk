'''
constants.py

D. Mentiplay, 2019.
'''

class Constants:
    '''
    This class contains physical constants in cgs units.

        Physical constants
        ------------------

        gravitationalConstant
        speedLight
        boltzmann
        massHydrogen
        gasConstant


        Time scales
        -----------

        seconds
        minutes
        hours
        days
        years


        Distance scales
        ---------------

        astronomicalUnit
        lightYear
        parsec
        kilometer
        micron
        millimeter


        Solar system
        ------------

        solarMass
        solarRadius

        jupiterMass
        saturnMass
        uranusMass
        neptuneMass

        earthMass
        earthRadius
    '''

    # pylint: disable=too-many-instance-attributes

    def __init__(self):

        #--- Physical constants

        self.gravitationalConstant = 6.67259e-8
        self.speedLight            = 2.997924e10
        self.boltzmann             = 1.380658e-16
        self.massHydrogen          = 1.6733e-24
        self.gasConstant           = 8.314e7

        #--- Time scales

        self.seconds               = 1.e0
        self.minutes               = 6.0e1
        self.hours                 = 3.6e3
        self.days                  = 8.64e4
        self.years                 = 3.1556926e7

        #--- Distance scales

        self.astronomicalUnit      = 1.496e13
        self.lightYear             = 9.4605e17
        self.parsec                = 3.086e18
        self.kilometer             = 1.e5
        self.micron                = 1.e-4
        self.millimeter            = 1.e-1

        #--- Solar system

        self.solarMass             = 1.99e33
        self.solarRadius           = 6.959500e10

        self.jupiterMass           = 1.898e30
        self.saturnMass            = 5.683e29
        self.uranusMass            = 8.681e29
        self.neptuneMass           = 1.024e29

        self.earthMass             = 5.972e27
        self.earthRadius           = 6.371315e8

# Instantiate Constants class to create global constants object.
constants = Constants()
