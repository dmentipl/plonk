'''
units['py']

D. Mentiplay, 2019.
'''

class Units:
    '''
    This class represents the dump file units.

    Each value is the cgs value of the underlying physical unit.

    Arguments:
        udist : dist unit [cgs]
        umass : mass unit [cgs]
        utime : time unit [cgs]
    '''

    # pylint: disable=too-many-instance-attributes

    def __init__(self, udist=None, umass=None, utime=None):

        self.units = dict()

        self.units['distance'] = udist
        self.units['mass']     = umass
        self.units['time']     = utime

        self.units['frequency']        = 1 / self.units['time']

        self.units['velocity']         = self.units['distance'] \
                                       / self.units['time']

        self.units['momentum']         = self.units['mass'] \
                                       * self.units['distance'] \
                                       / self.units['time']

        self.units['force']            = self.units['mass'] \
                                       * self.units['distance'] \
                                       / self.units['time']**2

        self.units['pressure']         = self.units['mass'] \
                                       / (self.units['distance'] \
                                       * self.units['time']**2)

        self.units['energy']           = self.units['mass'] \
                                       * self.units['distance'] \
                                       / self.units['time']**2

        self.units['density']          = self.units['mass'] \
                                       / self.units['distance']**3

        self.units['surfaceDensity']   = self.units['mass'] \
                                       / self.units['distance']**2

        self.units['angularMomentum']  = self.units['mass'] \
                                       * self.units['distance']**2 \
                                       / self.units['time']

        self.units['torque']           = self.units['mass'] \
                                       * self.units['distance']**2 \
                                       / self.units['time']**2

def convert(quantity, unitFrom, unitTo):
    '''
    Convert a quantity from one unit to another:

        newValue = quantity * unitFrom / unitTo

    e.g. massInEarthMass = massInCodeUnits * units['mass'] / constants.earthMass

    Arguments:
        quantity : quantity to be converted
        unitFrom : current value of unit
        unitTo   : unit to convert quantity to
    '''

    return quantity * unitFrom / unitTo
