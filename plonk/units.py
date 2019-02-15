'''
units.py

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

    def __init__(self, udist=1., umass=1., utime=1.):

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

        self.units['surface_density']   = self.units['mass'] \
                                       / self.units['distance']**2

        self.units['angular_momentum']  = self.units['mass'] \
                                       * self.units['distance']**2 \
                                       / self.units['time']

        self.units['torque']           = self.units['mass'] \
                                       * self.units['distance']**2 \
                                       / self.units['time']**2

def convert(quantity, unit_from, unit_to):
    '''
    Convert a quantity from one unit to another:

        new_value = quantity * unit_from / unit_to

    e.g. mass_in_earth_mass = \\
            mass_in_code_units * units['mass'] / constants.earth_mass

    Arguments:
        quantity  : quantity to be converted
        unit_from : current value of unit
        unit_to   : unit to convert quantity to
    '''

    return quantity * unit_from / unit_to
