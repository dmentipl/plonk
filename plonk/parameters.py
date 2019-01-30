'''
parameters.py

Daniel Mentiplay, 2019.
'''

from .units import Units

class Parameters:
    '''
    Phantom simulation parameters.
    '''

    def __init__(self):

        dust = dict()

        dust['nDustSmall'] = 0
        dust['nDustLarge'] = 0
        dust['grainSize'] = list()
        dust['grainDens'] = list()

        eos = dict()

        eos['ieos'] = 3
        eos['gamma'] = 1.
        eos['polyk'] = 1.
        eos['qfacdisc'] = 0.5

        numerical = dict()

        numerical['tolh'] = 0.0001
        numerical['C_cour'] = 0.3
        numerical['C_force'] = 0.25
        numerical['alpha'] = 0.1

        self.dust      = dust
        self.eos       = eos
        self.numerical = numerical
        self.units     = Units()
