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

        class Dust:
            '''
            Dust parameters.
            '''
            def __init__(self):
                self.ndustsmall = 0
                self.ndustlarge = 0
                self.grainSize = list()
                self.grainDens = list()

        class EOS:
            '''
            Equation of state parameters.
            '''
            def __init__(self):
                self.ieos = 3
                self.gamma = 1.
                self.polyk = 1.
                self.qfacdisc = 0.5

        class Numerical:
            '''
            Numerical parameters.
            '''
            def __init__(self):
                self.tolh = 0.0001
                self.C_cour = 0.3
                self.C_force = 0.25
                self.alpha = 0.1

        self.dust      = Dust()
        self.eos       = EOS()
        self.numerical = Numerical()
        self.units     = Units()
