'''
parameters.py

Daniel Mentiplay, 2019.
'''

class Parameters:
    '''
    Phantom simulation parameters.
    '''

    def __init__(self):

        class Particles:
            '''
            Particle parameters.
            '''
            def __init__(self):
                self.ntypes = 17
                self.npartoftype = list()
                self.massoftype = list()

        class Dust:
            '''
            Dust parameters.
            '''
            def __init__(self):
                self.ndustsmall = 0
                self.ndustlarge = 0
                self.grainSize = list()
                self.grainDens = list()

        class Sink:
            '''
            Sink parameters.
            '''
            def __init__(self):
                self.nSinks = 0

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
                self.alphau = 1.0

        self.particles = Particles()
        self.dust      = Dust()
        self.sink      = Sink()
        self.eos       = EOS()
        self.numerical = Numerical()
