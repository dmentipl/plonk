'''
sinks.py

Daniel Mentiplay, 2019.
'''

class Sinks:
    '''
    Sinks class represents Phantom sink particles with the following properties:

        - mass
        - accretionRadius
        - spin
        - position
        - velocity
    '''

    def __init__(self, position=None, velocity=None, accretionRadius=None,
                 mass=None, spin=None):

        self.position = position
        self.velocity = velocity
        self.accretionRadius = accretionRadius
        self.mass = mass
        self.spin = spin
