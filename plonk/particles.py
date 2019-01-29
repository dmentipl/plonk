'''
particles.py

Daniel Mentiplay, 2019.
'''

class Particles:
    '''
    Particles class represents Phantom particle arrays with the following
    properties:

        - position
        - velocity
        - smoothing length
        - dust fraction
    '''

    def __init__(self, position, velocity, smoothingLength, dustFrac):

        self.position = position
        self.velocity = velocity
        self.smoothingLength = smoothingLength
        self.dustFrac = dustFrac
