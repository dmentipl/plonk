'''
energy.py

Daniel Mentiplay, 2019.
'''

class Energy:
    '''
    Energy class represents Phantom simulation energy with the following
    properties:

        - kinetic
        - internal
        - potential
    '''

    def __init__(self, kinetic=None, internal=None, potential=None):

        self.kinetic = kinetic
        self.internal = internal
        self.potential = potential
