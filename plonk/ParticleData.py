'''
ParticleData.py

Daniel Mentiplay, 2019.
'''

import numpy as np

iGas  = 1
iDust = 7

class ParticleData:
    '''
    Generic SPH particles.
    '''

    def __init__(self):

        self.position        = None
        self.velocity        = None
        self.smoothingLength = None
        self.dustFrac        = None
        self.particleMass    = None
        self.particleType    = None

def density_from_smoothing_length(smoothingLength, particleMass, hfact=1.2):
    '''
    Calculate density from particle mass and smoothing length.
    '''

    return particleMass * (hfact/np.abs(smoothingLength))**3
