'''
utils.py

Daniel Mentiplay, 2019.
'''

import numpy as np

hfact = 1.2

def density_from_smoothing_length(smoothingLength, particleMass):
    '''
    Calculate density from particle mass and smoothing length.
    '''

    return particleMass * (hfact/np.abs(smoothingLength))**3
