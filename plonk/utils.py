'''
utils.py

Daniel Mentiplay, 2019.
'''

import numpy as np

def density_from_smoothing_length(smoothingLength, particleMass, hfact=1.2):
    '''
    Calculate density from particle mass and smoothing length.
    '''

    return particleMass * (hfact/np.abs(smoothingLength))**3
