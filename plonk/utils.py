'''
utils.py

Daniel Mentiplay, 2019.
'''

import numpy as np

def print_warning(message):
    '''
    Print formatted warning message.
    '''
    l = len(message)
    l = min(l, 80)
    message = 'WARNING: ' + message
    print( '\n' \
         + l*'-' + '\n' \
         + message + '\n' \
         + l*'-' )

def print_error(message):
    '''
    Print formatted error message.
    '''
    l = len(message)
    l = min(l, 80)
    message = 'ERROR: ' + message
    print( '\n' \
         + l*'-' + '\n' \
         + message + '\n' \
         + l*'-' )

def rotate_vector_arbitrary_axis(u, v, theta):
    '''
    Rotate a 3d vector (u) around an axis defined by another 3d vector (v)
    by an angle (theta) using the Rodrigues rotation formula.

    u can also be a numpy array of ndim 2, each row of which is a 3d vector.
    '''

    norm = np.linalg.norm(v)
    if np.isclose(norm, 0.):
        k = v
    else:
        k = v / norm

    w = np.cross(k, u)

    k = np.stack(u.shape[0]*[k])
    dot_kv = np.sum(k*u, axis=1)
    dot_kv = np.stack(3*[dot_kv]).T

    return u*np.cos(theta) + w*np.sin(theta) + k*dot_kv*(1-np.cos(theta))
