"""
utils.py

Daniel Mentiplay, 2019.
"""

import numpy as np


def print_warning(message):
    """Print formatted warning message."""
    message_length = len(message)
    message_length = min(message_length, 80)
    message = 'WARNING: ' + message
    print('\n' + message_length*'-'
               + '\n'
               + message
               + '\n'
               + message_length*'-')


def print_error(message):
    """Print formatted error message."""
    message_length = len(message)
    message_length = min(message_length, 80)
    message = 'ERROR: ' + message
    print('\n' + message_length*'-'
               + '\n'
               + message + '\n'
               + message_length*'-')


def rotate_vector_arbitrary_axis(u, v, theta):
    """Rotate a 3d vector (u) around an axis defined by another 3d
    vector (v) by an angle (theta) using the Rodrigues rotation formula.

    u can also be a numpy array of ndim 2, each row of which is a 3d
    vector.
    """

    norm = np.linalg.norm(v)
    if np.isclose(norm, 0.):
        k = v
    else:
        k = v / norm

    w = np.cross(k, u)
    dot_kv = np.sum(k*u, axis=1)
    dot_kv = np.stack(3*[dot_kv]).T

    return u*np.cos(theta) + w*np.sin(theta) + k*dot_kv*(1-np.cos(theta))


def normalize_vector(v):
    """Normalize a single vector."""

    norm = np.linalg.norm(v)
    if np.isclose(norm, 0.):
        return v
    return v / norm
