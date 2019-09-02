"""
Miscellaneous utilites.

For example, transformations on vectors.
"""

import numpy as np


def rotate_vector_arbitrary_axis(u, v, theta):
    """
    Rotate a vector around an arbitrary axis.

    Rotate a 3d vector (u) around an axis defined by another 3d
    vector (v) by an angle (theta) using the Rodrigues rotation formula.

    Parameters
    ----------
    u : (3,) or (3,N) ndarray
        The vector, or collection of vectors, to rotate.
    v : (3,) ndarray
        The axis around which to rotate.
    theta : float
        The angle of the rotation.

    Returns
    -------
    (3,) or (3,N) ndarray
        The rotated vector.
    """

    norm = np.linalg.norm(v)
    if np.isclose(norm, 0.0):
        k = v
    else:
        k = v / norm

    w = np.cross(k, u)

    if u.ndim == 2:
        dot_kv = np.sum(k * u, axis=1)
    elif u.ndim == 1:
        dot_kv = np.sum(k * u)
    else:
        raise ValueError('u has incorrect dimensions')

    dot_kv = np.stack(3 * [dot_kv]).T

    return u * np.cos(theta) + w * np.sin(theta) + k * dot_kv * (1 - np.cos(theta))


def normalize_vector(v):
    """
    Normalize a single vector.

    Parameters
    ----------
    u : (N,) ndarray
        The vector to normalize.

    Returns
    -------
    (N,) ndarray
        The normalized vector.
    """

    norm = np.linalg.norm(v)
    if np.isclose(norm, 0.0):
        return v
    return v / norm
