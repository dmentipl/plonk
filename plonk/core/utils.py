"""
utils.py

Daniel Mentiplay, 2019.
"""

import os

import numpy as np


def rotate_vector_arbitrary_axis(u, v, theta):
    """Rotate a 3d vector (u) around an axis defined by another 3d
    vector (v) by an angle (theta) using the Rodrigues rotation formula.

    u can also be a numpy array of ndim 2, each row of which is a 3d
    vector.
    """

    norm = np.linalg.norm(v)
    if np.isclose(norm, 0.0):
        k = v
    else:
        k = v / norm

    w = np.cross(k, u)
    dot_kv = np.sum(k * u, axis=1)
    dot_kv = np.stack(3 * [dot_kv]).T

    return (
        u * np.cos(theta) + w * np.sin(theta) + k * dot_kv * (1 - np.cos(theta))
    )


def normalize_vector(v):
    """Normalize a single vector."""

    norm = np.linalg.norm(v)
    if np.isclose(norm, 0.0):
        return v
    return v / norm


def file_size(filepath):
    """
    Get file size.

    Parameters
    ----------
    filepath : str
        Path to file. Can be relative.

    Returns
    -------
    size : int
        Size of file in units specified by `unit`.
    unit : str
        Unit of `size` return.
    """

    size_in_bytes = os.stat(filepath).st_size

    return nbytes_to_human(size_in_bytes)


def nbytes_to_human(size_in_bytes):
    """
    Convert a number of bytes to human readable form.

    Parameters
    ----------
    size_in_bytes : int
        Number of bytes.

    Returns
    -------
    size : int
        Size of file in units specified by `unit`.
    unit : str
        Unit of `size` return.
    """

    size_names = ('B', 'KB', 'MB', 'GB', 'TB', 'PB')
    size_maxes = (1e3, 1e6, 1e9, 1e12, 1e15)

    for size_name, size_max in zip(size_names, size_maxes):
        if size_in_bytes < size_max:
            return (int(size_in_bytes * 1e3 / size_max), size_name)

    return int(size_in_bytes / size_maxes[-1]), size_names[-1]
