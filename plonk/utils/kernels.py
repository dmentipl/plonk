"""SPH kernels."""

import numba
import numpy as np


@numba.njit
def kernel_cubic(q):
    """Cubic kernel function.

    The form of this function includes the "C_norm" factor. I.e.
    C_norm * f(q).

    Parameters
    ----------
    q
        The particle separation in units of smoothing length, i.e. r/h.

    Returns
    -------
    float
        C_norm * f(q) for the cubic kernel.
    """
    if q < 1:
        return (0.75 * q ** 3 - 1.5 * q ** 2 + 1) / np.pi
    elif q < 2:
        return (-0.25 * (q - 2) ** 3) / np.pi
    else:
        return 0.0


@numba.njit
def kernel_gradient_cubic(q):
    """Cubic kernel gradient function.

    The form of this function includes the "C_norm" factor. I.e.
    C_norm * f'(q).

    Parameters
    ----------
    q
        The particle separation in units of smoothing length, i.e. r/h.

    Returns
    -------
    float
        C_norm * f'(q) for the cubic kernel.
    """
    if q < 1:
        return q * (2.25 * q - 3.0) / np.pi
    elif q < 2:
        return -0.75 * (q - 2.0) ** 2 / np.pi
    else:
        return 0.0


@numba.njit
def kernel_quintic(q):
    """Quintic kernel function.

    The form of this function includes the "C_norm" factor. I.e.
    C_norm * f(q).

    Parameters
    ----------
    q
        The particle separation in units of smoothing length, i.e. r/h.

    Returns
    -------
    float
        C_norm * f(q) for the quintic kernel.
    """
    if q < 1:
        return (-10 * q ** 5 + 30 * q ** 4 - 60 * q ** 2 + 66) / (120 * np.pi)
    elif q < 2:
        return (-((q - 3) ** 5) + 6 * (q - 2) ** 5) / (120 * np.pi)
    elif q < 3:
        return (-((q - 3) ** 5)) / (120 * np.pi)
    else:
        return 0.0


@numba.njit
def kernel_gradient_quintic(q):
    """Quintic kernel gradient function.

    The form of this function includes the "C_norm" factor. I.e.
    C_norm * f'(q).

    Parameters
    ----------
    q
        The particle separation in units of smoothing length, i.e. r/h.

    Returns
    -------
    float
        C_norm * f'(q) for the quintic kernel.
    """
    if q < 1:
        return (q * (-50 * q ** 3 + 120 * q ** 2 - 120)) / (120 * np.pi)
    elif q < 2:
        return (-5 * (q - 3) ** 4 + 30 * (q - 2.0) ** 4) / (120 * np.pi)
    elif q < 3:
        return (-5 * (q - 3) ** 4) / (120 * np.pi)
    else:
        return 0.0


@numba.njit
def kernel_wendland_c4(q):
    """Wendland C4 kernel function.

    The form of this function includes the "C_norm" factor. I.e.
    C_norm * f(q).

    Parameters
    ----------
    q
        The particle separation in units of smoothing length, i.e. r/h.

    Returns
    -------
    float
        C_norm * f(q) for the Wendland C4 kernel.
    """
    if q < 2:
        return (
            ((-q / 2 + 1) ** 6 * (35 * q ** 2 / 12 + 3 * q + 1)) * 495 / (256 * np.pi)
        )
    else:
        return 0.0


@numba.njit
def kernel_gradient_wendland_c4(q):
    """Wendland C4 kernel gradient function.

    The form of this function includes the "C_norm" factor. I.e.
    C_norm * f'(q).

    Parameters
    ----------
    q
        The particle separation in units of smoothing length, i.e. r/h.

    Returns
    -------
    float
        C_norm * f'(q) for the Wendland C4 kernel.
    """
    if q < 2:
        return (
            (
                11.6666666666667 * q ** 2 * (0.5 * q - 1) ** 5
                + 4.66666666666667 * q * (0.5 * q - 1) ** 5
            )
            * 495
            / (256 * np.pi)
        )
    else:
        return 0.0


kernel_names = (
    'cubic',
    'quintic',
    'Wendland C4',
)

kernel_radius = {
    'cubic': 2.0,
    'quintic': 3.0,
    'Wendland C4': 2.0,
}

kernel_function = {
    'cubic': kernel_cubic,
    'quintic': kernel_quintic,
    'Wendland C4': kernel_wendland_c4,
}

kernel_gradient_function = {
    'cubic': kernel_gradient_cubic,
    'quintic': kernel_gradient_quintic,
    'Wendland C4': kernel_gradient_wendland_c4,
}
