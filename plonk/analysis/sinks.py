"""Calculate quantities on the sinks.

The following functions are available:

- kinetic_energy
- gravitational_potential_energy
- specific_orbital_energy
- specific_angular_momentum
- specific_angular_momentum
- eccentricity
- semi_major_axis
- inclination
- orbital_period
- mean_motion
- Roche_sphere
- Hill_radius
"""

import numpy as np


def kinetic_energy(m, v):
    """Calculate the kinetic energy.

    Parameters
    ----------
    m
        The mass of the body.
    v
        The velocity of the body.

    Returns
    -------
    The kinetic energy.
    """
    return 1 / 2 * m * _norm(v) ** 2


def gravitational_potential_energy(m1, m2, x1, x2, G):
    """Calculate the gravitational potential energy.

    Parameters
    ----------
    m1
        The mass of the first body.
    m2
        The mass of the second body.
    x1
        The position of the first body.
    x2
        The position of the second body.
    v1
        The velocity of the first body.
    v2
        The velocity of the second body.
    G
        The gravitational constant (in appropriate units).

    Returns
    -------
    The gravitational potential energy.
    """
    r = _norm(x1 - x2)
    return -G * m1 * m2 / r


def specific_orbital_energy(m1, m2, x1, x2, v1, v2, G):
    """Calculate the specific orbital energy.

    Parameters
    ----------
    m1
        The mass of the first body.
    m2
        The mass of the second body.
    x1
        The position of the first body.
    x2
        The position of the second body.
    v1
        The velocity of the first body.
    v2
        The velocity of the second body.
    G
        The gravitational constant (in appropriate units).

    Returns
    -------
    The specific orbital energy.
    """
    mu = m1 * m2 / (m1 + m2)
    ke = kinetic_energy(m1, v1) + kinetic_energy(m2, v2)
    pe = gravitational_potential_energy(m1, m2, x1, x2, G)
    return (ke + pe) / mu


def specific_angular_momentum(m1, m2, x1, x2, v1, v2):
    """Calculate the specific orbital energy.

    Parameters
    ----------
    m1
        The mass of the first body.
    m2
        The mass of the second body.
    x1
        The position of the first body.
    x2
        The position of the second body.
    v1
        The velocity of the first body.
    v2
        The velocity of the second body.

    Returns
    -------
    The specific orbital energy.
    """
    r = x1 - x2
    v = v1 - v2
    return np.cross(r, v)


def eccentricity(m1, m2, x1, x2, v1, v2, G):
    """Calculate the eccentricity.

    Parameters
    ----------
    m1
        The mass of the first body.
    m2
        The mass of the second body.
    x1
        The position of the first body.
    x2
        The position of the second body.
    v1
        The velocity of the first body.
    v2
        The velocity of the second body.
    G
        The gravitational constant (in appropriate units).

    Returns
    -------
    The eccentricity.
    """
    mu = G * (m1 + m2)
    eps = specific_orbital_energy(m1, m2, x1, x2, v1, v2, G)
    h = specific_angular_momentum(m1, m2, x1, x2, v1, v2)
    h_mag = _norm(h)
    return np.sqrt(1 + 2 * eps * h_mag ** 2 / mu ** 2)


def semi_major_axis(m1, m2, x1, x2, v1, v2, G):
    """Calculate the semi-major axis.

    Parameters
    ----------
    m1
        The mass of the first body.
    m2
        The mass of the second body.
    x1
        The position of the first body.
    x2
        The position of the second body.
    v1
        The velocity of the first body.
    v2
        The velocity of the second body.
    G
        The gravitational constant (in appropriate units).

    Returns
    -------
    The semi-major axis.
    """
    mu = G * (m1 + m2)
    h = specific_angular_momentum(m1, m2, x1, x2, v1, v2)
    h_mag = _norm(h)
    e = eccentricity(m1, m2, x1, x2, v1, v2, G)
    return h_mag ** 2 / (mu * (1 - e ** 2))


def inclination(m1, m2, x1, x2, v1, v2):
    """Calculate the inclination.

    Parameters
    ----------
    m1
        The mass of the first body.
    m2
        The mass of the second body.
    x1
        The position of the first body.
    x2
        The position of the second body.
    v1
        The velocity of the first body.
    v2
        The velocity of the second body.

    Returns
    -------
    The inclination.
    """
    h = specific_angular_momentum(m1, m2, x1, x2, v1, v2)
    h_z = h[..., 2]
    h_mag = _norm(h)
    return np.arccos(h_z / h_mag)


def orbital_period(m1, m2, x1, x2, v1, v2, G):
    """Calculate the orbital period.

    Parameters
    ----------
    m1
        The mass of the first body.
    m2
        The mass of the second body.
    x1
        The position of the first body.
    x2
        The position of the second body.
    v1
        The velocity of the first body.
    v2
        The velocity of the second body.
    G
        The gravitational constant (in appropriate units).

    Returns
    -------
    The orbital period.
    """
    a = semi_major_axis(m1, m2, x1, x2, v1, v2, G)
    mu = G * (m1 + m2)
    return 2 * np.pi * np.sqrt(a ** 3 / mu)


def mean_motion(m1, m2, x1, x2, v1, v2, G):
    """Calculate the mean motion.

    Parameters
    ----------
    m1
        The mass of the first body.
    m2
        The mass of the second body.
    x1
        The position of the first body.
    x2
        The position of the second body.
    v1
        The velocity of the first body.
    v2
        The velocity of the second body.
    G
        The gravitational constant (in appropriate units).

    Returns
    -------
    The mean motion.
    """
    P = orbital_period(m1, m2, x1, x2, v1, v2, G)
    return 2 * np.pi / P


def Roche_sphere(m1, m2, x1, x2):
    """Calculate an estimate of the Roche sphere.

    Uses the formula from Eggleton (1983) ApJ 268, 368-369.

    Parameters
    ----------
    m1
        The mass of the body around which to calculate the Roche sphere.
    m2
        The mass of the second body.
    x1
        The position of the first body.
    x2
        The position of the second body.

    Returns
    -------
    The Roche sphere radius.
    """
    separation = _norm(x1 - x2)
    q = m1 / m2
    return (
        separation
        * 0.49
        * q ** (2 / 3)
        / (0.6 * q ** (2 / 3) + np.log(1.0 + q ** (1 / 3)))
    )


def Hill_radius(M, m, X, x):
    """Calculate the Hill radius.

    This calculation assumes zero eccentricity, i.e. a circular orbit.

    Parameters
    ----------
    M
        The mass of the larger body.
    m
        The mass of the smaller body.
    X
        The position of the larger body.
    x
        The position of the smaller body.

    Returns
    -------
    The Hill sphere radius.
    """
    a = _norm(X - x)
    return a * (m / (3 * M)) ** (1 / 3)


def _norm(x):
    return np.sqrt(x[..., 0] ** 2 + x[..., 1] ** 2 + x[..., 2] ** 2)
