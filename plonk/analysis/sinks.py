"""Sink quantities.

Calculate various quantities relating to sink particles.
"""

from typing import Union

import numpy as np

from ..snap.snap import Snap, SubSnap

SnapLike = Union[Snap, SubSnap]


def Roche_sphere(m1: float, m2: float, separation: float):
    """Calculate an estimate of the Roche sphere.

    Uses the formula from Eggleton (1983) ApJ 268, 368-369.

    Parameters
    ----------
    m1
        The mass of the body around which to calculate the Roche sphere.
    m2
        The mass of the second body.

    """
    q = m1 / m2
    return (
        separation
        * 0.49
        * q ** (2 / 3)
        / (0.6 * q ** (2 / 3) + np.log(1.0 + q ** (1 / 3)))
    )
