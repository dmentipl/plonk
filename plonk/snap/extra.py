"""Add extra quantities to Snap."""

from ..analysis import particles


def extra_quantities(snap):
    """Make extra quantities available.

    Parameters
    ----------
    snap
        The Snap object to add extra quantities to.
    """
    for array in particles.array_requires:
        _add_array(snap, array)


def _can_add_array(snap, name):
    return set(particles.array_requires[name]).issubset(snap.available_arrays())


def _vector(name):
    return True if name in particles.vector_arrays else False


def _dust(name):
    return True if name in particles.dust_arrays else False


def _add_array(snap, name):
    if _can_add_array(snap, name):
        vector = _vector(name)
        dust = _dust(name)
        snap._array_registry[name] = getattr(particles, name)
        if vector is True:
            snap._vector_arrays.add(name)
        if dust is True:
            snap._dust_arrays.add(name)
