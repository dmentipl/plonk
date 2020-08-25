"""Add extra quantities to Snap."""

from .. import analysis

MODULES = {'common': analysis.particles, 'disc': analysis.discs}


def add_quantities(snap, category: str = 'common'):
    """Make extra quantities available on Snap.

    Parameters
    ----------
    snap
        The Snap object to add extra quantities to.
    category
        The category from which to add extra quantities. Options include:

        - 'common': adds common quantities appropriate for most
          simulations, such as angular momentum.
        - 'disc': adds accretion disc quantities, such as eccentricity
          or Keplerian frequency.

        Default is 'common'.
    """
    if category not in MODULES:
        raise ValueError(f'category must be in {list(MODULES)}')
    module = MODULES[category]

    for array in module.array_requires:  # type: ignore
        _add_array(snap, array, module)


def _can_add_array(snap, name, module):
    return set(module.array_requires[name]).issubset(snap.available_arrays())


def _vector(name, module):
    return True if name in module.vector_arrays else False


def _dust(name, module):
    return True if name in module.dust_arrays else False


def _add_array(snap, name, module):
    if _can_add_array(snap, name, module):
        vector = _vector(name, module)
        dust = _dust(name, module)
        snap._array_registry[name] = getattr(module, name)
        if vector is True:
            snap._vector_arrays.add(name)
        if dust is True:
            snap._dust_arrays.add(name)
