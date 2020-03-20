"""Physical units on snapshots."""

from .. import units as plonk_units
from .snap import SnapLike

# Add units
plonk_units.define('solar_mass = 1.9891e33 g')
plonk_units.define('solar_radius = 6.959500e10 cm')
plonk_units.define('earth_mass = 5.979e27 g')
plonk_units.define('earth_radius = 6.371315e8 cm')
plonk_units.define('jupiter_mass = 1.89813e30 g')


def generate_units_dictionary(length, mass, time, magnetic_field):
    """Generate units dictionary.

    Parameters
    ----------
    length
        Length unit as a Pint quantity.
    mass
        Mass unit as a Pint quantity.
    time
        Time unit as a Pint quantity.
    magnetic_field
        Magnetic field unit as a Pint quantity.

    Returns
    -------
    units
        A dictionary of units as Pint quantities.
    """
    units = {}

    units['dimensionless'] = length / length
    units['length'] = length
    units['time'] = time
    units['mass'] = mass
    units['magnetic_field'] = magnetic_field
    units['frequency'] = 1 / time
    units['velocity'] = length / time
    units['momentum'] = mass * length / time
    units['angular_momentum'] = mass * length ** 2 / time
    units['specific_angular_momentum'] = length ** 2 / time
    units['density'] = mass / length ** 3
    units['acceleration'] = length / time ** 2
    units['force'] = mass * length / time ** 2
    units['energy'] = mass * length ** 2 / time ** 2
    units['specific_energy'] = length ** 2 / time ** 2
    units['pressure'] = mass / time ** 2 / length
    units['temperature'] = plonk_units.kelvin

    return units


def gravitational_constant_in_code_units(snap: SnapLike) -> float:
    """Gravitational constant in code units.

    Parameters
    ----------
    snap
        The Snap object.

    Returns
    -------
    float
        The gravitational constant in code units.
    """
    G = plonk_units.newtonian_constant_of_gravitation
    G_units = snap.units['length'] ** 3 / snap.units['mass'] / snap.units['time'] ** 2
    G = (G / G_units).to_base_units().magnitude
    return G
