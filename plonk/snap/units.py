"""Physical units on snapshots."""


def generate_units_dictionary(length, mass, time):
    """Generate units dictionary.

    Parameters
    ----------
    length
        Length unit as a Pint quantity.
    mass
        Mass unit as a Pint quantity.
    time
        Time unit as a Pint quantity.

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
    units['frequency'] = 1 / time
    units['velocity'] = length / time
    units['momentum'] = mass * length / time
    units['angular_momentum'] = mass * length ** 2 / time
    units['density'] = mass / length ** 3
    units['acceleration'] = length / time ** 2
    units['force'] = mass * length / time ** 2
    units['energy'] = mass * length ** 2 / time ** 2

    return units
