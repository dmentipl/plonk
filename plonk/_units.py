"""Units."""

import copy

import pint

units = pint.UnitRegistry()
Quantity = units.Quantity

# Add units
units.define('solar_mass = 1.9891e30 kg')
units.define('solar_radius = 6.959500e8 m')
units.define('earth_mass = 5.979e24 kg')
units.define('earth_radius = 6.371315e6 m')
units.define('jupiter_mass = 1.89813e27 kg')

# Units dictionary
#   key is the array name
#   val is a tuple with
#     - the array type, e.g. 'length' or 'density'
#     - the array default unit, e.g. 'm' or 'kg / m ** 3'
units_dict = {
    'accretion_radius': ('length', 'm'),
    'alpha_viscosity_numerical': ('dimensionless', ''),
    'angular_momentum': ('angular_momentum', 'kg * m ** 2 / s'),
    'angular_velocity': ('velocity', 'm / s'),
    'azimuthal_angle': ('angle', 'rad'),
    'density': ('density', 'kg / m ** 3'),
    'differential_velocity': ('velocity', 'm / s'),
    'dust_density': ('density', 'kg / m ** 3'),
    'dust_fraction': ('dimensionless', ''),
    'dust_mass': ('mass', 'kg'),
    'dust_to_gas_ratio': ('dimensionless', ''),
    'eccentricity': ('dimensionless', ''),
    'gas_density': ('density', 'kg / m ** 3'),
    'gas_fraction': ('dimensionless', ''),
    'gas_mass': ('mass', 'kg'),
    'gravitational_potential': ('energy', 'J'),
    'inclination': ('angle', 'rad'),
    'internal_energy': ('energy_specific', 'J / kg'),
    'keplerian_frequency': ('frequency', 'Hz'),
    'kinetic_energy': ('energy', 'J'),
    'last_injection_time': ('time', 's'),
    'magnetic_field': ('magnetic_field', 'T'),
    'mass': ('mass', 'kg'),
    'mass_accreted': ('mass', 'kg'),
    'momentum': ('momentum', 'kg * m / s'),
    'polar_angle': ('angle', 'rad'),
    'position': ('length', 'm'),
    'pressure': ('pressure', 'Pa'),
    'projection': ('length', 'm'),
    'radius_cylindrical': ('length', 'm'),
    'radius_spherical': ('length', 'm'),
    'semi_major_axis': ('length', 'm'),
    'smoothing_length': ('length', 'm'),
    'softening_radius': ('length', 'm'),
    'sound_speed': ('velocity', 'm / s'),
    'specific_angular_momentum': ('angular_momentum_specific', 'm ** 2 / s'),
    'specific_kinetic_energy': ('energy_specific', 'J / kg'),
    'spin': ('angular_momentum', 'kg * m ** 2 / s'),
    'stokes_number': ('dimensionless', ''),
    'stopping_time': ('time', 's'),
    'sub_type': ('dimensionless', ''),
    'temperature': ('temperature', 'K'),
    'timestep': ('time', 's'),
    'type': ('dimensionless', ''),
    'velocity': ('velocity', 'm / s'),
    'velocity_divergence': ('frequency', 'Hz'),
    'velocity_radial_cylindrical': ('velocity', 'm / s'),
    'velocity_radial_spherical': ('velocity', 'm / s'),
}


def units_defaults():
    """Return a dictionary of arrays with unit strings.

    Like the following:

        `{'density': 'kg / m ** 3', ... 'position': 'm', ... }`

    This is useful for setting units for plots.
    """
    units_defaults_dict = {key: val[1] for key, val in units_dict.items()}
    return copy.copy(units_defaults_dict)


def generate_array_units_dict(units_dictionary):
    """Generate units array dictionary.

    Parameters
    ----------
    units_dictionary
        TODO: See generate_code_units_dict.

    Returns
    -------
    units
        A dictionary of units as Pint quantities.
    """
    _units = dict()
    units_type_dict = {key: val[0] for key, val in units_dict.items()}
    for arr, unit in units_type_dict.items():
        _units[arr] = units_dictionary[unit]
    return _units


def generate_code_units_dict(length, mass, time, temperature, magnetic_field):
    """Generate units dictionary.

    Parameters
    ----------
    length
        Length unit as a Pint quantity.
    mass
        Mass unit as a Pint quantity.
    time
        Time unit as a Pint quantity.
    temperature
        Temperature unit as a Pint quantity.
    magnetic_field
        Magnetic field unit as a Pint quantity.

    Returns
    -------
    units
        A dictionary of units as Pint quantities.
    """
    _units = dict()

    _units['acceleration'] = length / time ** 2
    _units['angle'] = units('radian')
    _units['angular_momentum'] = mass * length ** 2 / time
    _units['angular_momentum_specific'] = length ** 2 / time
    _units['density'] = mass / length ** 3
    _units['dimensionless'] = units('dimensionless')
    _units['energy'] = (mass * length ** 2 / time ** 2).to('joule')
    _units['energy_specific'] = (length ** 2 / time ** 2).to('joule/kg')
    _units['entropy'] = (mass * length ** 2 / time ** 2).to('joule') / temperature
    _units['force'] = (mass * length / time ** 2).to('newton')
    _units['frequency'] = (1 / time).to('hertz')
    _units['length'] = length
    _units['magnetic_field'] = magnetic_field.to('tesla')
    _units['mass'] = mass
    _units['momentum'] = mass * length / time
    _units['power'] = (mass * length ** 2 / time ** 3).to('watt')
    _units['pressure'] = (mass / time ** 2 / length).to('pascal')
    _units['temperature'] = 1.0 * temperature
    _units['time'] = time
    _units['velocity'] = length / time

    return _units
