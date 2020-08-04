"""Add extra quantities to Snap."""

from numpy import ndarray

from ..analysis import particles


def _add_array(snap, name):
    return set(particles.array_requires[name]).issubset(snap.available_arrays())


def extra_quantities(snap):
    """Make extra quantities available.

    Parameters
    ----------
    snap
        The Snap object to add extra quantities to.
    """
    if _add_array(snap, 'momentum'):
        rotatable = particles.array_rotatable['momentum']
        dust = particles.array_dust['momentum']

        @snap.add_array(rotatable=rotatable, dust=dust)
        def momentum(snap) -> ndarray:
            """Momentum."""
            return particles.momentum(snap=snap)

    if _add_array(snap, 'angular_momentum'):
        rotatable = particles.array_rotatable['angular_momentum']
        dust = particles.array_dust['angular_momentum']

        @snap.add_array(rotatable=rotatable, dust=dust)
        def angular_momentum(snap) -> ndarray:
            """Angular momentum."""
            origin = (
                snap.translation if snap.translation is not None else (0.0, 0.0, 0.0)
            )
            return particles.angular_momentum(snap=snap, origin=origin)

    if _add_array(snap, 'specific_angular_momentum'):
        rotatable = particles.array_rotatable['specific_angular_momentum']
        dust = particles.array_dust['specific_angular_momentum']

        @snap.add_array(rotatable=rotatable, dust=dust)
        def specific_angular_momentum(snap) -> ndarray:
            """Specific angular momentum."""
            origin = (
                snap.translation if snap.translation is not None else (0.0, 0.0, 0.0)
            )
            return particles.specific_angular_momentum(snap=snap, origin=origin)

    if _add_array(snap, 'kinetic_energy'):
        rotatable = particles.array_rotatable['kinetic_energy']
        dust = particles.array_dust['kinetic_energy']

        @snap.add_array(rotatable=rotatable, dust=dust)
        def kinetic_energy(snap) -> ndarray:
            """Kinetic energy."""
            return particles.kinetic_energy(snap=snap)

    if _add_array(snap, 'specific_kinetic_energy'):
        rotatable = particles.array_rotatable['specific_kinetic_energy']
        dust = particles.array_dust['specific_kinetic_energy']

        @snap.add_array(rotatable=rotatable, dust=dust)
        def specific_kinetic_energy(snap) -> ndarray:
            """Specific kinetic energy."""
            return particles.specific_kinetic_energy(snap=snap)

    if _add_array(snap, 'keplerian_frequency'):
        rotatable = particles.array_rotatable['keplerian_frequency']
        dust = particles.array_dust['keplerian_frequency']

        @snap.add_array(rotatable=rotatable, dust=dust)
        def keplerian_frequency(snap) -> ndarray:
            """Keplerian orbital frequency."""
            gravitational_parameter = snap.properties.get('gravitational_parameter')
            if gravitational_parameter is None:
                raise ValueError(
                    'To get Keplerian frequency, first set the gravitational parameter\n'
                    'via snap.set_gravitational_parameter.'
                )
            origin = (
                snap.translation if snap.translation is not None else (0.0, 0.0, 0.0)
            )
            return particles.keplerian_frequency(
                snap=snap,
                gravitational_parameter=gravitational_parameter,
                origin=origin,
            )

    if _add_array(snap, 'semi_major_axis'):
        rotatable = particles.array_rotatable['semi_major_axis']
        dust = particles.array_dust['semi_major_axis']

        @snap.add_array(rotatable=rotatable, dust=dust)
        def semi_major_axis(snap) -> ndarray:
            """Semi-major axis."""
            gravitational_parameter = snap.properties.get('gravitational_parameter')
            if gravitational_parameter is None:
                raise ValueError(
                    'To get semi-major axis, first set the gravitational parameter\n'
                    'via snap.set_gravitational_parameter.'
                )
            origin = (
                snap.translation if snap.translation is not None else (0.0, 0.0, 0.0)
            )
            return particles.semi_major_axis(
                snap=snap,
                gravitational_parameter=gravitational_parameter,
                origin=origin,
            )

    if _add_array(snap, 'eccentricity'):
        rotatable = particles.array_rotatable['eccentricity']
        dust = particles.array_dust['eccentricity']

        @snap.add_array(rotatable=rotatable, dust=dust)
        def eccentricity(snap) -> ndarray:
            """Eccentricity."""
            gravitational_parameter = snap.properties.get('gravitational_parameter')
            if gravitational_parameter is None:
                raise ValueError(
                    'To get eccentricity, first set the gravitational parameter\n'
                    'via snap.set_gravitational_parameter.'
                )
            origin = (
                snap.translation if snap.translation is not None else (0.0, 0.0, 0.0)
            )
            return particles.eccentricity(
                snap=snap,
                gravitational_parameter=gravitational_parameter,
                origin=origin,
            )

    if _add_array(snap, 'inclination'):
        rotatable = particles.array_rotatable['inclination']
        dust = particles.array_dust['inclination']

        @snap.add_array(rotatable=rotatable, dust=dust)
        def inclination(snap) -> ndarray:
            """Inclination."""
            return particles.inclination(snap=snap)

    if _add_array(snap, 'radial_distance'):
        rotatable = particles.array_rotatable['radial_distance']
        dust = particles.array_dust['radial_distance']

        @snap.add_array(rotatable=rotatable, dust=dust)
        def radius_cylindrical(snap) -> ndarray:
            """Cylindrical radius."""
            return particles.radial_distance(snap=snap, coordinates='cylindrical')

        @snap.add_array(rotatable=rotatable, dust=dust)
        def radius_spherical(snap) -> ndarray:
            """Spherical radius."""
            return particles.radial_distance(snap=snap, coordinates='spherical')

    if _add_array(snap, 'azimuthal_angle'):
        rotatable = particles.array_rotatable['azimuthal_angle']
        dust = particles.array_dust['azimuthal_angle']

        @snap.add_array(rotatable=rotatable, dust=dust)
        def azimuthal_angle(snap) -> ndarray:
            """Azimuthal angle."""
            return particles.azimuthal_angle(snap=snap)

    if _add_array(snap, 'polar_angle'):
        rotatable = particles.array_rotatable['polar_angle']
        dust = particles.array_dust['polar_angle']

        @snap.add_array(rotatable=rotatable, dust=dust)
        def polar_angle(snap) -> ndarray:
            """Polar angle."""
            return particles.polar_angle(snap=snap)

    if _add_array(snap, 'radial_velocity'):
        rotatable = particles.array_rotatable['radial_velocity']
        dust = particles.array_dust['radial_velocity']

        @snap.add_array(rotatable=rotatable, dust=dust)
        def velocity_radial_cylindrical(snap) -> ndarray:
            """Cylindrical radial velocity."""
            return particles.radial_velocity(snap=snap, coordinates='cylindrical')

        @snap.add_array(rotatable=rotatable, dust=dust)
        def velocity_radial_spherical(snap) -> ndarray:
            """Spherical radial velocity."""
            return particles.radial_velocity(snap=snap, coordinates='spherical')

    if _add_array(snap, 'angular_velocity'):
        rotatable = particles.array_rotatable['angular_velocity']
        dust = particles.array_dust['angular_velocity']

        @snap.add_array(rotatable=rotatable, dust=dust)
        def angular_velocity(snap) -> ndarray:
            """Angular velocity."""
            return particles.angular_velocity(snap=snap)

    if _add_array(snap, 'temperature'):
        rotatable = particles.array_rotatable['temperature']
        dust = particles.array_dust['temperature']

        @snap.add_array(rotatable=rotatable, dust=dust)
        def temperature(snap) -> ndarray:
            """Temperature."""
            molecular_weight = snap.properties.get('molecular_weight')
            if molecular_weight is None:
                raise ValueError(
                    'To get temperature, first set the molecular weight parameter\n'
                    'via snap.set_molecular_weight method.'
                )
            return particles.temperature(snap=snap, molecular_weight=molecular_weight)

    if _add_array(snap, 'gas_fraction'):
        rotatable = particles.array_rotatable['gas_fraction']
        dust = particles.array_dust['gas_fraction']

        @snap.add_array(rotatable=rotatable, dust=dust)
        def gas_fraction(snap) -> ndarray:
            """Gas fraction."""
            return particles.gas_fraction(snap=snap)

    if _add_array(snap, 'gas_mass'):
        rotatable = particles.array_rotatable['gas_mass']
        dust = particles.array_dust['gas_mass']

        @snap.add_array(rotatable=rotatable, dust=dust)
        def gas_mass(snap) -> ndarray:
            """Gas mass."""
            return particles.gas_mass(snap=snap)

    if _add_array(snap, 'gas_density'):
        rotatable = particles.array_rotatable['gas_density']
        dust = particles.array_dust['gas_density']

        @snap.add_array(rotatable=rotatable, dust=dust)
        def gas_density(snap) -> ndarray:
            """Gas density."""
            return particles.gas_density(snap=snap)

    if _add_array(snap, 'dust_fraction'):
        rotatable = particles.array_rotatable['dust_fraction']
        dust = particles.array_dust['dust_fraction']

        if snap.properties['dust_method'] == 'dust as separate sets of particles':

            @snap.add_array(rotatable=rotatable, dust=dust)
            def dust_fraction(snap) -> ndarray:
                """Dust fraction."""
                return particles.dust_fraction(snap=snap)

    if _add_array(snap, 'dust_mass'):
        rotatable = particles.array_rotatable['dust_mass']
        dust = particles.array_dust['dust_mass']

        @snap.add_array(rotatable=rotatable, dust=dust)
        def dust_mass(snap) -> ndarray:
            """Dust mass."""
            return particles.dust_mass(snap=snap)

    if _add_array(snap, 'dust_density'):
        rotatable = particles.array_rotatable['dust_density']
        dust = particles.array_dust['dust_density']

        @snap.add_array(rotatable=rotatable, dust=dust)
        def dust_density(snap) -> ndarray:
            """Dust density."""
            return particles.dust_density(snap=snap)

    if _add_array(snap, 'stokes_number'):
        rotatable = particles.array_rotatable['stokes_number']
        dust = particles.array_dust['stokes_number']

        @snap.add_array(rotatable=rotatable, dust=dust)
        def stokes_number(snap) -> ndarray:
            """Stokes number."""
            gravitational_parameter = snap.properties.get('gravitational_parameter')
            if gravitational_parameter is None:
                raise ValueError(
                    'To get eccentricity, first set the gravitational parameter\n'
                    'via snap.set_gravitational_parameter.'
                )
            origin = (
                snap.translation if snap.translation is not None else (0.0, 0.0, 0.0)
            )
            return particles.stokes_number(
                snap=snap,
                gravitational_parameter=gravitational_parameter,
                origin=origin,
            )
