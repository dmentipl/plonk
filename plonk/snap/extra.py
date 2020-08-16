"""Add extra quantities to Snap."""

from numpy import ndarray

from .._units import units as plonk_units
from ..analysis import particles

ORIGIN = (0, 0, 0) * plonk_units.au


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
        vector = particles.vector_arrays['momentum']
        dust = particles.dust_arrays['momentum']

        @snap.add_array(vector=vector, dust=dust)
        def momentum(snap) -> ndarray:
            """Momentum."""
            return particles.momentum(snap=snap)

    if _add_array(snap, 'angular_momentum'):
        vector = particles.vector_arrays['angular_momentum']
        dust = particles.dust_arrays['angular_momentum']

        @snap.add_array(vector=vector, dust=dust)
        def angular_momentum(snap) -> ndarray:
            """Angular momentum."""
            origin = snap.translation if snap.translation is not None else ORIGIN
            return particles.angular_momentum(snap=snap, origin=origin)

    if _add_array(snap, 'specific_angular_momentum'):
        vector = particles.vector_arrays['specific_angular_momentum']
        dust = particles.dust_arrays['specific_angular_momentum']

        @snap.add_array(vector=vector, dust=dust)
        def specific_angular_momentum(snap) -> ndarray:
            """Specific angular momentum."""
            origin = snap.translation if snap.translation is not None else ORIGIN
            return particles.specific_angular_momentum(snap=snap, origin=origin)

    if _add_array(snap, 'kinetic_energy'):
        vector = particles.vector_arrays['kinetic_energy']
        dust = particles.dust_arrays['kinetic_energy']

        @snap.add_array(vector=vector, dust=dust)
        def kinetic_energy(snap) -> ndarray:
            """Kinetic energy."""
            return particles.kinetic_energy(snap=snap)

    if _add_array(snap, 'specific_kinetic_energy'):
        vector = particles.vector_arrays['specific_kinetic_energy']
        dust = particles.dust_arrays['specific_kinetic_energy']

        @snap.add_array(vector=vector, dust=dust)
        def specific_kinetic_energy(snap) -> ndarray:
            """Specific kinetic energy."""
            return particles.specific_kinetic_energy(snap=snap)

    if _add_array(snap, 'keplerian_frequency'):
        vector = particles.vector_arrays['keplerian_frequency']
        dust = particles.dust_arrays['keplerian_frequency']

        @snap.add_array(vector=vector, dust=dust)
        def keplerian_frequency(snap) -> ndarray:
            """Keplerian orbital frequency."""
            gravitational_parameter = snap.properties.get('gravitational_parameter')
            if gravitational_parameter is None:
                raise ValueError(
                    'To get Keplerian frequency, first set the gravitational parameter\n'
                    'via snap.set_gravitational_parameter.'
                )
            origin = snap.translation if snap.translation is not None else ORIGIN
            return particles.keplerian_frequency(
                snap=snap,
                gravitational_parameter=gravitational_parameter,
                origin=origin,
            )

    if _add_array(snap, 'semi_major_axis'):
        vector = particles.vector_arrays['semi_major_axis']
        dust = particles.dust_arrays['semi_major_axis']

        @snap.add_array(vector=vector, dust=dust)
        def semi_major_axis(snap) -> ndarray:
            """Semi-major axis."""
            gravitational_parameter = snap.properties.get('gravitational_parameter')
            if gravitational_parameter is None:
                raise ValueError(
                    'To get semi-major axis, first set the gravitational parameter\n'
                    'via snap.set_gravitational_parameter.'
                )
            origin = snap.translation if snap.translation is not None else ORIGIN
            return particles.semi_major_axis(
                snap=snap,
                gravitational_parameter=gravitational_parameter,
                origin=origin,
            )

    if _add_array(snap, 'eccentricity'):
        vector = particles.vector_arrays['eccentricity']
        dust = particles.dust_arrays['eccentricity']

        @snap.add_array(vector=vector, dust=dust)
        def eccentricity(snap) -> ndarray:
            """Eccentricity."""
            gravitational_parameter = snap.properties.get('gravitational_parameter')
            if gravitational_parameter is None:
                raise ValueError(
                    'To get eccentricity, first set the gravitational parameter\n'
                    'via snap.set_gravitational_parameter.'
                )
            origin = snap.translation if snap.translation is not None else ORIGIN
            return particles.eccentricity(
                snap=snap,
                gravitational_parameter=gravitational_parameter,
                origin=origin,
            )

    if _add_array(snap, 'inclination'):
        vector = particles.vector_arrays['inclination']
        dust = particles.dust_arrays['inclination']

        @snap.add_array(vector=vector, dust=dust)
        def inclination(snap) -> ndarray:
            """Inclination."""
            return particles.inclination(snap=snap)

    if _add_array(snap, 'radial_distance'):
        vector = particles.vector_arrays['radial_distance']
        dust = particles.dust_arrays['radial_distance']

        @snap.add_array(vector=vector, dust=dust)
        def radius_cylindrical(snap) -> ndarray:
            """Cylindrical radius."""
            return particles.radial_distance(snap=snap, coordinates='cylindrical')

        @snap.add_array(vector=vector, dust=dust)
        def radius_spherical(snap) -> ndarray:
            """Spherical radius."""
            return particles.radial_distance(snap=snap, coordinates='spherical')

    if _add_array(snap, 'azimuthal_angle'):
        vector = particles.vector_arrays['azimuthal_angle']
        dust = particles.dust_arrays['azimuthal_angle']

        @snap.add_array(vector=vector, dust=dust)
        def azimuthal_angle(snap) -> ndarray:
            """Azimuthal angle."""
            return particles.azimuthal_angle(snap=snap)

    if _add_array(snap, 'polar_angle'):
        vector = particles.vector_arrays['polar_angle']
        dust = particles.dust_arrays['polar_angle']

        @snap.add_array(vector=vector, dust=dust)
        def polar_angle(snap) -> ndarray:
            """Polar angle."""
            return particles.polar_angle(snap=snap)

    if _add_array(snap, 'radial_velocity'):
        vector = particles.vector_arrays['radial_velocity']
        dust = particles.dust_arrays['radial_velocity']

        @snap.add_array(vector=vector, dust=dust)
        def velocity_radial_cylindrical(snap) -> ndarray:
            """Cylindrical radial velocity."""
            return particles.radial_velocity(snap=snap, coordinates='cylindrical')

        @snap.add_array(vector=vector, dust=dust)
        def velocity_radial_spherical(snap) -> ndarray:
            """Spherical radial velocity."""
            return particles.radial_velocity(snap=snap, coordinates='spherical')

    if _add_array(snap, 'angular_velocity'):
        vector = particles.vector_arrays['angular_velocity']
        dust = particles.dust_arrays['angular_velocity']

        @snap.add_array(vector=vector, dust=dust)
        def angular_velocity(snap) -> ndarray:
            """Angular velocity."""
            return particles.angular_velocity(snap=snap)

    if _add_array(snap, 'temperature'):
        vector = particles.vector_arrays['temperature']
        dust = particles.dust_arrays['temperature']

        @snap.add_array(vector=vector, dust=dust)
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
        vector = particles.vector_arrays['gas_fraction']
        dust = particles.dust_arrays['gas_fraction']

        @snap.add_array(vector=vector, dust=dust)
        def gas_fraction(snap) -> ndarray:
            """Gas fraction."""
            return particles.gas_fraction(snap=snap)

    if _add_array(snap, 'gas_mass'):
        vector = particles.vector_arrays['gas_mass']
        dust = particles.dust_arrays['gas_mass']

        @snap.add_array(vector=vector, dust=dust)
        def gas_mass(snap) -> ndarray:
            """Gas mass."""
            return particles.gas_mass(snap=snap)

    if _add_array(snap, 'gas_density'):
        vector = particles.vector_arrays['gas_density']
        dust = particles.dust_arrays['gas_density']

        @snap.add_array(vector=vector, dust=dust)
        def gas_density(snap) -> ndarray:
            """Gas density."""
            return particles.gas_density(snap=snap)

    if _add_array(snap, 'dust_fraction'):
        vector = particles.vector_arrays['dust_fraction']
        dust = particles.dust_arrays['dust_fraction']

        if snap.properties['dust_method'] == 'dust as separate sets of particles':

            @snap.add_array(vector=vector, dust=dust)
            def dust_fraction(snap) -> ndarray:
                """Dust fraction."""
                return particles.dust_fraction(snap=snap)

    if _add_array(snap, 'dust_mass'):
        vector = particles.vector_arrays['dust_mass']
        dust = particles.dust_arrays['dust_mass']

        @snap.add_array(vector=vector, dust=dust)
        def dust_mass(snap) -> ndarray:
            """Dust mass."""
            return particles.dust_mass(snap=snap)

    if _add_array(snap, 'dust_density'):
        vector = particles.vector_arrays['dust_density']
        dust = particles.dust_arrays['dust_density']

        @snap.add_array(vector=vector, dust=dust)
        def dust_density(snap) -> ndarray:
            """Dust density."""
            return particles.dust_density(snap=snap)

    if _add_array(snap, 'stokes_number'):
        vector = particles.vector_arrays['stokes_number']
        dust = particles.dust_arrays['stokes_number']

        @snap.add_array(vector=vector, dust=dust)
        def stokes_number(snap) -> ndarray:
            """Stokes number."""
            gravitational_parameter = snap.properties.get('gravitational_parameter')
            if gravitational_parameter is None:
                raise ValueError(
                    'To get eccentricity, first set the gravitational parameter\n'
                    'via snap.set_gravitational_parameter.'
                )
            origin = snap.translation if snap.translation is not None else ORIGIN
            return particles.stokes_number(
                snap=snap,
                gravitational_parameter=gravitational_parameter,
                origin=origin,
            )
