"""Add extra quantities to Snap."""

from numpy import ndarray

from .._units import units as plonk_units
from ..analysis import particles

ORIGIN = (0, 0, 0) * plonk_units.au


def _add_array(snap, name):
    return set(particles.array_requires[name]).issubset(snap.available_arrays())


def _vector(name):
    return True if name in particles.vector_arrays else False


def _dust(name):
    return True if name in particles.dust_arrays else False


def extra_quantities(snap):
    """Make extra quantities available.

    Parameters
    ----------
    snap
        The Snap object to add extra quantities to.
    """
    if _add_array(snap, 'momentum'):
        vector = _vector('momentum')
        dust = _dust('momentum')

        @snap.add_array(vector=vector, dust=dust)
        def momentum(snap) -> ndarray:
            """Momentum."""
            return particles.momentum(snap=snap)

    if _add_array(snap, 'angular_momentum'):
        vector = _vector('angular_momentum')
        dust = _dust('angular_momentum')

        @snap.add_array(vector=vector, dust=dust)
        def angular_momentum(snap) -> ndarray:
            """Angular momentum."""
            origin = snap.translation if snap.translation is not None else ORIGIN
            return particles.angular_momentum(snap=snap, origin=origin)

    if _add_array(snap, 'specific_angular_momentum'):
        vector = _vector('specific_angular_momentum')
        dust = _dust('specific_angular_momentum')

        @snap.add_array(vector=vector, dust=dust)
        def specific_angular_momentum(snap) -> ndarray:
            """Specific angular momentum."""
            origin = snap.translation if snap.translation is not None else ORIGIN
            return particles.specific_angular_momentum(snap=snap, origin=origin)

    if _add_array(snap, 'kinetic_energy'):
        vector = _vector('kinetic_energy')
        dust = _dust('kinetic_energy')

        @snap.add_array(vector=vector, dust=dust)
        def kinetic_energy(snap) -> ndarray:
            """Kinetic energy."""
            return particles.kinetic_energy(snap=snap)

    if _add_array(snap, 'specific_kinetic_energy'):
        vector = _vector('specific_kinetic_energy')
        dust = _dust('specific_kinetic_energy')

        @snap.add_array(vector=vector, dust=dust)
        def specific_kinetic_energy(snap) -> ndarray:
            """Specific kinetic energy."""
            return particles.specific_kinetic_energy(snap=snap)

    if _add_array(snap, 'keplerian_frequency'):
        vector = _vector('keplerian_frequency')
        dust = _dust('keplerian_frequency')

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
        vector = _vector('semi_major_axis')
        dust = _dust('semi_major_axis')

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
        vector = _vector('eccentricity')
        dust = _dust('eccentricity')

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
        vector = _vector('inclination')
        dust = _dust('inclination')

        @snap.add_array(vector=vector, dust=dust)
        def inclination(snap) -> ndarray:
            """Inclination."""
            return particles.inclination(snap=snap)

    if _add_array(snap, 'radial_distance'):
        vector = _vector('radial_distance')
        dust = _dust('radial_distance')

        @snap.add_array(vector=vector, dust=dust)
        def radius_cylindrical(snap) -> ndarray:
            """Cylindrical radius."""
            return particles.radial_distance(snap=snap, coordinates='cylindrical')

        @snap.add_array(vector=vector, dust=dust)
        def radius_spherical(snap) -> ndarray:
            """Spherical radius."""
            return particles.radial_distance(snap=snap, coordinates='spherical')

    if _add_array(snap, 'azimuthal_angle'):
        vector = _vector('azimuthal_angle')
        dust = _dust('azimuthal_angle')

        @snap.add_array(vector=vector, dust=dust)
        def azimuthal_angle(snap) -> ndarray:
            """Azimuthal angle."""
            return particles.azimuthal_angle(snap=snap)

    if _add_array(snap, 'polar_angle'):
        vector = _vector('polar_angle')
        dust = _dust('polar_angle')

        @snap.add_array(vector=vector, dust=dust)
        def polar_angle(snap) -> ndarray:
            """Polar angle."""
            return particles.polar_angle(snap=snap)

    if _add_array(snap, 'radial_velocity'):
        vector = _vector('radial_velocity')
        dust = _dust('radial_velocity')

        @snap.add_array(vector=vector, dust=dust)
        def velocity_radial_cylindrical(snap) -> ndarray:
            """Cylindrical radial velocity."""
            return particles.radial_velocity(snap=snap, coordinates='cylindrical')

        @snap.add_array(vector=vector, dust=dust)
        def velocity_radial_spherical(snap) -> ndarray:
            """Spherical radial velocity."""
            return particles.radial_velocity(snap=snap, coordinates='spherical')

    if _add_array(snap, 'angular_velocity'):
        vector = _vector('angular_velocity')
        dust = _dust('angular_velocity')

        @snap.add_array(vector=vector, dust=dust)
        def angular_velocity(snap) -> ndarray:
            """Angular velocity."""
            return particles.angular_velocity(snap=snap)

    if _add_array(snap, 'temperature'):
        vector = _vector('temperature')
        dust = _dust('temperature')

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
        vector = _vector('gas_fraction')
        dust = _dust('gas_fraction')

        @snap.add_array(vector=vector, dust=dust)
        def gas_fraction(snap) -> ndarray:
            """Gas fraction."""
            return particles.gas_fraction(snap=snap)

    if _add_array(snap, 'gas_mass'):
        vector = _vector('gas_mass')
        dust = _dust('gas_mass')

        @snap.add_array(vector=vector, dust=dust)
        def gas_mass(snap) -> ndarray:
            """Gas mass."""
            return particles.gas_mass(snap=snap)

    if _add_array(snap, 'gas_density'):
        vector = _vector('gas_density')
        dust = _dust('gas_density')

        @snap.add_array(vector=vector, dust=dust)
        def gas_density(snap) -> ndarray:
            """Gas density."""
            return particles.gas_density(snap=snap)

    if _add_array(snap, 'dust_fraction'):
        vector = _vector('dust_fraction')
        dust = _dust('dust_fraction')

        if snap.properties['dust_method'] == 'dust as separate sets of particles':

            @snap.add_array(vector=vector, dust=dust)
            def dust_fraction(snap) -> ndarray:
                """Dust fraction."""
                return particles.dust_fraction(snap=snap)

    if _add_array(snap, 'dust_mass'):
        vector = _vector('dust_mass')
        dust = _dust('dust_mass')

        @snap.add_array(vector=vector, dust=dust)
        def dust_mass(snap) -> ndarray:
            """Dust mass."""
            return particles.dust_mass(snap=snap)

    if _add_array(snap, 'dust_density'):
        vector = _vector('dust_density')
        dust = _dust('dust_density')

        @snap.add_array(vector=vector, dust=dust)
        def dust_density(snap) -> ndarray:
            """Dust density."""
            return particles.dust_density(snap=snap)

    if _add_array(snap, 'stokes_number'):
        vector = _vector('stokes_number')
        dust = _dust('stokes_number')

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
