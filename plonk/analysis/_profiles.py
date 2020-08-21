"""Extra pre-defined profiles."""

from functools import partial

import numpy as np

from .._units import Quantity
from .._units import units as plonk_units

G = (1 * plonk_units.newtonian_constant_of_gravitation).to_base_units()


def extra_profiles(profile, num_mixture_dust_species: int = 0):
    """Make extra profiles available.

    Parameters
    ----------
    profile
        The profile object to add extra profiles to.
    num_mixture_dust_species
        The number of "mixture" dust species.
    """

    @profile.add_profile
    def mass(prof) -> Quantity:
        """Mass profile."""
        M = prof.snap['mass']
        return prof.particles_to_binned_quantity('sum', M)

    @profile.add_profile
    def surface_density(prof) -> Quantity:
        """Surface density profile.

        Units are [mass / length ** ndim], which depends on ndim of profile.
        """
        return prof['mass'] / prof['size']

    @profile.add_profile
    def scale_height(prof) -> Quantity:
        """Scale height profile."""
        z = prof.snap['z']
        return prof.particles_to_binned_quantity('std', z)

    @profile.add_profile
    def aspect_ratio(prof) -> Quantity:
        """Aspect ratio profile."""
        H = prof['scale_height']
        R = prof['radius']
        return H / R

    @profile.add_profile
    def angular_momentum_theta(prof) -> Quantity:
        """Angle between specific angular momentum and xy-plane."""
        angular_momentum_z = prof['angular_momentum_z']
        angular_momentum_magnitude = prof['angular_momentum_mag']

        return np.arccos(angular_momentum_z / angular_momentum_magnitude)

    @profile.add_profile
    def angular_momentum_phi(prof) -> Quantity:
        """Angle between specific angular momentum and x-axis in xy-plane."""
        angular_momentum_x = prof['angular_momentum_x']
        angular_momentum_y = prof['angular_momentum_y']
        return np.arctan2(angular_momentum_y, angular_momentum_x)

    @profile.add_profile
    def toomre_Q(prof) -> Quantity:
        """Toomre Q parameter."""
        return (
            prof['sound_speed']
            * prof['keplerian_frequency']
            / (np.pi * G * prof['surface_density'])
        )

    if num_mixture_dust_species > 0:

        @profile.add_profile
        def gas_mass(prof) -> Quantity:
            """Gas mass profile."""
            M = prof.snap['gas_mass']
            return prof.particles_to_binned_quantity('sum', M)

        @profile.add_profile
        def gas_surface_density(prof) -> Quantity:
            """Gas surface density profile.

            Units are [mass / length ** ndim], which depends on ndim of profile.
            """
            return prof['gas_mass'] / prof['size']

        for idx in range(num_mixture_dust_species):

            def dust_mass(idx, prof) -> Quantity:
                """Dust mass profile."""
                M = prof.snap[f'dust_mass_{idx+1:03}']
                return prof.particles_to_binned_quantity('sum', M)

            def dust_surface_density(idx, prof) -> Quantity:
                """Dust surface density profile.

                Units are [mass / length ** ndim], which depends on ndim of profile.
                """
                return prof[f'dust_mass_{idx+1:03}'] / prof['size']

            profile._profile_functions[f'dust_mass_{idx+1:03}'] = partial(
                dust_mass, idx
            )

            profile._profile_functions[f'dust_surface_density_{idx+1:03}'] = partial(
                dust_surface_density, idx
            )
