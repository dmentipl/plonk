"""Extra pre-defined profiles."""

from functools import partial

from .._logging import logger
from . import profiles


def extra_profiles(profile, num_separate_dust: int = 0, num_mixture_dust: int = 0):
    """Make extra profiles available.

    Parameters
    ----------
    profile
        The profile object to add extra profiles to.
    num_separate_dust
        The number of "separate sets of particles" dust species.
    num_mixture_dust
        The number of "mixture" dust species.
    """
    num_dust_species = num_mixture_dust + num_separate_dust

    for prof, ndims in profiles.appropriate_ndim.items():
        if profile.ndim in ndims:
            if prof not in profiles.dust_profiles:
                profile._profile_functions[prof] = getattr(profiles, prof)
            else:
                if profiles.dust_profiles[prof] == 'both':
                    for idx in range(num_dust_species):
                        profile._profile_functions[f'{prof}_{idx+1:03}'] = partial(
                            getattr(profiles, prof), idx
                        )
                elif profiles.dust_profiles[prof] == 'mixture':
                    for idx in range(num_mixture_dust):
                        profile._profile_functions[f'{prof}_{idx+1:03}'] = partial(
                            getattr(profiles, prof), idx
                        )
                elif profiles.dust_profiles[prof] == 'mixture (gas)':
                    if num_mixture_dust > 0:
                        profile._profile_functions[prof] = getattr(profiles, prof)
                else:
                    logger.debug(f'{prof} - cannot add profile')
