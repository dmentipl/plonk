"""
disc.py

Analysis for dusty discs.

Daniel Mentiplay, 2019.
"""

import numpy as np
from numpy.linalg import norm
import pandas as pd

from ..constants import constants


# ---------------------------------------------------------------------------- #

def disc_analysis(radius_in=None,
                  radius_out=None,
                  number_radial_bins=None,
                  dump=None,
                  min_particle_average=None):
    """
    Perform disc analysis.

    TODO: add more docs

    Arguments:
        radius_in          : Inner disc radius for radial binning
        radius_out         : Outer disc radius for radial binning
        number_radial_bins : Number of radial bins
        dump               : Dump object

    Optional:
        min_particle_average : Minimum number of particles to compute averages
    """

    # TODO: add to docs

# --- Calculate extra quantities on particles

    _calculate_extra_quantities(dump)

# --- Calculate radially binned quantities

    radial_data = _calculate_radially_binned_quantities(number_radial_bins,
                                                        radius_in,
                                                        radius_out,
                                                        dump.particles,
                                                        min_particle_average)

# --- Dust

    n_dust_small = dump.parameters['ndustsmall']
    n_dust_large = dump.parameters['ndustlarge']

    contains_small_dust = bool(n_dust_small > 0)
    contains_large_dust = bool(n_dust_large > 0)
    contains_dust = bool(contains_small_dust or contains_large_dust)

    if contains_dust:
        grain_size = dump.parameters['grainsize']
        grain_dens = dump.parameters['graindens']

    gamma = dump.parameters['gamma']

# --- Stokes

    # for idxi in range(len(radial_bins_disc)):
    #     for idxj in range(n_dust_large):

    #         Stokes[idxj][idxi] = \
    #             np.sqrt(gamma*np.pi/8) * grain_dens[idxj] * grain_size[idxj] \
    #             / (scale_height_gas[idxi] \
    #             * (midplane_density_gas[idxi] + midplane_density_dust[idxj][idxi]))

# --- Return

    return radial_data


# ---------------------------------------------------------------------------- #

def _calculate_extra_quantities(dump):
    """
    Calculate extra quantities on particles and sinks appropriate for a disc.
    """

    print('Calculating extra quantities on the particles and sinks...')
    print('And adding them to the DataFrame\n')
    particles = dump.particles
    sinks = dump.sinks
    parameters = dump.parameters
    units = dump.units

# --- Dump type

    is_full_dump = 'vx' in particles.columns

# --- Units

    u_dist = units['distance']
    u_time = units['time']
    u_mass = units['mass']

# --- Spherical and cylindrical distance

    particles['r'] = norm(particles[['x', 'y', 'z']], axis=1)
    particles['R'] = norm(particles[['x', 'y']], axis=1)

# --- Momentum and angular momentum

    if is_full_dump:

        particles['|v|'] = norm(particles[['vx', 'vy', 'vz']], axis=1)

        particles['px'] = particles['m'] * particles['vx']
        particles['py'] = particles['m'] * particles['vy']
        particles['pz'] = particles['m'] * particles['vz']

        angular_momentum = np.cross(particles[['x', 'y', 'z']],
                                    particles[['px', 'py', 'pz']])

        particles['Lx'] = angular_momentum[:, 0]
        particles['Ly'] = angular_momentum[:, 1]
        particles['Lz'] = angular_momentum[:, 2]
        particles['|L|'] = norm(angular_momentum, axis=1)

    else:

        particles['|v|'] = np.nan
        particles['px'] = np.nan
        particles['py'] = np.nan
        particles['pz'] = np.nan
        particles['Lx'] = np.nan
        particles['Ly'] = np.nan
        particles['Lz'] = np.nan
        particles['|L|'] = np.nan

# --- Sinks

    n_sinks = parameters['nptmass']
    contains_sinks = bool(n_sinks > 0)

    if contains_sinks:

        sinks['r'] = norm(sinks[['x', 'y', 'z']], axis=1)
        sinks['R'] = norm(sinks[['x', 'y']], axis=1)

        # TODO: check if sink[0] is really the star; check if binary
        stellar_mass = sinks['m'][0]
        print('Assuming the first sink particle is the central star')

        gravitational_parameter = \
            constants.gravitational_constant \
            / (u_dist**3 / u_time**2 / u_mass) * stellar_mass
    else:

        raise Exception('Must have at least one sink to do disc analysis')

    if is_full_dump:

        sinks['|v|'] = norm(sinks[['vx', 'vy', 'vz']], axis=1)

        sinks['px'] = sinks['m'] * sinks['vx']
        sinks['py'] = sinks['m'] * sinks['vy']
        sinks['pz'] = sinks['m'] * sinks['vz']

        angular_momentum = np.cross(sinks[['x', 'y', 'z']],
                                    sinks[['px', 'py', 'pz']])

        sinks['Lx'] = angular_momentum[:, 0]
        sinks['Ly'] = angular_momentum[:, 1]
        sinks['Lz'] = angular_momentum[:, 2]
        sinks['|L|'] = norm(angular_momentum, axis=1)

    else:

        sinks['|v|'] = np.nan
        sinks['px'] = np.nan
        sinks['py'] = np.nan
        sinks['pz'] = np.nan
        sinks['Lx'] = np.nan
        sinks['Ly'] = np.nan
        sinks['Lz'] = np.nan
        sinks['|L|'] = np.nan

# --- Eccentricity

    if is_full_dump:

        kinetic_energy = 1/2 * particles['|v|']**2
        gravitational_energy = - gravitational_parameter / particles['r']

        energy = kinetic_energy + gravitational_energy
        term = 2 * energy * (particles['|L|']/particles['m'])**2 \
            / gravitational_parameter**2

        particles['e'] = np.sqrt(1 + term)

        kinetic_energy = 1/2 * sinks['|v|']**2
        gravitational_energy = - gravitational_parameter / sinks['r']

        energy = kinetic_energy + gravitational_energy
        term = 2 * energy * (sinks['|L|']/sinks['m'])**2 \
            / gravitational_parameter**2

        sinks['e'] = np.sqrt(1 + term)

    else:

        particles['e'] = np.nan
        sinks['e'] = np.nan


# ---------------------------------------------------------------------------- #

def _calculate_radially_binned_quantities(number_radial_bins=None,
                                          radius_in=None,
                                          radius_out=None,
                                          particles=None,
                                          min_particle_average=None):
    """
    Calculate averaged radially binned quantities:
        - radial bins
        - surface density
        - midplane density
        - scale height
        - smoothing length
        - angular momentum
        - tilt
        - twist
        - psi
        - eccentricity
    """

    if number_radial_bins is None:
        raise ValueError('Need number_radial_bins')

    if radius_in is None:
        raise ValueError('Need radius_in')

    if radius_out is None:
        raise ValueError('Need radius_out')

    if particles is None:
        raise ValueError('Need particles')

    if min_particle_average is None:
        min_particle_average = 5

    radial_bin_width = (radius_out - radius_in) / (number_radial_bins - 1)
    radial_bins = np.linspace(radius_in, radius_out, number_radial_bins)

    radial_averages = pd.DataFrame(radial_bins, columns=['R'])
    radial_averages['area'] = \
        np.pi * ((radial_bins + radial_bin_width/2)**2
                 (radial_bins - radial_bin_width/2)**2)

    radial_averages = radial_averages.reindex(
        columns=radial_averages.columns.tolist()
        + ['sigma', 'h', 'H', 'Lx', 'Ly', 'Lz', '|L|', 'tilt', 'twist', 'e'])

    for index, radius in enumerate(radial_bins):

        radius_left = radius - radial_bin_width/2
        radius_right = radius + radial_bin_width/2

        particles_in_bin = \
            particles.loc[(particles['R'] >= radius_left)
                          & (particles['R'] <= radius_right)]

        radial_averages['sigma'].iloc[index] = \
            particles_in_bin['m'].sum() / radial_averages['area'].iloc[index]

        if len(particles_in_bin) > min_particle_average:

            radial_averages['h'].iloc[index] = particles_in_bin['h'].mean()
            radial_averages['H'].iloc[index] = particles_in_bin['z'].std()

            radial_averages['Ly'].iloc[index] = particles_in_bin['Ly'].mean()
            radial_averages['Lx'].iloc[index] = particles_in_bin['Lx'].mean()
            radial_averages['Lz'].iloc[index] = particles_in_bin['Lz'].mean()
            radial_averages['|L|'].iloc[index] = particles_in_bin['|L|'].mean()

            radial_averages['e'].iloc[index] = particles_in_bin['e'].mean()

    radial_averages['tilt'] = \
        np.arccos(radial_averages['Lz'] / radial_averages['|L|'])

    radial_averages['twist'] = \
        np.arctan2(radial_averages['Ly'] / radial_averages['|L|'],
                   radial_averages['Lx'] / radial_averages['|L|'])

    radial_averages['rho'] = \
        radial_averages['sigma'] / radial_averages['H'] / np.sqrt(2*np.pi)

    lx = radial_averages['Lx'] / radial_averages['|L|']
    ly = radial_averages['Ly'] / radial_averages['|L|']
    lz = radial_averages['Lz'] / radial_averages['|L|']

    radial_averages['psi'] = radial_bins * np.sqrt(
        np.gradient(lx, radial_bin_width)**2 +
        np.gradient(ly, radial_bin_width)**2 +
        np.gradient(lz, radial_bin_width)**2)

    return radial_averages
