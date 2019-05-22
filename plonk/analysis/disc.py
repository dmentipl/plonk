import warnings

import numpy as np
import pandas as pd

from ..core.constants import constants
from ..core.particles import (
    _angular_momentum,
    _cylindrical_radius,
    _eccentricity,
)

_MIN_PARTICLES_IN_BIN = 5


def disc_analysis(
    dump,
    radius_in=None,
    radius_out=None,
    number_radial_bins=None,
    sink_index=None,
):
    """
    Disc analysis.

    Disc analysis takes a dump file with a single gas disc around a sink
    and computes azimuthal and vertical averages in radial bins.

    Parameters
    ----------
    dump : Dump
        Dump object.

    radius_in : float
        Inner disc radius for radial binning.

    radius_out : float
        Outer disc radius for radial binning.

    number_radial_bins : int (default 200)
        Number of radial bins.

    sink_index : int (default 1)
        Index (counting from 1) of the sink representing the central
        star.

    Returns
    -------
    averages : pandas.DataFrame
        The DataFrame contains the following quantities averaged over
        the disc:
            cyclindrical radius ['R']
            number of particles ['npart']
            surface density ['sigma']
            midplane density ['rho']
            scale height ['H']
            average smoothing length ['h']
            angular momentum ['Lx', 'Ly', 'Lz', '|L|']
            eccentricity ['e']
            tilt ['tilt']
            twist ['twist']
            psi ['psi']
    """

    if radius_in is None:
        raise ValueError('Need radius_in')

    if radius_out is None:
        raise ValueError('Need radius_out')

    if number_radial_bins is None:
        number_radial_bins = 200

    if sink_index is None:
        sink_index = 1

    if np.unique(dump.particles.arrays['itype']).size > 1:
        warnings.warn('Disc analysis uses all particle itypes.')

    return _calculate_radially_binned_quantities(
        dump, number_radial_bins, radius_in, radius_out, sink_index
    )


def _calculate_radially_binned_quantities(
    dump, number_radial_bins, radius_in, radius_out, sink_index
):

    particles = dump.particles.to_structured_array()
    particle_masses = dump.particles.mass

    radius = _cylindrical_radius(particles['xyz'])

    radial_bin_width = (radius_out - radius_in) / (number_radial_bins - 1)
    radial_bins = np.linspace(radius_in, radius_out, number_radial_bins)

    gravitational_constant = constants.gravitational_constant / (
        dump.header['udist'] ** 3
        / dump.header['umass']
        / dump.header['utime'] ** 2
    )
    stellar_mass = dump.sinks.arrays['m'][sink_index - 1]
    gravitational_parameter = gravitational_constant * stellar_mass
    eccentricity = _eccentricity(
        particles['xyz'], particles['vxyz'], gravitational_parameter
    )

    averages = pd.DataFrame(radial_bins, columns=['R'])

    surface_area = np.pi * (
        (radial_bins + radial_bin_width / 2) ** 2
        - (radial_bins - radial_bin_width / 2) ** 2
    )

    averages = averages.reindex(
        columns=averages.columns.tolist()
        + [
            'npart',
            'sigma',
            'h',
            'H',
            'Lx',
            'Ly',
            'Lz',
            '|L|',
            'tilt',
            'twist',
            'e',
        ]
    )

    radius_left = radial_bins - radial_bin_width / 2
    radius_right = radial_bins + radial_bin_width / 2

    masks = [
        (radius >= radius_left) & (radius <= radius_right)
        for radius_left, radius_right in zip(radius_left, radius_right)
    ]

    for index, mask in enumerate(masks):

        part = particles[mask]
        pmass = particle_masses[mask]
        ecc = eccentricity[mask]

        n_particles_in_bin = len(part)
        averages['npart'].iloc[index] = n_particles_in_bin

        if n_particles_in_bin < _MIN_PARTICLES_IN_BIN:
            continue

        averages['sigma'].iloc[index] = (
            np.sum(n_particles_in_bin * pmass) / surface_area[index]
        )

        averages['h'].iloc[index] = part['h'].mean()
        averages['H'].iloc[index] = part['xyz'][:, 2].std()

        angular_momentum = _angular_momentum(
            part['xyz'], part['vxyz'], pmass
        ).mean(axis=0)

        averages['Lx'].iloc[index] = angular_momentum[0]
        averages['Ly'].iloc[index] = angular_momentum[1]
        averages['Lz'].iloc[index] = angular_momentum[2]
        averages['|L|'].iloc[index] = np.linalg.norm(angular_momentum)

        averages['e'].iloc[index] = ecc.mean()

    lx = averages['Lx'] / averages['|L|']
    ly = averages['Ly'] / averages['|L|']
    lz = averages['Lz'] / averages['|L|']

    averages['tilt'] = np.arccos(lz)
    averages['twist'] = np.arctan2(ly, lx)
    averages['rho'] = averages['sigma'] / averages['H'] / np.sqrt(2 * np.pi)
    averages['psi'] = radial_bins * np.sqrt(
        np.gradient(lx, radial_bin_width) ** 2
        + np.gradient(ly, radial_bin_width) ** 2
        + np.gradient(lz, radial_bin_width) ** 2
    )

    return averages
