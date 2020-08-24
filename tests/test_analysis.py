"""Test analysis functions."""

from pathlib import Path

import pytest

import plonk
from plonk import analysis

from .data.phantom import adiabatic, dustmixture, dustseparate, mhd

SNAPTYPES = [adiabatic, dustmixture, dustseparate, mhd]
DIR = Path(__file__).parent / 'data/phantom'
AU = plonk.units('au')


@pytest.mark.parametrize('snaptype', SNAPTYPES)
def test_particles(snaptype):
    """Test particles analysis functions."""
    filename = DIR / snaptype.filename
    snap = plonk.load_snap(filename)

    snap.set_molecular_weight(2.381)
    if snap.num_sinks > 0:
        snap.set_gravitational_parameter(0)
    mu = snap.properties.get('gravitational_parameter')

    _test_particles(snap=snap, ignore=False, mu=mu)
    _test_particles(snap=snap, ignore=True, mu=mu)

    snap.close_file()


def _test_particles(snap, ignore, mu):

    analysis.particles.angular_momentum(snap=snap, ignore_accreted=ignore)
    analysis.particles.angular_velocity(snap=snap, ignore_accreted=ignore)
    analysis.particles.azimuthal_angle(snap=snap, ignore_accreted=ignore)
    analysis.particles.inclination(snap=snap, ignore_accreted=ignore)
    analysis.particles.kinetic_energy(snap=snap, ignore_accreted=ignore)
    analysis.particles.momentum(snap=snap, ignore_accreted=ignore)
    analysis.particles.polar_angle(snap=snap, ignore_accreted=ignore)
    analysis.particles.radius_cylindrical(snap=snap, ignore_accreted=ignore)
    analysis.particles.radius_spherical(snap=snap, ignore_accreted=ignore)
    analysis.particles.specific_angular_momentum(snap=snap, ignore_accreted=ignore)
    analysis.particles.specific_kinetic_energy(snap=snap, ignore_accreted=ignore)
    analysis.particles.temperature(snap=snap, ignore_accreted=ignore)

    if snap.num_dust_species > 0:
        if snap.properties['dust_method'] == 'dust/gas mixture':
            analysis.particles.dust_density(snap=snap, ignore_accreted=ignore)
            analysis.particles.dust_mass(snap=snap, ignore_accreted=ignore)
            analysis.particles.gas_density(snap=snap, ignore_accreted=ignore)
            analysis.particles.gas_fraction(snap=snap, ignore_accreted=ignore)
            analysis.particles.gas_mass(snap=snap, ignore_accreted=ignore)

    if mu is not None:
        analysis.particles.eccentricity(
            snap=snap, gravitational_parameter=mu, ignore_accreted=ignore
        )
        analysis.particles.keplerian_frequency(
            snap=snap, gravitational_parameter=mu, ignore_accreted=ignore
        )
        analysis.particles.semi_major_axis(
            snap=snap, gravitational_parameter=mu, ignore_accreted=ignore
        )
        analysis.particles.stokes_number(
            snap=snap, gravitational_parameter=mu, ignore_accreted=ignore
        )


@pytest.mark.parametrize('snaptype', SNAPTYPES)
def test_total(snaptype):
    """Test total analysis functions."""
    filename = DIR / snaptype.filename
    snap = plonk.load_snap(filename)

    analysis.total.accreted_mass(snap=snap)
    analysis.total.angular_momentum(snap=snap)
    analysis.total.center_of_mass(snap=snap)
    analysis.total.inclination(snap=snap)
    analysis.total.kinetic_energy(snap=snap)
    analysis.total.mass(snap=snap)
    analysis.total.momentum(snap=snap)
    analysis.total.position_angle(snap=snap)
    analysis.total.specific_angular_momentum(snap=snap)
    analysis.total.specific_kinetic_energy(snap=snap)

    snap.close_file()


@pytest.mark.parametrize('snaptype', SNAPTYPES)
def test_filters(snaptype):
    """Test particle filter functions."""
    filename = DIR / snaptype.filename
    snap = plonk.load_snap(filename)

    xwidth, ywidth, zwidth = 10 * AU, 10 * AU, 10 * AU
    height = 10 * AU
    radius = 100 * AU
    radius_min = 10 * AU
    radius_max = 20 * AU

    analysis.filters.annulus(
        snap=snap, radius_min=radius_min, radius_max=radius_max, height=height
    )
    analysis.filters.box(snap=snap, xwidth=xwidth, ywidth=ywidth, zwidth=zwidth)
    analysis.filters.cylinder(snap=snap, radius=radius, height=height)
    analysis.filters.shell(snap=snap, radius_min=radius_min, radius_max=radius_max)
    analysis.filters.sphere(snap=snap, radius=radius)

    snap.close_file()


@pytest.mark.parametrize('snaptype', [dustmixture, dustseparate])
def test_discs(snaptype):
    """Test discs analysis functions."""
    filename = DIR / snaptype.filename
    snap = plonk.load_snap(filename)

    analysis.discs.normal(snap)
    analysis.discs.rotate_edge_on(snap)
    analysis.discs.rotate_face_on(snap)

    snap.close_file()
