"""Test analysis functions."""

from pathlib import Path

import pytest

import plonk
from plonk.analysis import discs, filters, particles, total

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

    _test_particles(snap=snap, ignore=False)
    _test_particles(snap=snap, ignore=True)

    snap.close_file()


def _test_particles(snap, ignore):

    particles.angular_momentum(snap=snap, ignore_accreted=ignore)
    particles.angular_velocity(snap=snap, ignore_accreted=ignore)
    particles.azimuthal_angle(snap=snap, ignore_accreted=ignore)
    particles.kinetic_energy(snap=snap, ignore_accreted=ignore)
    particles.momentum(snap=snap, ignore_accreted=ignore)
    particles.polar_angle(snap=snap, ignore_accreted=ignore)
    particles.radius_cylindrical(snap=snap, ignore_accreted=ignore)
    particles.radius_spherical(snap=snap, ignore_accreted=ignore)
    particles.specific_angular_momentum(snap=snap, ignore_accreted=ignore)
    particles.specific_kinetic_energy(snap=snap, ignore_accreted=ignore)
    particles.temperature(snap=snap, ignore_accreted=ignore)

    if snap.num_dust_species > 0:
        if snap.properties['dust_method'] == 'dust/gas mixture':
            particles.dust_density(snap=snap, ignore_accreted=ignore)
            particles.dust_mass(snap=snap, ignore_accreted=ignore)
            particles.gas_density(snap=snap, ignore_accreted=ignore)
            particles.gas_fraction(snap=snap, ignore_accreted=ignore)
            particles.gas_mass(snap=snap, ignore_accreted=ignore)


@pytest.mark.parametrize('snaptype', SNAPTYPES)
def test_total(snaptype):
    """Test total analysis functions."""
    filename = DIR / snaptype.filename
    snap = plonk.load_snap(filename)

    total.accreted_mass(snap=snap)
    total.angular_momentum(snap=snap)
    total.center_of_mass(snap=snap)
    total.kinetic_energy(snap=snap)
    total.mass(snap=snap)
    total.momentum(snap=snap)
    total.specific_angular_momentum(snap=snap)
    total.specific_kinetic_energy(snap=snap)

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

    filters.annulus(
        snap=snap, radius_min=radius_min, radius_max=radius_max, height=height
    )
    filters.box(snap=snap, xwidth=xwidth, ywidth=ywidth, zwidth=zwidth)
    filters.cylinder(snap=snap, radius=radius, height=height)
    filters.shell(snap=snap, radius_min=radius_min, radius_max=radius_max)
    filters.sphere(snap=snap, radius=radius)

    snap.close_file()


@pytest.mark.parametrize('snaptype', [dustmixture, dustseparate])
def test_discs(snaptype):
    """Test discs analysis functions."""
    filename = DIR / snaptype.filename
    snap = plonk.load_snap(filename)

    if snap.num_sinks > 0:
        snap.set_central_body(0)

    discs.unit_normal(snap=snap)
    discs.rotate_edge_on(snap=snap)
    discs.rotate_face_on(snap=snap)

    discs.position_angle(snap=snap)
    discs.inclination_angle(snap=snap)

    discs.eccentricity(snap=snap)
    discs.inclination(snap=snap)
    discs.keplerian_frequency(snap=snap)
    discs.semi_major_axis(snap=snap)
    discs.stokes_number(snap=snap)

    snap.close_file()
