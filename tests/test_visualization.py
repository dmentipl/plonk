"""
Testing Visualization.
"""

import contextlib
import io
import pathlib
import unittest

import numpy as np

import plonk

I_GAS = 1
I_DUST = 7

TEST_FILE = pathlib.Path(__file__).parent / 'stubdata/phantom_00000.h5'

TEXT_TRAP = io.StringIO()


class TestInitializeVisualization(unittest.TestCase):
    """Test initialization of Visualization object."""

    def test_initialization(self):

        dump = plonk.Dump(TEST_FILE)

        with contextlib.redirect_stdout(TEXT_TRAP):
            plonk.Visualization(dump)
            plonk.Visualization(dump, render='density')
            plonk.Visualization(dump, vector='velocity')


class TestNumpy(unittest.TestCase):
    """Test rendering with NumPy array."""

    def test_initialization(self):

        dump = plonk.Dump(TEST_FILE)

        with contextlib.redirect_stdout(TEXT_TRAP):
            plonk.Visualization(dump, render=np.linspace(0, 1))


class TestExtraQuantity(unittest.TestCase):
    """Test rendering with extra quantity."""

    def test_extra_quantity(self):

        dump = plonk.Dump(TEST_FILE)

        with contextlib.redirect_stdout(TEXT_TRAP):
            plonk.Visualization(dump, render='r')


class TestSymbolic(unittest.TestCase):
    """Test rendering a symbolic expression."""

    def test_symbolic(self):

        dump = plonk.Dump(TEST_FILE)

        with contextlib.redirect_stdout(TEXT_TRAP):
            plonk.Visualization(dump, render='x**2 + y**2')


class TestUnits(unittest.TestCase):
    """Test using physical units."""

    def test_units(self):

        dump = plonk.Dump(TEST_FILE)
        units = plonk.Units(length='cm', mass='g', time='s')

        with contextlib.redirect_stdout(TEXT_TRAP):
            viz = plonk.Visualization(dump, render='density', units=units)

            viz = plonk.Visualization(dump, render='density')
            viz.set_units(units, integrated_z=1.0)


class TestRenderScale(unittest.TestCase):
    """Test changing the render scale."""

    def test_set_render_scale(self):

        dump = plonk.Dump(TEST_FILE)

        with contextlib.redirect_stdout(TEXT_TRAP):
            viz = plonk.Visualization(dump, render='density')
            viz.set_render_scale('linear')

        self.assertEqual(viz._render_scale, 'linear')


class TestRenderRange(unittest.TestCase):
    """Test changing the render range."""

    def test_set_render_range(self):

        dump = plonk.Dump(TEST_FILE)

        with contextlib.redirect_stdout(TEXT_TRAP):
            viz = plonk.Visualization(dump, render='density')
            viz.set_render_range(vmin=0.0, vmax=1.0)

        self.assertEqual(viz._options.render.render_min, 0.0)
        self.assertEqual(viz._options.render.render_max, 1.0)


class TestColormap(unittest.TestCase):
    """Test changing the colormap."""

    def test_set_colormap(self):

        dump = plonk.Dump(TEST_FILE)

        with contextlib.redirect_stdout(TEXT_TRAP):
            viz = plonk.Visualization(dump, render='density')
            viz.set_colormap('viridis')

        self.assertEqual(viz._options.figure.colormap, 'viridis')


class TestParticleType(unittest.TestCase):
    """Test changing the particle type."""

    def test_set_particle_type(self):

        dump = plonk.Dump(TEST_FILE)

        with contextlib.redirect_stdout(TEXT_TRAP):
            viz = plonk.Visualization(dump, render='density')
            viz.set_particle_type(I_DUST)

        self.assertEqual(viz._particle_types, {I_DUST})

        with contextlib.redirect_stdout(TEXT_TRAP):
            viz.set_particle_type([I_GAS, I_DUST])

        self.assertEqual(viz._particle_types, {I_GAS, I_DUST})


class TestRotate(unittest.TestCase):
    """Test rotating the frame."""

    def test_rotate_frame(self):

        dump = plonk.Dump(TEST_FILE)

        with contextlib.redirect_stdout(TEXT_TRAP):
            viz = plonk.Visualization(dump, render='density')
            viz.rotate_frame([1, 0, 0], np.pi / 2)

        np.testing.assert_array_almost_equal(
            viz._options.rotation.rotation_axis, [1, 0, 0]
        )
        self.assertEqual(viz._options.rotation.rotation_angle, np.pi / 2)

        with contextlib.redirect_stdout(TEXT_TRAP):
            viz = plonk.Visualization(
                dump,
                render='density',
                rotation_axis=[1, 0, 0],
                rotation_angle=np.pi / 2,
            )

        np.testing.assert_array_almost_equal(
            viz._options.rotation.rotation_axis, [1, 0, 0]
        )
        self.assertEqual(viz._options.rotation.rotation_angle, np.pi / 2)

        with contextlib.redirect_stdout(TEXT_TRAP):
            viz = plonk.Visualization(
                dump, render='density', position_angle=np.pi / 2, inclination=np.pi / 2
            )

        np.testing.assert_array_almost_equal(
            viz._options.rotation.rotation_axis, [0, 1, 0]
        )
        self.assertEqual(viz._options.rotation.rotation_angle, np.pi / 2)


class TestImageSize(unittest.TestCase):
    """Test setting the image size."""

    def test_image_size(self):

        dump = plonk.Dump(TEST_FILE)

        with contextlib.redirect_stdout(TEXT_TRAP):
            viz = plonk.Visualization(dump, render='density')
            viz.set_image_size(extent=(-100, 100, -100, 100))

        self.assertEqual(viz._extent, (-100, 100, -100, 100))

        with contextlib.redirect_stdout(TEXT_TRAP):
            viz.set_image_size(extent=(-100, 100, -100, 100))
            viz.set_image_size(size=200)

        self.assertEqual(viz._extent, (-200, 200, -200, 200))

        with contextlib.redirect_stdout(TEXT_TRAP):
            viz.set_image_size(size=200)


if __name__ == '__main__':
    unittest.main(verbosity=2)
