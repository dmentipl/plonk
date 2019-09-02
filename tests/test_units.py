"""
Testing units.
"""

import unittest

import plonk


class TestInitializeUnits(unittest.TestCase):
    """Test initialize Units."""

    def test_initialize_units(self):

        plonk.Units()
        plonk.Units(1.0, 1.0, 1.0)
        plonk.Units('cm', 'g', 's')
        plonk.Units(momentum=1.0, E=1.0)
        self.assertRaises(ValueError, plonk.Units, non_existent_unit=1.0)
        self.assertRaises(ValueError, plonk.Units, length='non_existent_length')


class TestUnitConversion(unittest.TestCase):
    """Test converting units."""

    def test_unit_conversion(self):

        units = plonk.Units('au', 'g', 's')
        au = units.convert_quantity_to_cgs(1.0, 'L')
        self.assertEqual(au, 1.496e13)

        units = plonk.Units('au', 'g', 'yr')
        au_per_year = units.convert_quantity_to_cgs(1.0, 'L T^-1')
        self.assertEqual(au_per_year, 474063.91864657536)


class TestDimensionComparison(unittest.TestCase):
    """Test comparing dimensions."""

    def test_dimension_comparison(self):

        self.assertTrue(plonk.core.units.is_dimension_same('L T^-1', 'T^-1 L'))
        self.assertTrue(plonk.core.units.is_dimension_same('L T^-1', 'velocity'))
        self.assertTrue(plonk.core.units.is_dimension_same('L T^-1', 'v'))
        self.assertTrue(plonk.core.units.is_dimension_same('velocity', 'L T^-1'))
        self.assertTrue(plonk.core.units.is_dimension_same('v', 'L T^-1'))


if __name__ == '__main__':
    unittest.main(verbosity=2)
