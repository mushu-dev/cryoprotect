"""
Tests for the property transformer module.
"""

import unittest
from typing import Dict, Any, List

from ..transforms.property_transform import PropertyTransformer


class TestPropertyTransformer(unittest.TestCase):
    """Tests for the PropertyTransformer class."""

    def test_standardize_property_name(self):
        """Test standardization of property names."""
        transformer = PropertyTransformer()
        
        # Test various property names
        self.assertEqual(transformer.standardize_property_name('logp'), 'LogP')
        self.assertEqual(transformer.standardize_property_name('xlogp3'), 'LogP')
        self.assertEqual(transformer.standardize_property_name('molecular_weight'), 'MolecularWeight')
        self.assertEqual(transformer.standardize_property_name('mol_weight'), 'MolecularWeight')
        self.assertEqual(transformer.standardize_property_name('topological_polar_surface_area'), 'TPSA')
        self.assertEqual(transformer.standardize_property_name('h_bond_donor_count'), 'HBondDonorCount')
        
        # Test non-standard name (should capitalize words)
        self.assertEqual(transformer.standardize_property_name('custom_property'), 'CustomProperty')

    def test_standardize_unit(self):
        """Test standardization of units."""
        transformer = PropertyTransformer()
        
        # Test various units
        self.assertEqual(transformer.standardize_unit('g/mol'), 'g/mol')
        self.assertEqual(transformer.standardize_unit('daltons'), 'g/mol')
        self.assertEqual(transformer.standardize_unit('Da'), 'g/mol')
        self.assertEqual(transformer.standardize_unit('a^2'), 'Å²')
        self.assertEqual(transformer.standardize_unit('angstrom^2'), 'Å²')
        self.assertEqual(transformer.standardize_unit('mg/ml'), 'mg/mL')
        
        # Test case insensitivity
        self.assertEqual(transformer.standardize_unit('G/MOL'), 'g/mol')
        
        # Test non-standard unit (should return as is)
        self.assertEqual(transformer.standardize_unit('custom_unit'), 'custom_unit')

    def test_convert_unit(self):
        """Test unit conversion."""
        transformer = PropertyTransformer()

        # Test temperature conversions
        self.assertAlmostEqual(transformer.convert_unit(32, '°F', '°C'), 0)
        self.assertAlmostEqual(transformer.convert_unit(100, '°C', '°F'), 212)
        self.assertAlmostEqual(transformer.convert_unit(0, '°C', 'K'), 273.15)

        # Test area conversions
        self.assertAlmostEqual(transformer.convert_unit(100, 'Å²', 'nm²'), 1)
        self.assertAlmostEqual(transformer.convert_unit(1, 'nm²', 'Å²'), 100)

        # Test molecular weight dependent conversion
        # mol/L to mg/mL conversion: 1 mol/L * 180.16 g/mol = 180.16 g/L = 180.16 mg/mL
        molecule_data = {'molecular_weight': 180.16}  # glucose
        result = transformer.convert_unit(1, 'mol/L', 'mg/mL', molecule_data)
        self.assertAlmostEqual(result, 180.16)

        # Test molecular weight conversions
        self.assertAlmostEqual(transformer.convert_unit(180.16, 'g/mol', 'kDa'), 0.18016)
        self.assertAlmostEqual(transformer.convert_unit(0.18016, 'kDa', 'g/mol'), 180.16)
        
        # Test unsupported conversion (should return None)
        self.assertIsNone(transformer.convert_unit(1, 'custom_unit', 'g/mol'))

    def test_standardize_property(self):
        """Test standardization of individual property dictionaries."""
        transformer = PropertyTransformer()
        
        # Test basic property standardization
        prop = {
            'property_name': 'molecular_weight',
            'numeric_value': 180.16,
            'unit': 'daltons'
        }
        
        result = transformer.standardize_property(prop)
        self.assertEqual(result['property_name'], 'MolecularWeight')
        self.assertEqual(result['unit'], 'g/mol')
        self.assertEqual(result['numeric_value'], 180.16)
        self.assertEqual(result['property_type'], 'physicochemical')
        self.assertEqual(result['source'], 'Unknown')
        
        # Test property with unit conversion
        prop = {
            'property_name': 'melting_point',
            'numeric_value': 32,
            'unit': 'F'
        }
        
        result = transformer.standardize_property(prop)
        self.assertEqual(result['property_name'], 'MeltingPoint')
        self.assertEqual(result['unit'], '°C')
        self.assertAlmostEqual(result['numeric_value'], 0)
        
        # Test property with missing unit
        prop = {
            'property_name': 'logp',
            'numeric_value': 2.5
        }
        
        result = transformer.standardize_property(prop)
        self.assertEqual(result['property_name'], 'LogP')
        self.assertEqual(result['unit'], 'log units')

    def test_config_options(self):
        """Test configuration options for PropertyTransformer."""
        # Test with standardize_units disabled
        config = {'standardize_units': False}
        transformer = PropertyTransformer(config=config)
        
        prop = {
            'property_name': 'molecular_weight',
            'numeric_value': 180.16,
            'unit': 'daltons'
        }
        
        result = transformer.standardize_property(prop)
        self.assertEqual(result['property_name'], 'MolecularWeight')
        self.assertEqual(result['unit'], 'daltons')  # Should not be standardized
        
        # Test with preferred units
        config = {
            'preferred_units': {
                'MolecularWeight': 'kDa'
            }
        }
        transformer = PropertyTransformer(config=config)
        
        prop = {
            'property_name': 'molecular_weight',
            'numeric_value': 180.16,
            'unit': 'g/mol'
        }
        
        # We now have a conversion function for g/mol to kDa in the updated implementation
        result = transformer.standardize_property(prop)
        self.assertEqual(result['unit'], 'kDa')
        self.assertAlmostEqual(result['numeric_value'], 0.18016)
        
        # Test with preserve_original_values option
        config = {
            'standardize_units': True,
            'preserve_original_values': True
        }
        transformer = PropertyTransformer(config=config)
        
        prop = {
            'property_name': 'melting_point',
            'numeric_value': 32,
            'unit': 'F'
        }
        
        result = transformer.standardize_property(prop)
        self.assertEqual(result['unit'], '°C')
        self.assertAlmostEqual(result['numeric_value'], 0)
        self.assertEqual(result['original_value'], 32)
        self.assertEqual(result['original_unit'], 'F')
        
        # Test with additional metadata
        config = {
            'additional_metadata': {
                'confidence': 'high',
                'method': 'calculated'
            }
        }
        transformer = PropertyTransformer(config=config)
        
        prop = {
            'property_name': 'logp',
            'numeric_value': 2.5
        }
        
        result = transformer.standardize_property(prop)
        self.assertEqual(result['confidence'], 'high')
        self.assertEqual(result['method'], 'calculated')

    def test_add_missing_properties(self):
        """Test adding missing properties from molecule data."""
        config = {'add_missing_properties': True}
        transformer = PropertyTransformer(config=config)
        
        # Create a molecule with properties
        molecule_data = {
            'name': 'Test Molecule',
            'properties': {
                'LogP': 2.5,
                'MolecularWeight': 180.16,
                'HBondDonorCount': 3
            }
        }
        
        # Some existing properties
        properties = [
            {
                'property_name': 'LogP',
                'numeric_value': 2.5,
                'source': 'Experimental'
            }
        ]
        
        # Standardize properties
        result = transformer.standardize_properties(properties, molecule_data)
        
        # Should have added the missing properties
        self.assertEqual(len(result), 3)  # LogP (existing) + MW and HBD (added)
        
        # Verify the added properties
        prop_names = [p['property_name'] for p in result]
        self.assertIn('LogP', prop_names)
        self.assertIn('MolecularWeight', prop_names)
        self.assertIn('HBondDonorCount', prop_names)
        
        # The original property should be preserved
        logp_prop = next(p for p in result if p['property_name'] == 'LogP')
        self.assertEqual(logp_prop['source'], 'Experimental')
        
        # The new properties should have the correct source
        mw_prop = next(p for p in result if p['property_name'] == 'MolecularWeight')
        self.assertEqual(mw_prop['source'], 'Derived')
        self.assertEqual(mw_prop['numeric_value'], 180.16)

    def test_merge_properties(self):
        """Test merging properties from different sources."""
        transformer = PropertyTransformer()
        
        # Primary properties
        primary = [
            {
                'property_name': 'LogP',
                'numeric_value': 2.5,
                'source': 'Experimental',
                'confidence_score': 0.9
            },
            {
                'property_name': 'MolecularWeight',
                'numeric_value': 180.16,
                'unit': 'g/mol',
                'source': 'Calculated'
            }
        ]
        
        # Secondary properties (some overlap, some new)
        secondary = [
            {
                'property_name': 'LogP',
                'numeric_value': 2.7,
                'source': 'Predicted',
                'confidence_score': 0.7
            },
            {
                'property_name': 'HBondDonorCount',
                'numeric_value': 3,
                'source': 'Calculated'
            }
        ]
        
        # Merge properties
        result = transformer.merge_properties(primary, secondary)
        
        # Should have all unique properties
        self.assertEqual(len(result), 3)
        
        # The LogP from primary should be preserved (higher confidence)
        logp_prop = next(p for p in result if p['property_name'] == 'LogP')
        self.assertEqual(logp_prop['numeric_value'], 2.5)
        self.assertEqual(logp_prop['source'], 'Experimental')
        
        # The new property should be added
        hbd_prop = next(p for p in result if p['property_name'] == 'HBondDonorCount')
        self.assertEqual(hbd_prop['numeric_value'], 3)
        
        # Test case where secondary has higher confidence
        primary[0]['confidence_score'] = 0.5  # Lower than secondary's 0.7
        
        result = transformer.merge_properties(primary, secondary)
        logp_prop = next(p for p in result if p['property_name'] == 'LogP')
        self.assertEqual(logp_prop['numeric_value'], 2.7)  # Should use secondary's value
        self.assertEqual(logp_prop['source'], 'Predicted')


if __name__ == '__main__':
    unittest.main()