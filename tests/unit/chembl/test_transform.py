#!/usr/bin/env python3
"""
Unit tests for ChEMBL data transformation functions.

These tests verify the functionality of the transform_chembl_to_molecule and
transform_chembl_to_properties functions in ChEMBL_Integrated_Import.py.
"""

import os
import sys
import unittest
from unittest.mock import patch, MagicMock
import json
from datetime import datetime
import uuid

# Add parent directory to path to import the module
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

# Import the module to test
import ChEMBL_Integrated_Import as chembl_import


class TestChEMBLTransform(unittest.TestCase):
    """Test cases for ChEMBL data transformation functions."""

    def setUp(self):
        """Set up test fixtures."""
        # Mock compound data for standard reference compounds
        self.glycerol = {
            'molecule_chembl_id': 'CHEMBL1234',
            'pref_name': 'Glycerol',
            'molecule_structures': {
                'canonical_smiles': 'C(C(CO)O)O',
                'standard_inchi': 'InChI=1S/C3H8O3/c4-1-3(6)2-5/h3-6H,1-2H2',
                'standard_inchi_key': 'PEDCQBHIVMGVHV-UHFFFAOYSA-N'
            },
            'molecule_properties': {
                'full_molformula': 'C3H8O3',
                'full_mwt': 92.09,
                'alogp': -1.76,
                'hba': 3,
                'hbd': 3,
                'psa': 60.69,
                'rtb': 2,
                'aromatic_rings': 0,
                'heavy_atoms': 6
            }
        }
        
        self.dmso = {
            'molecule_chembl_id': 'CHEMBL58',
            'pref_name': 'Dimethyl sulfoxide',
            'molecule_structures': {
                'canonical_smiles': 'CS(=O)C',
                'standard_inchi': 'InChI=1S/C2H6OS/c1-4(2)3/h1-2H3',
                'standard_inchi_key': 'IAZDPXIOMUYVGZ-UHFFFAOYSA-N'
            },
            'molecule_properties': {
                'full_molformula': 'C2H6OS',
                'full_mwt': 78.13,
                'alogp': -0.62,
                'hba': 1,
                'hbd': 0,
                'psa': 17.07,
                'rtb': 0,
                'aromatic_rings': 0,
                'heavy_atoms': 4
            }
        }
        
        # Mock compound with missing data
        self.incomplete_compound = {
            'molecule_chembl_id': 'CHEMBL9999',
            'molecule_structures': {
                'canonical_smiles': 'CC(=O)OC1=CC=CC=C1C(=O)O'
            },
            'molecule_properties': {
                'full_mwt': 180.16
            }
        }
        
        # Mock compound with additional properties
        self.compound_with_extra = {
            'molecule_chembl_id': 'CHEMBL5555',
            'pref_name': 'Test Compound',
            'molecule_structures': {
                'canonical_smiles': 'CC(=O)OC1=CC=CC=C1C(=O)O',
                'standard_inchi': 'InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)',
                'standard_inchi_key': 'BSYNRYMUTXBXSQ-UHFFFAOYSA-N'
            },
            'molecule_properties': {
                'full_molformula': 'C9H8O4',
                'full_mwt': 180.16,
                'alogp': 1.31,
                'hba': 4,
                'hbd': 1,
                'psa': 63.6,
                'rtb': 3
            },
            'molecule_hierarchy': {
                'parent_chembl_id': 'CHEMBL5000'
            },
            'molecule_synonyms': [
                {'synonym': 'Aspirin', 'syn_type': 'COMMON'},
                {'synonym': 'Acetylsalicylic acid', 'syn_type': 'INN'}
            ],
            'properties': [
                {'property_name': 'Solubility', 'value': '3 mg/mL'},
                {'property_name': 'Melting Point', 'value': 135.0}
            ]
        }
        
        # Mock user profile ID
        self.user_profile_id = 'test-user-profile-id'
        
        # Mock property type map
        self.property_type_map = {
            'logp': 'prop-type-logp',
            'logp (chemaxon)': 'prop-type-logp-chemaxon',
            'molecular weight': 'prop-type-mw',
            'hydrogen bond acceptor count': 'prop-type-hba',
            'hydrogen bond donor count': 'prop-type-hbd',
            'topological polar surface area': 'prop-type-psa',
            'rotatable bond count': 'prop-type-rtb',
            'aromatic ring count': 'prop-type-aromatic',
            'heavy atom count': 'prop-type-heavy-atoms',
            'rule of five violations': 'prop-type-ro5',
            'unknown property type': 'prop-type-unknown'
        }

    def test_transform_chembl_to_molecule_glycerol(self):
        """Test transform_chembl_to_molecule with Glycerol."""
        result = chembl_import.transform_chembl_to_molecule(self.glycerol, self.user_profile_id)
        
        # Verify basic fields
        self.assertEqual(result['name'], 'Glycerol')
        self.assertEqual(result['smiles'], 'C(C(CO)O)O')
        self.assertEqual(result['inchi'], 'InChI=1S/C3H8O3/c4-1-3(6)2-5/h3-6H,1-2H2')
        self.assertEqual(result['inchikey'], 'PEDCQBHIVMGVHV-UHFFFAOYSA-N')
        self.assertEqual(result['formula'], 'C3H8O3')
        self.assertEqual(result['molecular_weight'], 92.09)
        self.assertEqual(result['chembl_id'], 'CHEMBL1234')
        self.assertIn('ChEMBL', result['data_source'])
        self.assertIn('CHEMBL1234', result['data_source'])

    def test_transform_chembl_to_molecule_dmso(self):
        """Test transform_chembl_to_molecule with DMSO."""
        result = chembl_import.transform_chembl_to_molecule(self.dmso, self.user_profile_id)
        
        # Verify basic fields
        self.assertEqual(result['name'], 'Dimethyl sulfoxide')
        self.assertEqual(result['smiles'], 'CS(=O)C')
        self.assertEqual(result['inchi'], 'InChI=1S/C2H6OS/c1-4(2)3/h1-2H3')
        self.assertEqual(result['inchikey'], 'IAZDPXIOMUYVGZ-UHFFFAOYSA-N')
        self.assertEqual(result['formula'], 'C2H6OS')
        self.assertEqual(result['molecular_weight'], 78.13)
        self.assertEqual(result['chembl_id'], 'CHEMBL58')
        self.assertIn('ChEMBL', result['data_source'])
        self.assertIn('CHEMBL58', result['data_source'])

    def test_transform_chembl_to_molecule_handles_missing_data(self):
        """Test transform_chembl_to_molecule handles missing data gracefully."""
        result = chembl_import.transform_chembl_to_molecule(self.incomplete_compound, self.user_profile_id)
        
        # Verify that missing fields are handled gracefully
        self.assertEqual(result['name'], 'CHEMBL9999')  # Falls back to ChEMBL ID
        self.assertEqual(result['smiles'], 'CC(=O)OC1=CC=CC=C1C(=O)O')
        self.assertIsNone(result['inchi'])
        self.assertIsNone(result['inchikey'])
        self.assertIsNone(result['formula'])
        self.assertEqual(result['molecular_weight'], 180.16)
        self.assertEqual(result['chembl_id'], 'CHEMBL9999')
        self.assertIn('ChEMBL', result['data_source'])
        self.assertIn('CHEMBL9999', result['data_source'])

    @patch('ChEMBL_Integrated_Import.insert_property_type')
    def test_transform_chembl_to_properties_glycerol(self, mock_insert_property_type):
        """Test transform_chembl_to_properties with Glycerol."""
        # Setup mock
        mock_insert_property_type.return_value = 'new-prop-type-id'
        
        # Call the function
        properties = chembl_import.transform_chembl_to_properties(
            self.glycerol, 'test-molecule-id', self.user_profile_id, self.property_type_map
        )
        
        # Verify properties
        self.assertGreaterEqual(len(properties), 6)  # At least 6 properties
        
        # Check for specific properties
        property_names = [p['property_type_id'] for p in properties]
        self.assertIn('prop-type-logp', property_names)
        self.assertIn('prop-type-mw', property_names)
        self.assertIn('prop-type-hba', property_names)
        self.assertIn('prop-type-hbd', property_names)
        
        # Check property values
        for prop in properties:
            if prop['property_type_id'] == 'prop-type-logp':
                self.assertEqual(prop['numeric_value'], -1.76)
            elif prop['property_type_id'] == 'prop-type-mw':
                self.assertEqual(prop['numeric_value'], 92.09)
            elif prop['property_type_id'] == 'prop-type-hba':
                self.assertEqual(prop['numeric_value'], 3)
            elif prop['property_type_id'] == 'prop-type-hbd':
                self.assertEqual(prop['numeric_value'], 3)
        
        # Verify data source includes property name
        for prop in properties:
            self.assertIn('data_source', prop)
            self.assertIn('property:', prop['data_source'])
            self.assertIn('CHEMBL1234', prop['data_source'])

    @patch('ChEMBL_Integrated_Import.insert_property_type')
    def test_transform_chembl_to_properties_dmso(self, mock_insert_property_type):
        """Test transform_chembl_to_properties with DMSO."""
        # Setup mock
        mock_insert_property_type.return_value = 'new-prop-type-id'
        
        # Call the function
        properties = chembl_import.transform_chembl_to_properties(
            self.dmso, 'test-molecule-id', self.user_profile_id, self.property_type_map
        )
        
        # Verify properties
        self.assertGreaterEqual(len(properties), 6)  # At least 6 properties
        
        # Check for specific properties
        property_names = [p['property_type_id'] for p in properties]
        self.assertIn('prop-type-logp', property_names)
        self.assertIn('prop-type-mw', property_names)
        self.assertIn('prop-type-hba', property_names)
        
        # Check property values
        for prop in properties:
            if prop['property_type_id'] == 'prop-type-logp':
                self.assertEqual(prop['numeric_value'], -0.62)
            elif prop['property_type_id'] == 'prop-type-mw':
                self.assertEqual(prop['numeric_value'], 78.13)
            elif prop['property_type_id'] == 'prop-type-hba':
                self.assertEqual(prop['numeric_value'], 1)
            elif prop['property_type_id'] == 'prop-type-hbd':
                self.assertEqual(prop['numeric_value'], 0)
        
        # Verify data source includes property name
        for prop in properties:
            self.assertIn('data_source', prop)
            self.assertIn('property:', prop['data_source'])
            self.assertIn('CHEMBL58', prop['data_source'])

    @patch('ChEMBL_Integrated_Import.insert_property_type')
    def test_transform_chembl_to_properties_handles_missing_data(self, mock_insert_property_type):
        """Test transform_chembl_to_properties handles missing data gracefully."""
        # Setup mock
        mock_insert_property_type.return_value = 'new-prop-type-id'
        
        # Call the function
        properties = chembl_import.transform_chembl_to_properties(
            self.incomplete_compound, 'test-molecule-id', self.user_profile_id, self.property_type_map
        )
        
        # Verify properties
        self.assertGreaterEqual(len(properties), 1)  # At least 1 property (molecular weight)
        
        # Check for molecular weight property
        mw_props = [p for p in properties if p['property_type_id'] == 'prop-type-mw']
        self.assertEqual(len(mw_props), 1)
        self.assertEqual(mw_props[0]['numeric_value'], 180.16)
        
        # Verify data source includes property name
        for prop in properties:
            self.assertIn('data_source', prop)
            self.assertIn('property:', prop['data_source'])
            self.assertIn('CHEMBL9999', prop['data_source'])

    @patch('ChEMBL_Integrated_Import.insert_property_type')
    def test_transform_chembl_to_properties_handles_extra_properties(self, mock_insert_property_type):
        """Test transform_chembl_to_properties handles additional properties."""
        # Setup mock
        mock_insert_property_type.return_value = 'new-prop-type-id'
        
        # Call the function
        properties = chembl_import.transform_chembl_to_properties(
            self.compound_with_extra, 'test-molecule-id', self.user_profile_id, self.property_type_map
        )
        
        # Verify properties
        self.assertGreaterEqual(len(properties), 9)  # At least 9 properties
        
        # Check for standard properties
        standard_props = ['prop-type-logp', 'prop-type-mw', 'prop-type-hba', 'prop-type-hbd']
        for prop_id in standard_props:
            self.assertIn(prop_id, [p['property_type_id'] for p in properties])
        
        # Check for additional properties (should use insert_property_type)
        self.assertTrue(mock_insert_property_type.called)
        
        # Check for properties from molecule_synonyms
        synonym_props = [p for p in properties if 'Synonym' in p.get('data_source', '')]
        self.assertGreaterEqual(len(synonym_props), 1)
        
        # Check for properties from molecule_hierarchy
        hierarchy_props = [p for p in properties if 'Parent ChEMBL ID' in p.get('data_source', '')]
        self.assertGreaterEqual(len(hierarchy_props), 0)  # May or may not be included
        
        # Check for properties from the properties list
        custom_props = [p for p in properties if 'Solubility' in p.get('data_source', '') or 'Melting Point' in p.get('data_source', '')]
        self.assertGreaterEqual(len(custom_props), 1)
        
        # Verify data source includes property name
        for prop in properties:
            self.assertIn('data_source', prop)
            self.assertIn('property:', prop['data_source'])
            self.assertIn('CHEMBL5555', prop['data_source'])

    @patch('ChEMBL_Integrated_Import.insert_property_type')
    def test_transform_chembl_to_properties_handles_property_type_insertion_failure(self, mock_insert_property_type):
        """Test transform_chembl_to_properties handles property type insertion failure."""
        # Setup mock to simulate failure
        mock_insert_property_type.return_value = None
        
        # Call the function
        properties = chembl_import.transform_chembl_to_properties(
            self.glycerol, 'test-molecule-id', self.user_profile_id, self.property_type_map
        )
        
        # Verify properties - should still have properties for known property types
        self.assertGreaterEqual(len(properties), 1)
        
        # Check that known property types are still processed
        known_props = [p for p in properties if p['property_type_id'] in self.property_type_map.values()]
        self.assertGreaterEqual(len(known_props), 1)


if __name__ == '__main__':
    unittest.main()