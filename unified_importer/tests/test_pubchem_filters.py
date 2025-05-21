"""
Tests for PubChem data source property filters.

This module tests the property filtering capabilities of the PubChemDataSource class.
"""

import unittest
import asyncio
import logging
from unittest.mock import MagicMock, AsyncMock, patch
from typing import Dict, List, Any, Optional

from ..sources.pubchem_source import PubChemDataSource
from ..core.database import DatabaseOperations


class TestPubChemFilters(unittest.TestCase):
    """Test cases for PubChem property filters."""

    def setUp(self):
        """Set up test fixtures."""
        # Create a mock database
        self.db_mock = MagicMock(spec=DatabaseOperations)
        
        # Set up logger
        self.logger = logging.getLogger('test_pubchem_filters')
        self.logger.setLevel(logging.DEBUG)
        
        # Create PubChemDataSource with property filter configuration
        self.config = {
            'sources': {
                'pubchem': {
                    'property_filters': [
                        {
                            'name': 'test_terms_filter',
                            'description': 'Test filter with terms',
                            'terms': ['glycerol', 'dmso']
                        },
                        {
                            'name': 'test_property_filter',
                            'description': 'Test filter with property constraints',
                            'molecular_weight_max': 200,
                            'logp_max': 1.0
                        },
                        {
                            'name': 'test_combined_filter',
                            'description': 'Test filter with terms and properties',
                            'terms': ['cryoprotect'],
                            'molecular_weight_max': 300,
                            'hbond_donors_min': 2
                        }
                    ]
                }
            }
        }
        
        # Create the data source
        self.pubchem = PubChemDataSource(
            db_operations=self.db_mock,
            config=self.config,
            logger=self.logger
        )
        
        # Create some test compound data
        self.test_compound_small = {
            'CID': 1001,
            'props': [
                {'urn': {'label': 'IUPAC Name'}, 'value': {'sval': 'Small Molecule'}},
                {'urn': {'label': 'Molecular Weight'}, 'value': {'fval': 150.2}},
                {'urn': {'label': 'XLogP'}, 'value': {'fval': 0.8}},
                {'urn': {'label': 'HBondDonorCount'}, 'value': {'ival': 3}},
                {'urn': {'label': 'HBondAcceptorCount'}, 'value': {'ival': 2}},
                {'urn': {'label': 'RotatableBondCount'}, 'value': {'ival': 4}}
            ]
        }
        
        self.test_compound_large = {
            'CID': 1002,
            'props': [
                {'urn': {'label': 'IUPAC Name'}, 'value': {'sval': 'Large Molecule'}},
                {'urn': {'label': 'Molecular Weight'}, 'value': {'fval': 450.6}},
                {'urn': {'label': 'XLogP'}, 'value': {'fval': 3.2}},
                {'urn': {'label': 'HBondDonorCount'}, 'value': {'ival': 1}},
                {'urn': {'label': 'HBondAcceptorCount'}, 'value': {'ival': 5}},
                {'urn': {'label': 'RotatableBondCount'}, 'value': {'ival': 8}}
            ]
        }
        
        self.test_compound_glycerol = {
            'CID': 1003,
            'props': [
                {'urn': {'label': 'IUPAC Name'}, 'value': {'sval': 'Glycerol'}},
                {'urn': {'label': 'Molecular Weight'}, 'value': {'fval': 92.09}},
                {'urn': {'label': 'XLogP'}, 'value': {'fval': -1.76}},
                {'urn': {'label': 'HBondDonorCount'}, 'value': {'ival': 3}},
                {'urn': {'label': 'HBondAcceptorCount'}, 'value': {'ival': 3}},
                {'urn': {'label': 'RotatableBondCount'}, 'value': {'ival': 2}}
            ],
            'synonyms': ['glycerol', 'glycerin', '1,2,3-Propanetriol']
        }
        
        self.test_compound_cryoprotect = {
            'CID': 1004,
            'props': [
                {'urn': {'label': 'IUPAC Name'}, 'value': {'sval': 'Cryoprotectant X'}},
                {'urn': {'label': 'Molecular Weight'}, 'value': {'fval': 250.3}},
                {'urn': {'label': 'XLogP'}, 'value': {'fval': 1.2}},
                {'urn': {'label': 'HBondDonorCount'}, 'value': {'ival': 4}},
                {'urn': {'label': 'HBondAcceptorCount'}, 'value': {'ival': 6}},
                {'urn': {'label': 'RotatableBondCount'}, 'value': {'ival': 5}}
            ],
            'synonyms': ['cryoprotectant', 'freezing protectant']
        }

    async def _get_properties_for_filter_mock(self, compound_data: Dict[str, Any]) -> Dict[str, float]:
        """Mock implementation of _get_properties_for_filter."""
        cid = compound_data.get('CID')
        
        # Define property mappings and test values
        properties = {}
        
        if cid == 1001:  # Small molecule
            properties = {
                'molecular_weight': 150.2,
                'logp': 0.8,
                'hbond_donors': 3,
                'hbond_acceptors': 2,
                'rotatable_bonds': 4,
                'tpsa': 40.0
            }
        elif cid == 1002:  # Large molecule
            properties = {
                'molecular_weight': 450.6,
                'logp': 3.2,
                'hbond_donors': 1,
                'hbond_acceptors': 5,
                'rotatable_bonds': 8,
                'tpsa': 85.0
            }
        elif cid == 1003:  # Glycerol
            properties = {
                'molecular_weight': 92.09,
                'logp': -1.76,
                'hbond_donors': 3,
                'hbond_acceptors': 3,
                'rotatable_bonds': 2,
                'tpsa': 60.7
            }
        elif cid == 1004:  # Cryoprotectant
            properties = {
                'molecular_weight': 250.3,
                'logp': 1.2,
                'hbond_donors': 4,
                'hbond_acceptors': 6,
                'rotatable_bonds': 5,
                'tpsa': 110.2
            }
            
        return properties

    def test_initialize_property_filters(self):
        """Test that property filters are correctly initialized from config."""
        # Check the number of filters loaded
        self.assertEqual(len(self.pubchem.property_filters), 3)
        
        # Check filter names
        filter_names = [f['name'] for f in self.pubchem.property_filters]
        self.assertIn('test_terms_filter', filter_names)
        self.assertIn('test_property_filter', filter_names)
        self.assertIn('test_combined_filter', filter_names)
        
        # Check filter content
        terms_filter = next(f for f in self.pubchem.property_filters if f['name'] == 'test_terms_filter')
        self.assertEqual(len(terms_filter['terms']), 2)
        self.assertIn('glycerol', terms_filter['terms'])
        self.assertIn('dmso', terms_filter['terms'])
        
        property_filter = next(f for f in self.pubchem.property_filters if f['name'] == 'test_property_filter')
        self.assertEqual(property_filter['molecular_weight_max'], 200)
        self.assertEqual(property_filter['logp_max'], 1.0)

    def test_validate_property_filter(self):
        """Test property filter validation."""
        # Valid filter with terms
        valid_terms_filter = {
            'name': 'valid_terms',
            'terms': ['term1', 'term2']
        }
        try:
            self.pubchem._validate_property_filter(valid_terms_filter)
        except ValueError:
            self.fail("Validation failed for valid terms filter")
            
        # Valid filter with properties
        valid_props_filter = {
            'name': 'valid_props',
            'molecular_weight_max': 200,
            'logp_min': -2.0
        }
        try:
            self.pubchem._validate_property_filter(valid_props_filter)
        except ValueError:
            self.fail("Validation failed for valid properties filter")
            
        # Invalid filter - no name
        invalid_no_name = {
            'terms': ['term1']
        }
        with self.assertRaises(ValueError):
            self.pubchem._validate_property_filter(invalid_no_name)
            
        # Invalid filter - empty name
        invalid_empty_name = {
            'name': '',
            'terms': ['term1']
        }
        with self.assertRaises(ValueError):
            self.pubchem._validate_property_filter(invalid_empty_name)
            
        # Invalid filter - non-list terms
        invalid_terms = {
            'name': 'invalid',
            'terms': 'not_a_list'
        }
        with self.assertRaises(ValueError):
            self.pubchem._validate_property_filter(invalid_terms)
            
        # Invalid filter - empty terms list
        invalid_empty_terms = {
            'name': 'invalid',
            'terms': []
        }
        with self.assertRaises(ValueError):
            self.pubchem._validate_property_filter(invalid_empty_terms)
            
        # Invalid filter - non-numeric property constraint
        invalid_constraint = {
            'name': 'invalid',
            'molecular_weight_max': 'not_a_number'
        }
        with self.assertRaises(ValueError):
            self.pubchem._validate_property_filter(invalid_constraint)
            
        # Invalid filter - no terms or constraints
        invalid_empty = {
            'name': 'invalid'
        }
        with self.assertRaises(ValueError):
            self.pubchem._validate_property_filter(invalid_empty)

    @patch('unified_importer.sources.pubchem_source.PubChemDataSource._get_properties_for_filter')
    async def test_filter_compound_by_properties(self, mock_get_properties):
        """Test filtering compounds by properties."""
        # Set up the mock
        mock_get_properties.side_effect = self._get_properties_for_filter_mock
        
        # Test with no filters
        self.pubchem.property_filters = []
        result = await self.pubchem.filter_compound_by_properties(self.test_compound_small)
        self.assertTrue(result)  # Should pass with no filters
        
        # Restore filters
        self.pubchem.property_filters = self.config['sources']['pubchem']['property_filters']
        
        # Test terms filter with glycerol
        result = await self.pubchem._matches_filter(
            self.test_compound_glycerol,
            self.pubchem.property_filters[0]  # test_terms_filter
        )
        self.assertTrue(result)
        
        # Test terms filter with non-matching compound
        result = await self.pubchem._matches_filter(
            self.test_compound_large,
            self.pubchem.property_filters[0]  # test_terms_filter
        )
        self.assertFalse(result)
        
        # Test property filter with small molecule
        result = await self.pubchem._matches_filter(
            self.test_compound_small,
            self.pubchem.property_filters[1]  # test_property_filter
        )
        self.assertTrue(result)
        
        # Test property filter with large molecule
        result = await self.pubchem._matches_filter(
            self.test_compound_large,
            self.pubchem.property_filters[1]  # test_property_filter
        )
        self.assertFalse(result)
        
        # Test combined filter with cryoprotectant
        result = await self.pubchem._matches_filter(
            self.test_compound_cryoprotect,
            self.pubchem.property_filters[2]  # test_combined_filter
        )
        self.assertTrue(result)
        
        # Test overall filter_compound_by_properties method
        
        # Small molecule should match property filter
        result = await self.pubchem.filter_compound_by_properties(self.test_compound_small)
        self.assertTrue(result)
        
        # Glycerol should match terms filter
        result = await self.pubchem.filter_compound_by_properties(self.test_compound_glycerol)
        self.assertTrue(result)
        
        # Large molecule should not match any filter
        result = await self.pubchem.filter_compound_by_properties(self.test_compound_large)
        self.assertFalse(result)

    @patch('unified_importer.sources.pubchem_source.PubChemDataSource._make_api_request')
    @patch('unified_importer.sources.pubchem_source.PubChemDataSource.get_compound_batch')
    @patch('unified_importer.sources.pubchem_source.PubChemDataSource._matches_filter')
    async def test_search_compounds_with_filters(self, mock_matches_filter, mock_get_batch, mock_api_request):
        """Test searching compounds with filters applied."""
        # Mock API responses
        mock_api_request.side_effect = [
            # count response
            {'PC_Count': 2},
            # cids response
            {'IdentifierList': {'CID': [1001, 1002]}}
        ]
        
        # Mock batch retrieval
        mock_get_batch.return_value = [self.test_compound_small, self.test_compound_large]
        
        # Mock filter matching
        mock_matches_filter.side_effect = [True, False]  # Small passes, large fails
        
        # Call search_compounds with filtering
        results = await self.pubchem.search_compounds('test query', apply_filters=True)
        
        # Should have one compound that passes the filter
        self.assertEqual(len(results), 1)
        self.assertEqual(results[0], '1001')
        
        # Test without filtering
        # Reset mocks
        mock_api_request.side_effect = [
            # count response
            {'PC_Count': 2},
            # cids response
            {'IdentifierList': {'CID': [1001, 1002]}}
        ]
        
        # Call without applying filters
        results = await self.pubchem.search_compounds('test query', apply_filters=False)
        
        # Should have both compounds
        self.assertEqual(len(results), 2)
        self.assertIn('1001', results)
        self.assertIn('1002', results)

    @patch('unified_importer.sources.pubchem_source.PubChemDataSource.search_compounds')
    @patch('unified_importer.sources.pubchem_source.PubChemDataSource.get_compound_batch')
    @patch('unified_importer.sources.pubchem_source.PubChemDataSource._matches_filter')
    async def test_search_compounds_by_filter(self, mock_matches_filter, mock_get_batch, mock_search):
        """Test searching compounds by a specific property filter."""
        # Mock search results for terms
        mock_search.side_effect = [
            ['1001', '1003'],  # Results for 'glycerol'
            ['1002', '1004']   # Results for 'dmso'
        ]
        
        # Mock batch retrieval
        mock_get_batch.return_value = [
            self.test_compound_small,
            self.test_compound_glycerol,
            self.test_compound_large,
            self.test_compound_cryoprotect
        ]
        
        # Mock filter matching
        mock_matches_filter.side_effect = [True, True, False, False]
        
        # Search using the terms filter
        results = await self.pubchem.search_compounds_by_filter('test_terms_filter')
        
        # Should call search_compounds for both terms
        self.assertEqual(mock_search.call_count, 2)
        
        # Should have called _matches_filter for property constraints
        self.assertTrue(mock_matches_filter.called)
        
        # Should return matching compounds
        self.assertEqual(len(results), 2)


if __name__ == '__main__':
    unittest.main()