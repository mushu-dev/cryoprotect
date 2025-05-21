"""
Integration tests for the MolecularImporter class.

These tests verify that the MolecularImporter correctly integrates with
the molecule and property transformers.
"""

import unittest
import tempfile
import json
import os
import asyncio
from typing import Dict, Any, List
from unittest.mock import patch

from .mock_config import load_config
from ..main import MolecularImporter


class TestMolecularImporterIntegration(unittest.TestCase):
    """Integration tests for the MolecularImporter class."""
    
    def setUp(self):
        """Set up test fixtures."""
        # Create a temporary configuration file
        self.config_data = {
            "database": {
                "url": None,  # No actual database connection in tests
                "key": None,
                "use_connection_pool": False
            },
            "sources": {
                "chembl": {
                    "test": True,
                    "api_url": "https://test-chembl-url",
                    "use_client": False
                }
            },
            "transforms": {
                "molecule_transform": {
                    "resolve_cross_references": False,
                    "handle_mixtures": True,
                    "add_standardized_structure": True
                },
                "property_transform": {
                    "standardize_units": True,
                    "add_missing_properties": True,
                    "preferred_units": {
                        "MolecularWeight": "kDa"
                    }
                }
            },
            "logging": {
                "level": "DEBUG"
            }
        }
        
        # Write the config to a temporary file
        self.config_file = tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False)
        json.dump(self.config_data, self.config_file)
        self.config_file.close()
        
        # Create the importer with a patch to use our mock config loader
        with patch('unified_importer.main.load_config', new=load_config):
            self.importer = MolecularImporter(config_file=self.config_file.name)
    
    def tearDown(self):
        """Clean up test fixtures."""
        # Remove the temporary config file
        os.unlink(self.config_file.name)
    
    def test_transformers_initialization(self):
        """Test that the transformers are properly initialized with config."""
        # Check MoleculeTransformer initialization
        molecule_transformer = self.importer.molecule_transformer
        self.assertFalse(molecule_transformer.resolve_cross_references)
        self.assertTrue(molecule_transformer.handle_mixtures)
        self.assertTrue(molecule_transformer.config.get('add_standardized_structure'))
        
        # Check PropertyTransformer initialization
        property_transformer = self.importer.property_transformer
        self.assertTrue(property_transformer.standardize_units)
        self.assertTrue(property_transformer.add_missing_properties)
        self.assertEqual(property_transformer.preferred_units.get('MolecularWeight'), 'kDa')
    
    def test_molecule_transformer_sync(self):
        """Test sync wrapper for standardize_molecule."""
        # Create a mock molecule
        molecule_data = {
            'name': 'Test Molecule',
            'smiles': 'CCO',  # Ethanol
            'data_source': 'Test'
        }
        
        # Use run_until_complete to run the async function synchronously
        loop = asyncio.new_event_loop()
        asyncio.set_event_loop(loop)
        
        try:
            result = loop.run_until_complete(
                self.importer.standardize_molecule(molecule_data, resolve_ids=False)
            )
            
            # Basic checks
            self.assertEqual(result['name'], 'Test Molecule')
            self.assertEqual(result['smiles'], 'CCO')
            self.assertEqual(result['data_source'], 'Test')
        except ImportError:
            # If RDKit is not available, some functionality won't work
            pass
        finally:
            loop.close()
    
    def test_config_loading(self):
        """Test that the configuration is properly loaded."""
        self.assertIn('database', self.importer.config)
        self.assertIn('sources', self.importer.config)
        self.assertIn('transforms', self.importer.config)
        
        # Check specific values
        self.assertFalse(self.importer.config['database']['use_connection_pool'])
        
        # Check that the transforms section is correctly loaded
        transforms = self.importer.config['transforms']
        self.assertIn('molecule_transform', transforms)
        self.assertIn('property_transform', transforms)
        
        # Check that the property_transform config has the preferred_units
        property_transform = transforms['property_transform']
        self.assertIn('preferred_units', property_transform)
        self.assertEqual(property_transform['preferred_units']['MolecularWeight'], 'kDa')


if __name__ == '__main__':
    unittest.main()