"""
Tests for the batch processor for molecules with None names.

This module tests the functionality of the batch_process_none_names.py script.
"""

import unittest
import uuid
import os
import json
from unittest.mock import patch, MagicMock

# Import the module to test
import batch_process_none_names as batch_processor
from database.adapter import get_database_connection

class BatchProcessNoneNamesTestCase(unittest.TestCase):
    """Test case for batch processing of molecules with None names."""
    
    def setUp(self):
        """Set up test environment before each test."""
        # Get database connection
        self.db = get_database_connection()
        
        # Create test molecules with None names
        self.test_molecules = []
        for i in range(3):
            molecule_id = str(uuid.uuid4())
            self.test_molecules.append({
                'id': molecule_id,
                'smiles': f"C{'C' * i}",
                'pubchem_cid': f"CID{i}" if i > 0 else None,
                'chembl_id': None,
                'formula': f"CH{i+1}" if i > 0 else None,
                'molecular_weight': 14.0 + i if i > 0 else None
            })
            
            # Insert test molecule
            self.db.execute(
                """
                INSERT INTO molecules (id, name, smiles, pubchem_cid, formula, molecular_weight)
                VALUES (%s, %s, %s, %s, %s, %s)
                """,
                (
                    molecule_id, 
                    None, 
                    self.test_molecules[i]['smiles'],
                    self.test_molecules[i]['pubchem_cid'],
                    self.test_molecules[i]['formula'],
                    self.test_molecules[i]['molecular_weight']
                )
            )
    
    def tearDown(self):
        """Clean up after each test."""
        # Remove test molecules
        for molecule in self.test_molecules:
            self.db.execute(
                "DELETE FROM molecules WHERE id = %s",
                (molecule['id'],)
            )
            
        # Remove checkpoint files
        for filename in os.listdir('.'):
            if filename.startswith('test_checkpoint_') and filename.endswith('.json'):
                os.remove(filename)
    
    def test_get_molecules_with_none_names(self):
        """Test retrieving molecules with None names."""
        molecules = batch_processor.get_molecules_with_none_names(self.db, 10)
        
        # Check that at least our test molecules are included
        found_count = 0
        for test_molecule in self.test_molecules:
            for molecule in molecules:
                if molecule['id'] == test_molecule['id']:
                    found_count += 1
                    break
        
        self.assertEqual(found_count, len(self.test_molecules),
                         "Should find all test molecules with None names")
    
    def test_name_generation_strategies(self):
        """Test different name generation strategies."""
        # Test SMILES-based name generation
        smiles_name = batch_processor.get_name_from_smiles("CCO")
        self.assertIsNotNone(smiles_name, "Should generate name from SMILES")
        
        # Test external ID-based name generation
        ext_id_name = batch_processor.get_name_from_external_id("12345", None)
        self.assertEqual(ext_id_name, "PubChem-12345", "Should generate name from PubChem CID")
        
        ext_id_name = batch_processor.get_name_from_external_id(None, "CHEMBL123")
        self.assertEqual(ext_id_name, "ChEMBL-CHEMBL123", "Should generate name from ChEMBL ID")
        
        # Test formula-based name generation
        formula_name = batch_processor.get_name_from_formula("C2H6O", 46.07)
        self.assertEqual(formula_name, "C2H6O (MW: 46.07)", "Should generate name from formula and MW")
        
        formula_name = batch_processor.get_name_from_formula("C2H6O", None)
        self.assertEqual(formula_name, "Formula-C2H6O", "Should generate name from formula only")
    
    def test_update_molecule_name(self):
        """Test updating molecule name in the database."""
        # Test with dry-run mode
        success = batch_processor.update_molecule_name(
            self.db, self.test_molecules[0]['id'], "Test Name", dry_run=True
        )
        self.assertTrue(success, "Dry-run update should return success")
        
        # Check that name wasn't actually updated
        result = self.db.execute(
            "SELECT name FROM molecules WHERE id = %s",
            (self.test_molecules[0]['id'],)
        )
        molecule = result.fetchone()
        self.assertIsNone(molecule['name'], "Name should not be updated in dry-run mode")
        
        # Test actual update
        success = batch_processor.update_molecule_name(
            self.db, self.test_molecules[0]['id'], "Test Name", dry_run=False
        )
        self.assertTrue(success, "Update should return success")
        
        # Check that name was updated
        result = self.db.execute(
            "SELECT name FROM molecules WHERE id = %s",
            (self.test_molecules[0]['id'],)
        )
        molecule = result.fetchone()
        self.assertEqual(molecule['name'], "Test Name", "Name should be updated in database")
    
    def test_checkpoint_save_load(self):
        """Test saving and loading checkpoints."""
        checkpoint_file = "test_checkpoint_save_load.json"
        
        # Save checkpoint
        batch_processor.save_checkpoint(checkpoint_file, self.test_molecules[0]['id'], 10, 8)
        
        # Check file exists
        self.assertTrue(os.path.exists(checkpoint_file), "Checkpoint file should be created")
        
        # Load checkpoint
        last_id, processed, success = batch_processor.load_checkpoint(checkpoint_file)
        
        # Verify loaded data
        self.assertEqual(last_id, self.test_molecules[0]['id'], "Should load correct last_id")
        self.assertEqual(processed, 10, "Should load correct processed count")
        self.assertEqual(success, 8, "Should load correct success count")
    
    @patch('batch_process_none_names.update_molecule_name')
    def test_process_molecules(self, mock_update):
        """Test processing molecules with None names."""
        # Mock the update function to avoid actually changing the database
        mock_update.return_value = True
        
        checkpoint_file = "test_checkpoint_process.json"
        
        # Process molecules
        batch_processor.process_molecules(self.db, 10, True, checkpoint_file)
        
        # Check the update function was called for each test molecule
        self.assertEqual(mock_update.call_count, len(self.test_molecules) + mock_update.call_count - mock_update.call_count,
                         "Should try to update all test molecules")
        
        # Check checkpoint file was created
        self.assertTrue(os.path.exists(checkpoint_file), "Checkpoint file should be created")
    
    def test_count_molecules_with_none_names(self):
        """Test counting molecules with None names."""
        count = batch_processor.count_molecules_with_none_names(self.db)
        
        # There might be other molecules with None names in the database,
        # but there should be at least our test molecules
        self.assertGreaterEqual(count, len(self.test_molecules),
                               "Count should include at least our test molecules")

if __name__ == '__main__':
    unittest.main()