"""
Tests to demonstrate mock objects.
"""

import pytest
from tests.fixtures.mocks.supabase import patch_supabase_client, MockSupabase
from tests.fixtures.mocks.rdkit import patch_rdkit, MockRDKit, MockMolecule

def test_supabase_mock():
    """Test the Supabase mock."""
    with patch_supabase_client() as mock_db:
        # Add test data
        mock_db.add_test_data('molecules', [
            {'id': 'mol-1', 'name': 'Test Molecule 1'},
            {'id': 'mol-2', 'name': 'Test Molecule 2'}
        ])
        
        # Test select
        result = mock_db.table('molecules').select().execute()
        assert len(result.data) == 2
        assert result.data[0]['name'] == 'Test Molecule 1'
        
        # Test insert
        mock_db.table('molecules').insert({'name': 'Test Molecule 3'})
        result = mock_db.table('molecules').select().execute()
        assert len(result.data) == 3
        
        # Test filter
        result = mock_db.table('molecules').select().eq('name', 'Test Molecule 2').execute()
        assert len(result.data) == 1
        assert result.data[0]['name'] == 'Test Molecule 2'
        
        # Test update
        mock_db.table('molecules').update({'formula': 'H2O'}).eq('name', 'Test Molecule 2').execute()
        result = mock_db.table('molecules').select().eq('name', 'Test Molecule 2').execute()
        assert result.data[0]['formula'] == 'H2O'
        
        # Test delete
        mock_db.table('molecules').delete().eq('name', 'Test Molecule 1').execute()
        result = mock_db.table('molecules').select().execute()
        assert len(result.data) == 2
        assert 'Test Molecule 1' not in [item['name'] for item in result.data]
        
        # Test RPC
        result = mock_db.rpc('has_table', {'table_name': 'molecules'}).execute()
        assert result.data[0] is True
        
        # Test SQL
        result = mock_db.sql('SELECT * FROM molecules').execute()
        assert len(result.data) == 2

def test_rdkit_mock():
    """Test the RDKit mock."""
    with patch_rdkit() as rdkit:
        # Create molecule from SMILES
        mol = rdkit.Chem.MolFromSmiles('CCO')
        assert mol is not None
        assert mol.smiles == 'CCO'
        
        # Test SMILES conversion
        smiles = rdkit.Chem.MolToSmiles(mol)
        assert smiles == 'CCO'
        
        # Test property calculation
        weight = rdkit.Descriptors.MolWt(mol)
        assert weight == 180.16
        
        logp = rdkit.Descriptors.MolLogP(mol)
        assert logp == -0.5
        
        # Test invalid input
        invalid_mol = rdkit.Chem.MolFromSmiles('invalid')
        assert invalid_mol is None

def test_rdkit_fixture(mock_rdkit, molecule_factory):
    """Test the RDKit fixture and molecule factory."""
    # Create molecule using the mock_rdkit fixture
    mol = mock_rdkit.Chem.MolFromSmiles('CCO')
    assert mol is not None
    
    # Create a custom molecule with the factory
    custom_mol = molecule_factory(
        smiles='CCC',
        mol_weight=44.1,
        formula='C2H6'
    )
    
    assert custom_mol.smiles == 'CCC'
    assert custom_mol.mol_weight == 44.1
    assert custom_mol.formula == 'C2H6'
    
    # Test descriptor calculation
    weight = mock_rdkit.Descriptors.MolWt(custom_mol)
    assert weight == 44.1

def test_combined_mocks():
    """Test using both mocks together."""
    with patch_supabase_client() as mock_db, patch_rdkit() as rdkit:
        # Add test molecules
        mock_db.add_test_data('molecules', [
            {'id': 'mol-1', 'name': 'Ethanol', 'smiles': 'CCO'}
        ])
        
        # Query molecule
        result = mock_db.table('molecules').select().eq('name', 'Ethanol').execute()
        molecule_data = result.data[0]
        
        # Use RDKit to process the molecule
        mol = rdkit.Chem.MolFromSmiles(molecule_data['smiles'])
        assert mol is not None
        
        # Calculate properties
        weight = rdkit.Descriptors.MolWt(mol)
        assert weight == 180.16
        
        # Use the result to update the database
        mock_db.table('molecules').update({
            'molecular_weight': weight
        }).eq('id', 'mol-1').execute()
        
        # Verify the update
        result = mock_db.table('molecules').select().eq('id', 'mol-1').execute()
        assert result.data[0]['molecular_weight'] == 180.16