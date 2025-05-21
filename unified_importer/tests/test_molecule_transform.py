"""
Test cases for the molecule_transform module.

This module tests the functionality of the molecule transformation utilities,
including standardization, property calculation, and cross-reference resolution.
"""

import pytest
import asyncio
import logging
from typing import Dict, Any

from ..transforms.molecule_transform import (
    MoleculeTransformer, 
    normalize_smiles,
    get_molecule_fingerprint, 
    calculate_similarity
)

# Test data
TEST_MOLECULE = {
    'name': 'Aspirin',
    'smiles': 'CC(=O)OC1=CC=CC=C1C(=O)O',
    'chembl_id': 'CHEMBL25',
    'data_source': 'ChEMBL'
}

TEST_MOLECULE_PUBCHEM = {
    'name': 'Acetylsalicylic acid',
    'smiles': 'CC(=O)OC1=CC=CC=C1C(=O)O',
    'pubchem_cid': '2244',
    'data_source': 'PubChem'
}

TEST_MIXTURE = {
    'name': 'Sodium Chloride Solution',
    'smiles': '[Na+].[Cl-].[H]O[H]',
    'data_source': 'Test'
}


class TestMoleculeTransform:
    """Test the molecule transformation functionality."""
    
    @pytest.fixture
    def transformer(self):
        """Create a MoleculeTransformer instance for testing."""
        logger = logging.getLogger('test')
        config = {
            'api_delay': 0.1,  # Fast for testing
            'resolve_cross_references': True,
            'handle_mixtures': True,
            'batch_size': 5
        }
        return MoleculeTransformer(logger=logger, config=config)
    
    @pytest.mark.asyncio
    async def test_standardize_molecule(self, transformer):
        """Test molecule standardization."""
        result = await transformer.standardize_molecule(TEST_MOLECULE, resolve_ids=False)
        
        # Check that basic fields are preserved
        assert result['name'] == 'Aspirin'
        assert result['chembl_id'] == 'CHEMBL25'
        assert result['data_source'] == 'ChEMBL'
        
        # Check that identifiers are added/standardized
        assert 'inchi' in result
        assert 'inchikey' in result
        assert 'smiles' in result
        assert 'formula' in result
        assert 'molecular_weight' in result
        
        # Check that properties are calculated
        assert 'properties' in result
        assert 'calculated_properties' in result
        assert len(result['calculated_properties']) > 0
        
        # Check specific properties
        assert 'LogP' in result['properties']
        assert 'TPSA' in result['properties']
        assert 'HBondDonorCount' in result['properties']
        assert 'HBondAcceptorCount' in result['properties']
    
    @pytest.mark.asyncio
    async def test_mixture_detection(self, transformer):
        """Test detection of mixtures."""
        result = await transformer.standardize_molecule(TEST_MIXTURE, resolve_ids=False)
        
        # Check that mixture is detected
        assert 'is_mixture' in result
        assert result['is_mixture'] is True
        
        # Check that components are extracted
        assert 'mixture_components' in result
        assert len(result['mixture_components']) == 3  # Na+, Cl-, H2O
        
        # Check that each component has basic fields
        for component in result['mixture_components']:
            assert 'smiles' in component
            assert 'molecular_weight' in component
            assert 'formula' in component
            assert 'name' in component
    
    @pytest.mark.asyncio
    async def test_merge_molecule_data(self, transformer):
        """Test merging of molecule data."""
        # Create test data with different properties
        primary = TEST_MOLECULE.copy()
        primary['formula'] = 'C9H8O4'
        primary['synonyms'] = ['Aspirin', 'Acetylsalicylic acid']
        
        secondary = TEST_MOLECULE_PUBCHEM.copy()
        secondary['formula'] = 'C9H8O4'
        secondary['pubchem_cid'] = '2244'
        secondary['synonyms'] = ['ASA', '2-acetoxybenzoic acid', 'Acetylsalicylic acid']
        
        # Merge the data
        result = await transformer.merge_molecule_data(primary, secondary)
        
        # Check that fields are properly merged
        assert result['name'] == 'Aspirin'  # Kept from primary
        assert result['pubchem_cid'] == '2244'  # Added from secondary
        assert result['chembl_id'] == 'CHEMBL25'  # Kept from primary
        
        # Check that synonyms are merged without duplicates
        assert 'synonyms' in result
        assert len(result['synonyms']) == 4
        assert 'Aspirin' in result['synonyms']
        assert 'Acetylsalicylic acid' in result['synonyms']
        assert 'ASA' in result['synonyms']
        assert '2-acetoxybenzoic acid' in result['synonyms']
    
    @pytest.mark.asyncio
    async def test_standardize_molecules_batch(self, transformer):
        """Test batch processing of molecules."""
        molecules = [TEST_MOLECULE.copy(), TEST_MOLECULE_PUBCHEM.copy(), TEST_MIXTURE.copy()]
        
        # Process the batch
        results = await transformer.standardize_molecules_batch(molecules, resolve_ids=False)
        
        # Check that all molecules are processed
        assert len(results) == 3
        
        # Check that each molecule has been standardized
        for result in results:
            assert 'properties' in result
            assert 'calculated_properties' in result
    
    @pytest.mark.asyncio
    async def test_cross_reference_resolution(self, transformer):
        """
        Test cross-reference resolution between ChEMBL and PubChem.
        
        Note: This test requires internet connection and may be slow or fail if
        the external APIs are unavailable. Consider mocking the API responses
        for reliable unit testing.
        """
        # Skip this test in offline environments
        pytest.skip("Skipping online API test - mock this for CI environments")
        
        # Test ChEMBL to PubChem resolution
        chembl_mol = TEST_MOLECULE.copy()
        result = await transformer.standardize_molecule(chembl_mol, resolve_ids=True)
        
        # Check that PubChem CID is resolved
        assert 'pubchem_cid' in result
        assert result['pubchem_cid'] is not None
        
        # Test PubChem to ChEMBL resolution
        pubchem_mol = TEST_MOLECULE_PUBCHEM.copy()
        result = await transformer.standardize_molecule(pubchem_mol, resolve_ids=True)
        
        # Check that ChEMBL ID is resolved
        assert 'chembl_id' in result
        assert result['chembl_id'] is not None


def test_normalize_smiles():
    """Test SMILES normalization."""
    # Test with valid SMILES
    smiles = 'CC(=O)OC1=CC=CC=C1C(=O)O'
    normalized = normalize_smiles(smiles)
    assert normalized is not None
    assert len(normalized) > 0
    
    # Test with alternative representations of the same molecule
    alt_smiles = 'OC(=O)c1ccccc1OC(=O)C'
    alt_normalized = normalize_smiles(alt_smiles)
    
    # Both should normalize to the same canonical form
    assert normalized == alt_normalized


def test_get_molecule_fingerprint():
    """Test fingerprint generation."""
    smiles = 'CC(=O)OC1=CC=CC=C1C(=O)O'
    fp = get_molecule_fingerprint(smiles)
    
    # Check that fingerprint is generated
    assert fp is not None
    assert len(fp) > 0


def test_calculate_similarity():
    """Test similarity calculation."""
    # Test with identical molecules
    smiles1 = 'CC(=O)OC1=CC=CC=C1C(=O)O'  # Aspirin
    similarity = calculate_similarity(smiles1, smiles1)
    assert similarity == 1.0  # Identical molecules have similarity 1.0
    
    # Test with similar molecules
    smiles2 = 'OC(=O)C1=CC=CC=C1O'  # Salicylic acid (similar to aspirin)
    similarity = calculate_similarity(smiles1, smiles2)
    assert 0 < similarity < 1  # Similar molecules have similarity between 0 and 1
    
    # Test with different fingerprint types
    similarity_maccs = calculate_similarity(smiles1, smiles2, fingerprint_type='maccs')
    similarity_topo = calculate_similarity(smiles1, smiles2, fingerprint_type='topological')
    
    # Different fingerprint types should give different similarity values
    assert similarity_maccs != similarity_topo