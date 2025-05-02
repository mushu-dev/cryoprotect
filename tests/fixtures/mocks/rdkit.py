"""
RDKit mock implementation.

This module provides a mock implementation of RDKit
for testing without requiring the actual RDKit package.
"""

import pytest
from unittest.mock import MagicMock, patch
from contextlib import contextmanager
from typing import Dict, Any, Generator, Optional

class MockDescriptors:
    """Mock RDKit descriptors module."""
    
    @staticmethod
    def ExactMolWt(mol):
        """Calculate exact molecular weight."""
        return mol.mol_weight
    
    @staticmethod
    def MolWt(mol):
        """Calculate molecular weight."""
        return mol.mol_weight
    
    @staticmethod
    def MolLogP(mol):
        """Calculate LogP."""
        return mol.logp
    
    @staticmethod
    def NumHDonors(mol):
        """Count hydrogen bond donors."""
        return mol.h_donors
    
    @staticmethod
    def NumHAcceptors(mol):
        """Count hydrogen bond acceptors."""
        return mol.h_acceptors
    
    @staticmethod
    def TPSA(mol):
        """Calculate topological polar surface area."""
        return mol.tpsa
    
    @staticmethod
    def NumRotatableBonds(mol):
        """Count rotatable bonds."""
        return mol.rotatable_bonds

class MockMolecule:
    """Mock RDKit molecule."""
    
    def __init__(
        self, 
        smiles: str = None, 
        inchi: str = None,
        mol_weight: float = 180.16,
        logp: float = -0.5,
        h_donors: int = 5,
        h_acceptors: int = 6,
        tpsa: float = 110.38,
        rotatable_bonds: int = 3,
        num_atoms: int = 12,
        num_bonds: int = 12,
        formula: str = 'C6H12O6'
    ):
        self.smiles = smiles or 'C'
        self.inchi = inchi or 'InChI=1S/CH4/h1H4'
        self.mol_weight = mol_weight
        self.logp = logp
        self.h_donors = h_donors
        self.h_acceptors = h_acceptors
        self.tpsa = tpsa
        self.rotatable_bonds = rotatable_bonds
        self.num_atoms = num_atoms
        self.num_bonds = num_bonds
        self.formula = formula
        
    def GetNumAtoms(self):
        """Get number of atoms."""
        return self.num_atoms
        
    def GetNumBonds(self):
        """Get number of bonds."""
        return self.num_bonds
        
    def GetNumHeavyAtoms(self):
        """Get number of heavy atoms."""
        return self.num_atoms - self.h_donors
        
    def GetFormula(self):
        """Get molecular formula."""
        return self.formula
        
    def GetMolWt(self):
        """Get molecular weight."""
        return self.mol_weight

class MockChem:
    """Mock RDKit Chem module."""
    
    descriptors = MockDescriptors()
    
    @staticmethod
    def MolFromSmiles(smiles, sanitize=True):
        """Create molecule from SMILES."""
        if not smiles:
            return None
            
        # Return None for invalid SMILES to simulate RDKit behavior
        invalid_smiles = ['invalid', '']
        if smiles in invalid_smiles:
            return None
            
        return MockMolecule(smiles=smiles)
        
    @staticmethod
    def MolFromInchi(inchi, sanitize=True):
        """Create molecule from InChI."""
        if not inchi:
            return None
            
        # Return None for invalid InChI to simulate RDKit behavior
        invalid_inchi = ['invalid', '']
        if inchi in invalid_inchi:
            return None
            
        return MockMolecule(inchi=inchi)
        
    @staticmethod
    def MolToSmiles(mol, isomericSmiles=True):
        """Convert molecule to SMILES."""
        if mol is None:
            return None
            
        return mol.smiles
        
    @staticmethod
    def MolToInchi(mol):
        """Convert molecule to InChI."""
        if mol is None:
            return None
            
        return mol.inchi
        
    @staticmethod
    def AddHs(mol):
        """Add hydrogens to molecule."""
        if mol is None:
            return None
            
        return mol
        
    @staticmethod
    def RemoveHs(mol):
        """Remove hydrogens from molecule."""
        if mol is None:
            return None
            
        return mol
        
    @staticmethod
    def MolFromMolBlock(molblock, sanitize=True):
        """Create molecule from MDL molblock."""
        if not molblock:
            return None
            
        return MockMolecule()
        
    @staticmethod
    def MolToMolBlock(mol):
        """Convert molecule to MDL molblock."""
        if mol is None:
            return None
            
        return "Mock Mol Block"

class MockDraw:
    """Mock RDKit Draw module."""
    
    @staticmethod
    def MolToImage(mol, size=(300, 300), kekulize=True, wedgeBonds=True, fitImage=False, options=None):
        """Convert molecule to image."""
        if mol is None:
            return None
            
        # Return a simple PIL-like image mock
        from unittest.mock import MagicMock
        image = MagicMock()
        image.size = size
        image.save = MagicMock()
        return image
        
    @staticmethod
    def MolToFile(mol, filename, size=(300, 300), kekulize=True, wedgeBonds=True, imageType=None, fitImage=False, options=None):
        """Save molecule as image file."""
        if mol is None:
            return None
            
        # Just pretend to save the file
        with open(filename, 'w') as f:
            f.write('Mock molecule image')
            
        return True

class MockRDKit:
    """Mock RDKit library."""
    
    def __init__(self):
        self.Chem = MockChem()
        self.Descriptors = MockDescriptors()
        self.Draw = MockDraw()

@contextmanager
def patch_rdkit(mock_rdkit: Optional[MockRDKit] = None) -> Generator[MockRDKit, None, None]:
    """
    Context manager to patch RDKit.
    
    Args:
        mock_rdkit: Optional mock RDKit instance to use
        
    Yields:
        The mock RDKit instance
    """
    if mock_rdkit is None:
        mock_rdkit = MockRDKit()
        
    # Modules to patch
    modules_to_patch = [
        'rdkit.Chem',
        'rdkit.Chem.Descriptors',
        'rdkit.Chem.Draw',
        'api.rdkit_utils.Chem',
        'api.rdkit_utils.Descriptors',
        'api.rdkit_utils.Draw'
    ]
    
    # Start all patches
    patches = []
    for module in modules_to_patch:
        if module.endswith('Descriptors'):
            p = patch(module, mock_rdkit.Descriptors)
        elif module.endswith('Draw'):
            p = patch(module, mock_rdkit.Draw)
        else:
            p = patch(module, mock_rdkit.Chem)
            
        patches.append(p)
        p.start()
    
    try:
        yield mock_rdkit
    finally:
        # Stop all patches
        for p in patches:
            p.stop()

# Create additional fixtures for making it easy to use the mocks in tests
@pytest.fixture
def mock_rdkit() -> Generator[MockRDKit, None, None]:
    """
    Provide a mock RDKit implementation.
    
    Yields:
        MockRDKit instance
    """
    with patch_rdkit() as rdkit:
        yield rdkit
        
@pytest.fixture
def molecule_factory() -> callable:
    """
    Provide a factory function for creating mock molecules.
    
    Returns:
        Factory function for creating MockMolecule instances
    """
    def create_molecule(
        smiles: str = 'CC',
        mol_weight: float = 180.16,
        logp: float = -0.5,
        h_donors: int = 5,
        h_acceptors: int = 6,
        formula: str = 'C6H12O6'
    ) -> MockMolecule:
        """Create a mock molecule with specified properties."""
        return MockMolecule(
            smiles=smiles,
            mol_weight=mol_weight,
            logp=logp,
            h_donors=h_donors,
            h_acceptors=h_acceptors,
            formula=formula
        )
        
    return create_molecule