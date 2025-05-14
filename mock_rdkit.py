"""
Mock implementation of RDKit for Heroku deployment.
This file provides stub implementations of RDKit classes and functions.
"""

import logging

logger = logging.getLogger(__name__)
logger.warning("Using mock RDKit implementation. Limited functionality available.")

class Chem:
    """Mock implementation of RDKit's Chem module."""
    
    @staticmethod
    def MolFromSmiles(smiles):
        """Mock implementation that returns a simple object with the SMILES string."""
        if not smiles or not isinstance(smiles, str):
            return None
        return MockMolecule(smiles)
    
    @staticmethod
    def MolToSmiles(mol, **kwargs):
        """Mock implementation that returns the SMILES string from the molecule."""
        if not mol:
            return None
        return mol.smiles
    
    @staticmethod
    def MolFromMolBlock(molblock, **kwargs):
        """Mock implementation that returns a simple object."""
        if not molblock:
            return None
        return MockMolecule("C")  # Default carbon atom
    
    @staticmethod
    def MolToMolBlock(mol, **kwargs):
        """Mock implementation that returns a simple molblock."""
        if not mol:
            return None
        return f"Mock Mol Block for {mol.smiles}"
    
    class rdMolDescriptors:
        """Mock implementation of RDKit's molecular descriptors."""
        
        @staticmethod
        def CalcExactMolWt(mol):
            """Mock implementation that returns a fixed value."""
            return 100.0
        
        @staticmethod
        def CalcNumHBD(mol):
            """Mock implementation that returns a fixed value."""
            return 2
        
        @staticmethod
        def CalcNumHBA(mol):
            """Mock implementation that returns a fixed value."""
            return 2
        
        @staticmethod
        def CalcTPSA(mol):
            """Mock implementation that returns a fixed value."""
            return 50.0
        
        @staticmethod
        def CalcFractionCSP3(mol):
            """Mock implementation that returns a fixed value."""
            return 0.5
    
    class Descriptors:
        """Mock implementation of RDKit's descriptors."""
        
        @staticmethod
        def MolWt(mol):
            """Mock implementation that returns a fixed value."""
            return 100.0
        
        @staticmethod
        def MolLogP(mol):
            """Mock implementation that returns a fixed value."""
            return 2.0
        
        @staticmethod
        def NumRotatableBonds(mol):
            """Mock implementation that returns a fixed value."""
            return 3
    
    class AllChem:
        """Mock implementation of RDKit's AllChem module."""
        
        @staticmethod
        def Compute2DCoords(mol):
            """Mock implementation that does nothing."""
            return True
        
        @staticmethod
        def EmbedMolecule(mol, **kwargs):
            """Mock implementation that does nothing."""
            return 0
        
        @staticmethod
        def MMFFOptimizeMolecule(mol, **kwargs):
            """Mock implementation that does nothing."""
            return 0

class Draw:
    """Mock implementation of RDKit's Draw module."""
    
    @staticmethod
    def MolToImage(mol, size=(300, 300), **kwargs):
        """Mock implementation that returns None."""
        return None
    
    @staticmethod
    def MolToFile(mol, filename, size=(300, 300), **kwargs):
        """Mock implementation that does nothing."""
        with open(filename, 'w') as f:
            f.write(f"Mock molecule image for {mol.smiles}")
        return True

class MockMolecule:
    """Mock implementation of RDKit's Mol class."""
    
    def __init__(self, smiles):
        self.smiles = smiles
        self.formula = "C8H10N4O2"  # Caffeine formula as default
    
    def GetNumAtoms(self):
        """Mock implementation that returns a fixed value."""
        return 14
    
    def GetNumBonds(self):
        """Mock implementation that returns a fixed value."""
        return 15
    
    def GetFormula(self):
        """Mock implementation that returns a fixed formula."""
        return self.formula

# Replace imports in other modules
import sys
sys.modules['rdkit'] = type('rdkit', (), {})
sys.modules['rdkit.Chem'] = type('Chem', (), {'Chem': Chem})
sys.modules['rdkit.Chem.AllChem'] = type('AllChem', (), {'AllChem': Chem.AllChem})
sys.modules['rdkit.Chem.Draw'] = type('Draw', (), {'Draw': Draw})
sys.modules['rdkit.Chem.rdMolDescriptors'] = type('rdMolDescriptors', (), {'rdMolDescriptors': Chem.rdMolDescriptors})
sys.modules['rdkit.Chem.Descriptors'] = type('Descriptors', (), {'Descriptors': Chem.Descriptors})