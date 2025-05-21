#!/usr/bin/env python3
"""
RDKit Wrapper - Provides a unified interface to RDKit functionality
with fallback to mock implementation when needed.

This module serves as a compatibility layer between the CryoProtect application
and RDKit, allowing seamless operation in both environments with and without
full RDKit availability.
"""

import logging
import os
import sys
from typing import Dict, Any, Optional, List, Tuple, Union

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Try to import RDKit
RDKIT_AVAILABLE = False
try:
    from rdkit import Chem, __version__ as rdkit_version
    from rdkit.Chem import Descriptors, Lipinski, MolSurf, AllChem, rdMolDescriptors
    try:
        from rdkit.Chem import Draw
        from rdkit.Chem.Draw import rdMolDraw2D
        VISUALIZATION_AVAILABLE = True
    except ImportError:
        logger.warning("RDKit visualization modules not available")
        VISUALIZATION_AVAILABLE = False
    
    RDKIT_AVAILABLE = True
    RDKIT_VERSION = rdkit_version
    logger.info(f"Using RDKit version {RDKIT_VERSION}")
except ImportError:
    logger.warning("RDKit not available, falling back to mock implementation")
    RDKIT_VERSION = None
    VISUALIZATION_AVAILABLE = False
    
    # Import mock implementation
    try:
        from mock_rdkit_formula import (
            calculate_molecular_formula,
            calculate_molecular_weight,
            calculate_logp,
            calculate_tpsa,
            calculate_h_donors,
            calculate_h_acceptors,
            calculate_rotatable_bonds,
            calculate_ring_count,
            calculate_aromatic_ring_count,
            calculate_heavy_atom_count
        )
        logger.info("Successfully imported mock_rdkit_formula for fallback calculations")
    except ImportError:
        logger.error("Failed to import mock_rdkit_formula, RDKit functionality will be limited")

# Unified interface functions that work with both real and mock RDKit

def create_molecule_from_smiles(smiles: str) -> Optional[Any]:
    """
    Create a molecule object from SMILES string.
    
    Args:
        smiles: SMILES string representation of the molecule
        
    Returns:
        RDKit molecule object or the SMILES string if using mock mode
    """
    if not smiles:
        return None
        
    if RDKIT_AVAILABLE:
        return Chem.MolFromSmiles(smiles)
    else:
        # For mock implementation, just return the SMILES as a placeholder
        # This allows functions to work with both real and mock implementations
        return smiles

def create_molecule_from_inchi(inchi: str) -> Optional[Any]:
    """
    Create a molecule object from InChI string.
    
    Args:
        inchi: InChI string representation of the molecule
        
    Returns:
        RDKit molecule object or None if using mock mode
    """
    if not inchi:
        return None
        
    if RDKIT_AVAILABLE:
        return Chem.MolFromInchi(inchi)
    else:
        logger.warning("InChI conversion not supported in mock mode")
        return None

def smiles_to_inchi(smiles: str) -> Optional[str]:
    """
    Convert SMILES to InChI.
    
    Args:
        smiles: SMILES string representation of the molecule
        
    Returns:
        InChI string or None if not available
    """
    if not smiles:
        return None
        
    if RDKIT_AVAILABLE:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            return Chem.MolToInchi(mol)
    return None

def smiles_to_inchi_key(smiles: str) -> Optional[str]:
    """
    Convert SMILES to InChI Key.
    
    Args:
        smiles: SMILES string representation of the molecule
        
    Returns:
        InChI Key string or None if not available
    """
    if not smiles:
        return None
        
    if RDKIT_AVAILABLE:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            return Chem.MolToInchiKey(mol)
    return None

def calculate_properties(mol_or_smiles: Union[str, Any]) -> Dict[str, Any]:
    """
    Calculate molecular properties for a molecule or SMILES string.
    
    Args:
        mol_or_smiles: RDKit molecule object or SMILES string
        
    Returns:
        Dictionary of calculated properties
    """
    properties = {}
    
    # Convert SMILES to molecule if needed
    if isinstance(mol_or_smiles, str):
        mol = create_molecule_from_smiles(mol_or_smiles)
        smiles = mol_or_smiles
    else:
        mol = mol_or_smiles
        smiles = Chem.MolToSmiles(mol) if RDKIT_AVAILABLE and mol is not None else str(mol)
    
    if not mol:
        return properties
    
    # Calculate properties using either real or mock RDKit
    if RDKIT_AVAILABLE:
        # Standard properties
        properties["molecular_formula"] = Chem.rdMolDescriptors.CalcMolFormula(mol)
        properties["molecular_weight"] = round(Descriptors.MolWt(mol), 2)
        properties["exact_mass"] = round(Descriptors.ExactMolWt(mol), 4)
        properties["logp"] = round(Descriptors.MolLogP(mol), 2)
        properties["tpsa"] = round(MolSurf.TPSA(mol), 2)
        properties["h_donors"] = Lipinski.NumHDonors(mol)
        properties["h_acceptors"] = Lipinski.NumHAcceptors(mol)
        properties["rotatable_bonds"] = Descriptors.NumRotatableBonds(mol)
        properties["ring_count"] = Chem.rdMolDescriptors.CalcNumRings(mol)
        properties["aromatic_ring_count"] = Chem.rdMolDescriptors.CalcNumAromaticRings(mol)
        properties["heavy_atom_count"] = mol.GetNumHeavyAtoms()
        
        # Additional sophisticated properties available in full RDKit
        properties["fraction_csp3"] = Descriptors.FractionCSP3(mol)
        properties["num_stereocenters"] = Chem.rdMolDescriptors.CalcNumAtomStereoCenters(mol)
        properties["qed"] = Descriptors.qed(mol)  # Quantitative Estimate of Drug-likeness
        
        # Add fingerprint information
        properties["morgan_fingerprint_available"] = True
    else:
        # Use mock implementation for core properties
        properties["molecular_formula"] = calculate_molecular_formula(smiles)
        properties["molecular_weight"] = calculate_molecular_weight(smiles)
        properties["logp"] = calculate_logp(smiles)
        properties["tpsa"] = calculate_tpsa(smiles)
        properties["h_donors"] = calculate_h_donors(smiles)
        properties["h_acceptors"] = calculate_h_acceptors(smiles)
        properties["rotatable_bonds"] = calculate_rotatable_bonds(smiles)
        properties["ring_count"] = calculate_ring_count(smiles)
        properties["aromatic_ring_count"] = calculate_aromatic_ring_count(smiles)
        properties["heavy_atom_count"] = calculate_heavy_atom_count(smiles)
        
        # Advanced properties not available in mock mode
        properties["morgan_fingerprint_available"] = False
    
    # Add implementation metadata
    properties["rdkit_available"] = RDKIT_AVAILABLE
    properties["rdkit_version"] = RDKIT_VERSION
    properties["calculation_method"] = "rdkit" if RDKIT_AVAILABLE else "mock_rdkit_formula"
    
    return properties

def generate_fingerprint(mol_or_smiles: Union[str, Any], radius: int = 2, n_bits: int = 2048) -> Optional[Any]:
    """
    Generate Morgan fingerprint for molecule.
    
    Args:
        mol_or_smiles: RDKit molecule object or SMILES string
        radius: Fingerprint radius
        n_bits: Number of bits in fingerprint
        
    Returns:
        RDKit fingerprint object or None if not available
    """
    if not RDKIT_AVAILABLE:
        logger.warning("Fingerprint generation not available in mock mode")
        return None
    
    # Convert SMILES to molecule if needed
    if isinstance(mol_or_smiles, str):
        mol = create_molecule_from_smiles(mol_or_smiles)
    else:
        mol = mol_or_smiles
    
    if not mol:
        return None
    
    try:
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=n_bits)
        return fp
    except Exception as e:
        logger.error(f"Error generating fingerprint: {e}")
        return None

def calculate_similarity(mol1_or_smiles1: Union[str, Any], mol2_or_smiles2: Union[str, Any]) -> Dict[str, float]:
    """
    Calculate similarity between two molecules.
    
    Args:
        mol1_or_smiles1: First molecule (RDKit mol or SMILES)
        mol2_or_smiles2: Second molecule (RDKit mol or SMILES)
        
    Returns:
        Dictionary with similarity metrics
    """
    result = {"tanimoto": 0.0}
    
    if not RDKIT_AVAILABLE:
        logger.warning("Similarity calculation not available in mock mode")
        return result
    
    # Get fingerprints
    fp1 = generate_fingerprint(mol1_or_smiles1)
    fp2 = generate_fingerprint(mol2_or_smiles2)
    
    if fp1 and fp2:
        from rdkit import DataStructs
        result["tanimoto"] = DataStructs.TanimotoSimilarity(fp1, fp2)
    
    return result

def perform_substructure_search(query: str, target: Union[str, Any], query_type: str = "smarts", target_type: str = "smiles") -> Dict[str, Any]:
    """
    Perform substructure search.
    
    Args:
        query: SMARTS or SMILES pattern
        target: Target molecule (RDKit mol or SMILES)
        query_type: Type of query ('smarts' or 'smiles')
        target_type: Type of target ('mol' or 'smiles')
        
    Returns:
        Dictionary with search results
    """
    result = {"match": False, "match_atoms": []}
    
    if not RDKIT_AVAILABLE:
        logger.warning("Substructure search not available in mock mode")
        return result
    
    # Process query
    query_mol = None
    if query_type.lower() == "smarts":
        query_mol = Chem.MolFromSmarts(query)
    elif query_type.lower() == "smiles":
        query_mol = Chem.MolFromSmiles(query)
    
    if not query_mol:
        logger.error(f"Failed to parse query {query}")
        return result
    
    # Process target
    target_mol = None
    if isinstance(target, str) and target_type.lower() == "smiles":
        target_mol = Chem.MolFromSmiles(target)
    elif not isinstance(target, str):  # Assume it's already a molecule object
        target_mol = target
    
    if not target_mol:
        logger.error(f"Failed to process target molecule")
        return result
    
    # Perform search
    match = target_mol.HasSubstructMatch(query_mol)
    result["match"] = match
    
    if match:
        # Get matching atoms
        match_atoms = target_mol.GetSubstructMatch(query_mol)
        result["match_atoms"] = list(match_atoms)
    
    return result

def generate_molecule_svg(mol_or_smiles: Union[str, Any], width: int = 300, height: int = 200, 
                          highlight_atoms: List[int] = None, highlight_bonds: List[int] = None) -> Optional[str]:
    """
    Generate SVG visualization of a molecule.
    
    Args:
        mol_or_smiles: RDKit molecule object or SMILES string
        width: Image width in pixels
        height: Image height in pixels
        highlight_atoms: List of atom indices to highlight
        highlight_bonds: List of bond indices to highlight
        
    Returns:
        SVG string or None if visualization not available
    """
    if not RDKIT_AVAILABLE or not VISUALIZATION_AVAILABLE:
        logger.warning("RDKit visualization not available")
        return None
    
    # Convert SMILES to molecule if needed
    if isinstance(mol_or_smiles, str):
        mol = create_molecule_from_smiles(mol_or_smiles)
    else:
        mol = mol_or_smiles
    
    if not mol:
        return None
    
    # Generate SVG
    try:
        drawer = rdMolDraw2D.MolDraw2DSVG(width, height)
        
        # Apply highlighting if specified
        if highlight_atoms or highlight_bonds:
            highlight_atm = highlight_atoms or []
            highlight_bnd = highlight_bonds or []
            drawer.DrawMolecule(mol, highlightAtoms=highlight_atm, highlightBonds=highlight_bnd)
        else:
            drawer.DrawMolecule(mol)
            
        drawer.FinishDrawing()
        return drawer.GetDrawingText()
    except Exception as e:
        logger.error(f"Error generating SVG: {e}")
        return None

def generate_molecule_3d_coordinates(mol_or_smiles: Union[str, Any], optimize: bool = True) -> Optional[Any]:
    """
    Generate 3D coordinates for a molecule.
    
    Args:
        mol_or_smiles: RDKit molecule object or SMILES string
        optimize: Whether to optimize the structure after embedding
        
    Returns:
        RDKit molecule with 3D coordinates or None if not available
    """
    if not RDKIT_AVAILABLE:
        logger.warning("3D coordinate generation not available in mock mode")
        return None
    
    # Convert SMILES to molecule if needed
    if isinstance(mol_or_smiles, str):
        mol = create_molecule_from_smiles(mol_or_smiles)
    else:
        mol = mol_or_smiles
    
    if not mol:
        return None
    
    try:
        # Add hydrogens (important for proper 3D structure)
        mol = Chem.AddHs(mol)
        
        # Generate 3D coordinates
        AllChem.EmbedMolecule(mol, randomSeed=42)
        
        # Optimize if requested
        if optimize:
            AllChem.UFFOptimizeMolecule(mol)
        
        return mol
    except Exception as e:
        logger.error(f"Error generating 3D coordinates: {e}")
        return None

def get_rdkit_status() -> Dict[str, Any]:
    """
    Get the status of RDKit integration.
    
    Returns:
        Dictionary with RDKit status information
    """
    status = {
        "rdkit_available": RDKIT_AVAILABLE,
        "rdkit_version": RDKIT_VERSION,
        "visualization_available": VISUALIZATION_AVAILABLE,
        "using_mock_implementation": not RDKIT_AVAILABLE,
        "mock_capabilities": [
            "molecular_formula",
            "molecular_weight",
            "logp",
            "tpsa",
            "h_donors",
            "h_acceptors",
            "rotatable_bonds",
            "ring_count",
            "aromatic_ring_count",
            "heavy_atom_count"
        ]
    }
    
    # Add Python information
    status["python_version"] = sys.version
    status["python_executable"] = sys.executable
    
    return status

# Function to test RDKit wrapper capabilities
def run_test() -> Dict[str, Any]:
    """
    Run a simple test of RDKit wrapper capabilities.
    
    Returns:
        Dictionary with test results
    """
    test_results = {
        "status": get_rdkit_status(),
        "tests": {}
    }
    
    # Test molecules
    test_molecules = [
        ("Ethanol", "CCO"),
        ("Aspirin", "CC(=O)Oc1ccccc1C(=O)O"),
        ("Caffeine", "CN1C=NC2=C1C(=O)N(C(=O)N2C)C")
    ]
    
    # Test functions
    for name, smiles in test_molecules:
        # Create molecule
        mol = create_molecule_from_smiles(smiles)
        test_results["tests"][name] = {
            "smiles": smiles,
            "molecule_created": mol is not None,
            "properties": calculate_properties(smiles)
        }
        
        # Test visualization if available
        if VISUALIZATION_AVAILABLE:
            test_results["tests"][name]["svg_generated"] = generate_molecule_svg(smiles) is not None
        
        # Test 3D coordinates if available
        if RDKIT_AVAILABLE:
            test_results["tests"][name]["3d_coords_generated"] = generate_molecule_3d_coordinates(smiles) is not None
    
    return test_results

if __name__ == "__main__":
    # Run standalone test if executed directly
    import json
    
    print("Testing RDKit Wrapper capabilities...")
    results = run_test()
    
    # Print results
    print(f"RDKit Available: {results['status']['rdkit_available']}")
    if results['status']['rdkit_available']:
        print(f"RDKit Version: {results['status']['rdkit_version']}")
    else:
        print("Using Mock Implementation")
    
    print("\nTest Results:")
    for molecule, result in results["tests"].items():
        print(f"\n{molecule}:")
        print(f"  SMILES: {result['smiles']}")
        print(f"  Molecule Created: {result['molecule_created']}")
        print(f"  Properties:")
        for prop, value in result['properties'].items():
            if prop not in ['rdkit_available', 'rdkit_version', 'calculation_method']:
                print(f"    {prop}: {value}")
    
    # Save results to file
    with open("rdkit_wrapper_test_results.json", "w") as f:
        json.dump(results, f, indent=2)
    
    print("\nTest results saved to rdkit_wrapper_test_results.json")