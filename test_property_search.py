#!/usr/bin/env python3
"""
Test script to verify property-based molecule searching with RDKit
"""

import sys
import logging
from typing import List, Dict, Any

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

def create_test_molecules() -> List[Dict[str, Any]]:
    """Create a set of test molecules with properties for searching"""
    
    from rdkit import Chem
    from rdkit.Chem import Descriptors, Lipinski
    
    # Test molecule SMILES strings representing common cryoprotectants
    test_smiles = [
        ("Ethanol", "CCO"),
        ("Glycerol", "C(C(CO)O)O"),
        ("DMSO", "CS(=O)C"),
        ("Ethylene glycol", "C(CO)O"),
        ("Propylene glycol", "CC(CO)O"),
        ("Glucose", "C(C1C(C(C(C(O1)O)O)O)O)O"),
        ("Methanol", "CO"),
        ("1,3-Propanediol", "C(CCO)O"),
        ("Formamide", "C(=O)N"),
        ("Trehalose", "C1C(C(C(C(C1O)OC2C(C(C(C(O2)CO)O)O)O)O)O)O")
    ]
    
    # Process molecules and calculate properties
    molecules = []
    
    for name, smiles in test_smiles:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            logger.warning(f"Failed to parse SMILES for {name}: {smiles}")
            continue
            
        # Calculate properties
        mw = Descriptors.MolWt(mol)
        logp = Descriptors.MolLogP(mol)
        tpsa = Descriptors.TPSA(mol)
        h_donors = Lipinski.NumHDonors(mol)
        h_acceptors = Lipinski.NumHAcceptors(mol)
        rotatable_bonds = Descriptors.NumRotatableBonds(mol)
        h_bond_count = h_donors + h_acceptors
        
        molecules.append({
            "name": name,
            "smiles": smiles,
            "mol": mol,
            "molecular_weight": mw,
            "logp": logp,
            "tpsa": tpsa,
            "h_donors": h_donors,
            "h_acceptors": h_acceptors,
            "h_bond_count": h_bond_count,
            "rotatable_bonds": rotatable_bonds
        })
    
    logger.info(f"Created {len(molecules)} test molecules")
    return molecules

def search_by_property(molecules: List[Dict[str, Any]], 
                       property_name: str, 
                       min_value: float = None, 
                       max_value: float = None) -> List[Dict[str, Any]]:
    """Search molecules by property within a range"""
    
    results = []
    
    for mol in molecules:
        if property_name not in mol:
            continue
            
        value = mol[property_name]
        
        # Check if value is within range
        if min_value is not None and value < min_value:
            continue
        if max_value is not None and value > max_value:
            continue
            
        results.append(mol)
    
    return results

def substructure_search(molecules: List[Dict[str, Any]], 
                       query_smarts: str) -> List[Dict[str, Any]]:
    """Search molecules by substructure using SMARTS pattern"""
    
    from rdkit import Chem
    
    results = []
    
    # Parse query pattern
    query_mol = Chem.MolFromSmarts(query_smarts)
    if query_mol is None:
        logger.error(f"Failed to parse SMARTS pattern: {query_smarts}")
        return []
    
    for mol_data in molecules:
        if mol_data["mol"].HasSubstructMatch(query_mol):
            results.append(mol_data)
    
    return results

def run_searches(molecules: List[Dict[str, Any]]):
    """Run various property-based and substructure searches"""
    
    # Search 1: Small molecules (MW < 100)
    small_molecules = search_by_property(molecules, "molecular_weight", max_value=100)
    logger.info(f"Found {len(small_molecules)} small molecules (MW < 100):")
    for mol in small_molecules:
        logger.info(f"  {mol['name']}: {mol['molecular_weight']:.2f}")
    
    # Search 2: Hydrophilic molecules (LogP < 0)
    hydrophilic = search_by_property(molecules, "logp", max_value=0)
    logger.info(f"Found {len(hydrophilic)} hydrophilic molecules (LogP < 0):")
    for mol in hydrophilic:
        logger.info(f"  {mol['name']}: {mol['logp']:.2f}")
    
    # Search 3: Many hydrogen bonds (h_bond_count > 5)
    h_bond_rich = search_by_property(molecules, "h_bond_count", min_value=5)
    logger.info(f"Found {len(h_bond_rich)} molecules with many H-bonds (> 5):")
    for mol in h_bond_rich:
        logger.info(f"  {mol['name']}: {mol['h_bond_count']}")
    
    # Search 4: Find molecules with alcohols using SMARTS
    alcohols = substructure_search(molecules, "[OX2H]")
    logger.info(f"Found {len(alcohols)} molecules with alcohol groups:")
    for mol in alcohols:
        logger.info(f"  {mol['name']}")
    
    # Search 5: Find molecules with sulfoxide group (like DMSO)
    sulfoxides = substructure_search(molecules, "[#16X3](=[OX1])")
    logger.info(f"Found {len(sulfoxides)} molecules with sulfoxide groups:")
    for mol in sulfoxides:
        logger.info(f"  {mol['name']}")

def main():
    logger.info("Testing property-based molecular searching with RDKit...")
    
    try:
        import rdkit
        logger.info(f"Using RDKit version: {rdkit.__version__}")
        
        # Create test molecules
        molecules = create_test_molecules()
        
        # Run searches
        run_searches(molecules)
        
        logger.info("Property-based searching test completed successfully")
        return True
        
    except ImportError as e:
        logger.error(f"RDKit import failed: {str(e)}")
        return False
    except Exception as e:
        logger.error(f"Error during property-based searching test: {str(e)}")
        return False

if __name__ == "__main__":
    success = main()
    if success:
        logger.info("Property-based searching test passed!")
        sys.exit(0)
    else:
        logger.error("Property-based searching test failed")
        sys.exit(1)