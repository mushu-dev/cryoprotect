#!/usr/bin/env python3
"""
Standalone RDKit test script for CryoProtect
This script tests core RDKit functionality without dependencies on other modules
"""

import sys
import logging
from typing import Dict, Any, List, Optional

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

def test_rdkit_import() -> bool:
    """Test importing RDKit modules"""
    logger.info("Testing RDKit import...")
    
    try:
        from rdkit import Chem
        from rdkit.Chem import Descriptors, Lipinski, AllChem, MolSurf
        
        try:
            # Try to import visualization modules
            from rdkit.Chem import Draw
            from rdkit.Chem.Draw import rdMolDraw2D
            logger.info("Visualization modules available")
        except ImportError as e:
            logger.warning(f"Visualization modules not available: {e}")
        
        try:
            # Try to import IPythonConsole module
            from rdkit.Chem.Draw import IPythonConsole
            logger.info("IPython is available for RDKit visualization")
        except ImportError:
            logger.warning("IPython is not available. Some visualization features may be limited.")
        
        logger.info("RDKit import successful")
        return True
    except ImportError as e:
        logger.error(f"RDKit import failed: {e}")
        return False

def create_test_molecules() -> List[Dict[str, Any]]:
    """Create test molecules for RDKit functionality tests"""
    from rdkit import Chem
    from rdkit.Chem import Descriptors, Lipinski
    
    # Test molecule SMILES
    test_smiles = [
        ("Ethanol", "CCO"),
        ("Glycerol", "C(C(CO)O)O"),
        ("DMSO", "CS(=O)C"),
        ("Acetaminophen", "CC(=O)NC1=CC=C(C=C1)O"),
        ("Caffeine", "CN1C=NC2=C1C(=O)N(C(=O)N2C)C")
    ]
    
    molecules = []
    
    for name, smiles in test_smiles:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            logger.warning(f"Failed to parse SMILES for {name}: {smiles}")
            continue
            
        # Store molecule data
        molecules.append({
            "name": name,
            "smiles": smiles,
            "mol": mol
        })
    
    logger.info(f"Created {len(molecules)} test molecules")
    return molecules

def test_property_calculation(molecules: List[Dict[str, Any]]) -> bool:
    """Test property calculation using RDKit"""
    logger.info("Testing property calculation...")
    
    try:
        from rdkit.Chem import Descriptors, Lipinski, MolSurf
        
        for mol_data in molecules:
            name = mol_data["name"]
            mol = mol_data["mol"]
            
            # Calculate properties
            mw = Descriptors.MolWt(mol)
            exact_mw = Descriptors.ExactMolWt(mol)
            logp = Descriptors.MolLogP(mol)
            tpsa = Descriptors.TPSA(mol)
            h_donors = Lipinski.NumHDonors(mol)
            h_acceptors = Lipinski.NumHAcceptors(mol)
            rotatable_bonds = Descriptors.NumRotatableBonds(mol)
            
            logger.info(f"Properties for {name}:")
            logger.info(f"  Molecular Weight: {mw:.2f}")
            logger.info(f"  Exact Mass: {exact_mw:.4f}")
            logger.info(f"  LogP: {logp:.2f}")
            logger.info(f"  TPSA: {tpsa:.2f}")
            logger.info(f"  H-Bond Donors: {h_donors}")
            logger.info(f"  H-Bond Acceptors: {h_acceptors}")
            logger.info(f"  Rotatable Bonds: {rotatable_bonds}")
        
        logger.info("Property calculation successful")
        return True
    except Exception as e:
        logger.error(f"Error testing property calculation: {e}")
        return False

def test_3d_generation(molecules: List[Dict[str, Any]]) -> bool:
    """Test 3D coordinate generation using RDKit"""
    logger.info("Testing 3D coordinate generation...")
    
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
        
        success_count = 0
        
        for mol_data in molecules:
            name = mol_data["name"]
            mol = mol_data["mol"]
            
            # Add hydrogens
            mol_h = Chem.AddHs(mol)
            
            # Generate 3D coordinates
            logger.info(f"Generating 3D coordinates for {name}...")
            result = AllChem.EmbedMolecule(mol_h)
            
            if result == 0:  # 0 means success in RDKit
                # Optimize the structure
                result = AllChem.UFFOptimizeMolecule(mol_h)
                if result == 0:
                    logger.info(f"Successfully generated and optimized 3D structure for {name}")
                    success_count += 1
                else:
                    logger.warning(f"Failed to optimize 3D structure for {name}")
            else:
                logger.warning(f"Failed to generate 3D coordinates for {name}")
        
        success_rate = (success_count / len(molecules)) * 100
        logger.info(f"3D coordinate generation success rate: {success_rate:.1f}% ({success_count}/{len(molecules)})")
        
        return success_count > 0
    except Exception as e:
        logger.error(f"Error testing 3D coordinate generation: {e}")
        return False

def test_molecular_fingerprints(molecules: List[Dict[str, Any]]) -> bool:
    """Test molecular fingerprint generation and similarity calculation"""
    logger.info("Testing molecular fingerprints and similarity...")
    
    try:
        from rdkit import Chem, DataStructs
        from rdkit.Chem import AllChem
        
        # Generate fingerprints for all molecules
        fingerprints = []
        
        for mol_data in molecules:
            mol = mol_data["mol"]
            fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
            fingerprints.append((mol_data["name"], fp))
        
        # Calculate similarity matrix
        logger.info("Similarity matrix (Tanimoto):")
        print("              ", end="")
        for name, _ in fingerprints:
            print(f"{name:12}", end="")
        print()
        
        for i, (name1, fp1) in enumerate(fingerprints):
            print(f"{name1:12}", end="")
            for j, (name2, fp2) in enumerate(fingerprints):
                sim = DataStructs.TanimotoSimilarity(fp1, fp2)
                print(f"{sim:12.2f}", end="")
            print()
        
        logger.info("Fingerprint generation and similarity calculation successful")
        return True
    except Exception as e:
        logger.error(f"Error testing molecular fingerprints: {e}")
        return False

def test_substructure_search(molecules: List[Dict[str, Any]]) -> bool:
    """Test substructure searching using RDKit"""
    logger.info("Testing substructure search...")
    
    try:
        from rdkit import Chem
        
        # Define some SMARTS patterns to search for
        patterns = [
            ("Hydroxyl", "[OH]"),
            ("Carbonyl", "[CX3]=[OX1]"),
            ("Amide", "[NX3][CX3](=[OX1])[#6]"),
            ("Aromatic", "a1aaaaa1")  # Six-membered aromatic ring
        ]
        
        # Test each pattern against all molecules
        for pattern_name, smarts in patterns:
            logger.info(f"Searching for {pattern_name} ({smarts}):")
            
            # Parse the SMARTS pattern
            query_mol = Chem.MolFromSmarts(smarts)
            if query_mol is None:
                logger.error(f"Failed to parse SMARTS pattern: {smarts}")
                continue
            
            # Search each molecule
            matches = []
            for mol_data in molecules:
                name = mol_data["name"]
                mol = mol_data["mol"]
                
                if mol.HasSubstructMatch(query_mol):
                    matches.append(name)
            
            if matches:
                logger.info(f"  Found in: {', '.join(matches)}")
            else:
                logger.info("  Not found in any test molecule")
        
        logger.info("Substructure search testing completed")
        return True
    except Exception as e:
        logger.error(f"Error testing substructure search: {e}")
        return False

def test_visualization(molecules: List[Dict[str, Any]]) -> bool:
    """Test molecular visualization using RDKit"""
    logger.info("Testing molecular visualization...")
    
    try:
        from rdkit.Chem import Draw
        from rdkit.Chem.Draw import rdMolDraw2D
        
        for mol_data in molecules:
            name = mol_data["name"]
            mol = mol_data["mol"]
            
            # Generate SVG
            drawer = rdMolDraw2D.MolDraw2DSVG(300, 200)
            drawer.DrawMolecule(mol)
            drawer.FinishDrawing()
            svg = drawer.GetDrawingText()
            
            # Verify SVG generation
            if svg and "<svg" in svg and "</svg>" in svg:
                logger.info(f"Successfully generated SVG for {name}")
            else:
                logger.error(f"Failed to generate valid SVG for {name}")
                return False
        
        logger.info("Visualization test successful")
        return True
    except ImportError:
        logger.warning("Visualization modules not available, skipping test")
        return True
    except Exception as e:
        logger.error(f"Error testing visualization: {e}")
        return False

def main():
    """Main function to run all tests"""
    logger.info("Running standalone RDKit tests...")
    
    # Test RDKit import
    if not test_rdkit_import():
        logger.error("RDKit import test failed, aborting remaining tests")
        return 1
    
    # Get RDKit version
    try:
        from rdkit import __version__
        logger.info(f"RDKit version: {__version__}")
    except (ImportError, AttributeError):
        logger.warning("Could not determine RDKit version")
    
    # Create test molecules
    molecules = create_test_molecules()
    if not molecules:
        logger.error("Failed to create test molecules, aborting remaining tests")
        return 1
    
    # Run all tests
    tests = [
        ("Property calculation", lambda: test_property_calculation(molecules)),
        ("3D coordinate generation", lambda: test_3d_generation(molecules)),
        ("Molecular fingerprints", lambda: test_molecular_fingerprints(molecules)),
        ("Substructure search", lambda: test_substructure_search(molecules)),
        ("Visualization", lambda: test_visualization(molecules))
    ]
    
    success = True
    for test_name, test_func in tests:
        logger.info(f"\n{'=' * 40}")
        logger.info(f"Running test: {test_name}")
        logger.info(f"{'=' * 40}")
        
        test_result = test_func()
        if test_result:
            logger.info(f"{test_name} test: PASSED")
        else:
            logger.error(f"{test_name} test: FAILED")
            success = False
    
    # Print summary
    logger.info("\nTest Summary:")
    if success:
        logger.info("All RDKit tests passed successfully!")
    else:
        logger.error("Some RDKit tests failed. See logs for details.")
    
    return 0 if success else 1

if __name__ == "__main__":
    sys.exit(main())