#!/usr/bin/env python3
"""
CryoProtect Analyzer - RDKit Verification Script

This script verifies that the RDKit integration works correctly by testing
key functionality and ensuring that the Docker container works properly.
"""

import os
import sys
import logging
import json
import argparse
import random
import datetime
from typing import Dict, Any, List, Optional

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

def check_rdkit_import() -> bool:
    """Check if RDKit can be imported."""
    logger.info("Checking RDKit import...")
    try:
        from rdkit import Chem
        from rdkit.Chem import Descriptors, Lipinski, MolSurf, AllChem, Fragments
        from rdkit.Chem.Scaffolds import MurckoScaffold
        from rdkit.Chem import Draw
        from rdkit.Chem.Draw import rdMolDraw2D
        
        # Try to import IPythonConsole, but don't fail if it's not available
        try:
            from rdkit.Chem.Draw import IPythonConsole
            logger.info("IPython is available for RDKit visualization")
        except ImportError:
            logger.warning("IPython is not available. Some visualization features may be limited.")
        
        logger.info("RDKit import successful")
        return True
    except ImportError as e:
        logger.error(f"RDKit import failed: {str(e)}")
        return False

def test_molecule_parsing() -> bool:
    """Test molecule parsing functionality."""
    logger.info("Testing molecule parsing...")
    try:
        from rdkit import Chem
        
        # Test valid SMILES
        ethanol_smiles = "CCO"
        mol = Chem.MolFromSmiles(ethanol_smiles)
        if mol is None:
            logger.error("Failed to parse valid SMILES")
            return False
        
        logger.info("Molecule parsing successful")
        return True
    except Exception as e:
        logger.error(f"Error testing molecule parsing: {str(e)}")
        return False

def test_property_calculation() -> bool:
    """Test property calculation functionality."""
    logger.info("Testing property calculation...")
    try:
        from api.rdkit_utils import calculate_all_properties
        
        # Test with ethanol
        ethanol_smiles = "CCO"
        properties = calculate_all_properties(ethanol_smiles)
        
        # Check that we have all the property categories
        required_properties = [
            "hydrogen_bonding", "logp", "tpsa", "molecular_properties",
            "functional_groups", "permeability", "smiles", "inchi", "inchi_key"
        ]
        
        for prop in required_properties:
            if prop not in properties:
                logger.error(f"Missing property: {prop}")
                return False
        
        logger.info("Property calculation successful")
        return True
    except Exception as e:
        logger.error(f"Error testing property calculation: {str(e)}")
        return False

def test_visualization(reference_molecules: bool = False) -> bool:
    """Test visualization functionality."""
    logger.info("Testing visualization...")
    try:
        from api.rdkit_utils import generate_molecule_svg
        
        # Test molecules
        if reference_molecules:
            logger.info("Testing visualization with reference molecules")
            # Common cryoprotectants
            test_molecules = [
                ("DMSO", "CS(=O)C"),
                ("Glycerol", "C(C(CO)O)O"),
                ("Ethylene glycol", "C(CO)O"),
                ("Propylene glycol", "CC(CO)O"),
                ("Trehalose", "C1C(C(C(C(C1O)OC2C(C(C(C(O2)CO)O)O)O)O)O)O")
            ]
            
            success_count = 0
            for name, smiles in test_molecules:
                logger.info(f"Testing visualization of {name} ({smiles})")
                svg = generate_molecule_svg(smiles)
                
                # Check that we got a valid SVG string
                if isinstance(svg, str) and svg and "<svg" in svg and "</svg>" in svg:
                    logger.info(f"Successfully generated SVG for {name}")
                    success_count += 1
                else:
                    logger.error(f"Failed to generate SVG for {name}")
            
            success_rate = success_count / len(test_molecules) * 100
            logger.info(f"Visualization success rate: {success_rate:.1f}% ({success_count}/{len(test_molecules)})")
            return success_count == len(test_molecules)
        else:
            # Test with ethanol
            ethanol_smiles = "CCO"
            svg = generate_molecule_svg(ethanol_smiles)
            
            # Check that we got an SVG string
            if not isinstance(svg, str) or not svg:
                logger.error("Failed to generate SVG")
                return False
            
            if "<svg" not in svg or "</svg>" not in svg:
                logger.error("Generated string is not a valid SVG")
                return False
            
            logger.info("Visualization successful")
            return True
    except Exception as e:
        logger.error(f"Error testing visualization: {str(e)}")
        return False

def test_substructure_search() -> bool:
    """Test substructure search functionality."""
    logger.info("Testing substructure search...")
    try:
        from api.rdkit_utils import perform_substructure_search
        
        # Search for alcohol group in ethanol
        ethanol_smiles = "CCO"
        result = perform_substructure_search("[OH]", ethanol_smiles, "smarts", "smiles")
        
        # Should find a match
        if not result.get("match", False):
            logger.error("Failed to find substructure match")
            return False
        
        logger.info("Substructure search successful")
        return True
    except Exception as e:
        logger.error(f"Error testing substructure search: {str(e)}")
        return False

def test_similarity_calculation() -> bool:
    """Test similarity calculation functionality."""
    logger.info("Testing similarity calculation...")
    try:
        from api.rdkit_utils import calculate_similarity
        
        # Compare ethanol to itself (should be identical)
        ethanol_smiles = "CCO"
        result = calculate_similarity(ethanol_smiles, ethanol_smiles)
        
        # Tanimoto similarity should be 1.0 for identical molecules
        if result.get("tanimoto", 0.0) != 1.0:
            logger.error(f"Unexpected similarity value: {result.get('tanimoto', 0.0)}")
            return False
        
        logger.info("Similarity calculation successful")
        return True
    except Exception as e:
        logger.error(f"Error testing similarity calculation: {str(e)}")
        return False

def test_3d_coords() -> bool:
    """Test 3D coordinate generation functionality."""
    logger.info("Testing 3D coordinate generation...")
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
        
        # Test molecules
        test_molecules = [
            ("DMSO", "CS(=O)C"),
            ("Glycerol", "C(C(CO)O)O"),
            ("Ethylene glycol", "C(CO)O"),
            ("Propylene glycol", "CC(CO)O"),
            ("Trehalose", "C1C(C(C(C(C1O)OC2C(C(C(C(O2)CO)O)O)O)O)O)O")
        ]
        
        success_count = 0
        for name, smiles in test_molecules:
            logger.info(f"Testing 3D coordinate generation for {name} ({smiles})")
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                logger.error(f"Failed to parse SMILES for {name}")
                continue
                
            # Add hydrogens
            mol = Chem.AddHs(mol)
            
            # Generate 3D coordinates
            success = AllChem.EmbedMolecule(mol, randomSeed=42)
            if success == 0:  # 0 means success in RDKit
                # Optimize the structure
                success = AllChem.UFFOptimizeMolecule(mol, maxIters=200)
                if success == 0:
                    logger.info(f"Successfully generated 3D coordinates for {name}")
                    success_count += 1
                else:
                    logger.error(f"Failed to optimize 3D structure for {name}")
            else:
                logger.error(f"Failed to generate 3D coordinates for {name}")
        
        success_rate = success_count / len(test_molecules) * 100
        logger.info(f"3D coordinate generation success rate: {success_rate:.1f}% ({success_count}/{len(test_molecules)})")
        return success_count == len(test_molecules)
    except Exception as e:
        logger.error(f"Error testing 3D coordinate generation: {str(e)}")
        return False

def test_random_samples(num_samples: int = 50) -> bool:
    """Test visualization with random samples from the database."""
    logger.info(f"Testing visualization with {num_samples} random samples from the database...")
    try:
        import psycopg2
        from api.rdkit_utils import generate_molecule_svg
        from config import get_db_config
        
        # Connect to the database
        db_config = get_db_config()
        conn = psycopg2.connect(
            host=db_config.get('host', 'localhost'),
            port=db_config.get('port', 5432),
            dbname=db_config.get('dbname', 'postgres'),
            user=db_config.get('user', 'postgres'),
            password=db_config.get('password', '')
        )
        
        # Get random molecules from the database
        cursor = conn.cursor()
        cursor.execute(f"SELECT id, name, smiles FROM molecules ORDER BY RANDOM() LIMIT {num_samples}")
        molecules = cursor.fetchall()
        
        if not molecules:
            logger.error("No molecules found in the database")
            return False
        
        # Test visualization for each molecule
        success_count = 0
        failed_molecules = []
        
        for mol_id, name, smiles in molecules:
            if not smiles:
                logger.warning(f"Molecule {mol_id} ({name}) has no SMILES string")
                continue
                
            logger.info(f"Testing visualization of {name} (ID: {mol_id})")
            try:
                svg = generate_molecule_svg(smiles)
                
                # Check that we got a valid SVG string
                if isinstance(svg, str) and svg and "<svg" in svg and "</svg>" in svg:
                    logger.info(f"Successfully generated SVG for {name}")
                    success_count += 1
                else:
                    logger.error(f"Failed to generate SVG for {name}")
                    failed_molecules.append((mol_id, name, smiles))
            except Exception as e:
                logger.error(f"Error generating SVG for {name}: {str(e)}")
                failed_molecules.append((mol_id, name, smiles))
        
        # Close database connection
        cursor.close()
        conn.close()
        
        # Calculate success rate
        total_tested = len(molecules)
        success_rate = (success_count / total_tested * 100) if total_tested > 0 else 0
        logger.info(f"Visualization success rate: {success_rate:.1f}% ({success_count}/{total_tested})")
        
        # Generate report
        report = {
            "timestamp": datetime.datetime.now().isoformat(),
            "test_type": "random_samples_visualization",
            "total_samples": total_tested,
            "successful": success_count,
            "failed": total_tested - success_count,
            "success_rate": success_rate,
            "failed_molecules": [
                {"id": mol_id, "name": name, "smiles": smiles}
                for mol_id, name, smiles in failed_molecules
            ]
        }
        
        # Save report to file
        os.makedirs("reports", exist_ok=True)
        report_path = "reports/visualization_validation_report.json"
        with open(report_path, 'w') as f:
            json.dump(report, f, indent=2)
        
        logger.info(f"Saved visualization validation report to {report_path}")
        
        return success_rate >= 90  # Consider test successful if at least 90% of molecules can be visualized
    except Exception as e:
        logger.error(f"Error testing random samples: {str(e)}")
        return False

def run_all_tests() -> Dict[str, bool]:
    """Run all tests and return results."""
    results = {
        "rdkit_import": check_rdkit_import(),
        "molecule_parsing": test_molecule_parsing(),
        "property_calculation": test_property_calculation(),
        "visualization": test_visualization(),
        "substructure_search": test_substructure_search(),
        "similarity_calculation": test_similarity_calculation()
    }
    
    return results

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="Verify RDKit integration")
    parser.add_argument("--test-visualization", action="store_true", help="Test visualization with reference molecules")
    parser.add_argument("--test-3d-coords", action="store_true", help="Test 3D coordinate generation")
    parser.add_argument("--test-random-samples", type=int, metavar="N", help="Test visualization with N random samples from the database")
    parser.add_argument("--all", action="store_true", help="Run all tests")
    return parser.parse_args()

def main():
    """Main function."""
    logger.info("Starting RDKit verification...")
    
    # Add the parent directory to the path so we can import the api package
    sys.path.insert(0, os.path.abspath(os.path.dirname(__file__)))
    
    # Parse command line arguments
    args = parse_arguments()
    
    # Determine which tests to run
    if args.all or (not args.test_visualization and not args.test_3d_coords and args.test_random_samples is None):
        # Run all basic tests
        results = run_all_tests()
        
        # Print results
        logger.info("Verification results:")
        all_passed = True
        for test_name, passed in results.items():
            status = "PASSED" if passed else "FAILED"
            logger.info(f"  {test_name}: {status}")
            if not passed:
                all_passed = False
        
        # Print overall result
        if all_passed:
            logger.info("All tests passed! RDKit integration is working correctly.")
        else:
            logger.error("Some tests failed. Please check the logs for details.")
            sys.exit(1)
    else:
        # Run specific tests based on arguments
        all_passed = True
        
        if args.test_visualization:
            logger.info("Running visualization test with reference molecules...")
            vis_passed = test_visualization(reference_molecules=True)
            status = "PASSED" if vis_passed else "FAILED"
            logger.info(f"Visualization test: {status}")
            all_passed = all_passed and vis_passed
        
        if args.test_3d_coords:
            logger.info("Running 3D coordinate generation test...")
            coords_passed = test_3d_coords()
            status = "PASSED" if coords_passed else "FAILED"
            logger.info(f"3D coordinate generation test: {status}")
            all_passed = all_passed and coords_passed
        
        if args.test_random_samples is not None:
            logger.info(f"Running visualization test with {args.test_random_samples} random samples...")
            samples_passed = test_random_samples(args.test_random_samples)
            status = "PASSED" if samples_passed else "FAILED"
            logger.info(f"Random samples visualization test: {status}")
            all_passed = all_passed and samples_passed
        
        # Print overall result
        if all_passed:
            logger.info("All requested tests passed! RDKit integration is working correctly.")
        else:
            logger.error("Some tests failed. Please check the logs for details.")
            sys.exit(1)

if __name__ == "__main__":
    main()