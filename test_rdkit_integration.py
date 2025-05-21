#!/usr/bin/env python3
"""
CryoProtect v2 - RDKit Integration Tests

This script tests the RDKit integration of the CryoProtect v2 application.
It verifies that RDKit is correctly installed and functioning properly.
"""

import os
import sys
import json
import logging
from datetime import datetime
from typing import Dict, List, Any, Optional, Tuple

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler("rdkit_integration_test.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

class RDKitIntegrationTest:
    """Test class for RDKit integration."""

    def __init__(self):
        """Initialize the test class."""
        self.test_results = {
            "status": "Not Started",
            "total_tests": 0,
            "passed_tests": 0,
            "failed_tests": 0,
            "skipped_tests": 0,
            "test_cases": []
        }
        self.test_molecules = {
            "ethanol": "CCO",
            "glycerol": "C(C(CO)O)O",
            "dmso": "CS(=O)C",
            "aspirin": "CC(=O)OC1=CC=CC=C1C(=O)O"
        }

    def run_tests(self) -> Dict[str, Any]:
        """Run all RDKit integration tests."""
        self.test_results["status"] = "Running"
        self.test_results["start_time"] = datetime.now().isoformat()

        # Run the tests
        self.test_rdkit_import()
        self.test_molecule_parsing()
        self.test_property_calculation()
        self.test_visualization()
        self.test_substructure_search()
        self.test_similarity_calculation()
        self.test_scientific_accuracy()

        # Calculate test results
        self.test_results["total_tests"] = len(self.test_results["test_cases"])
        self.test_results["passed_tests"] = sum(1 for tc in self.test_results["test_cases"] if tc["status"] == "Passed")
        self.test_results["failed_tests"] = sum(1 for tc in self.test_results["test_cases"] if tc["status"] == "Failed")
        self.test_results["skipped_tests"] = sum(1 for tc in self.test_results["test_cases"] if tc["status"] == "Skipped")
        
        if self.test_results["failed_tests"] == 0:
            self.test_results["status"] = "Passed"
        else:
            self.test_results["status"] = "Failed"
        
        self.test_results["end_time"] = datetime.now().isoformat()
        
        return self.test_results

    def add_test_result(self, test_id: str, test_name: str, status: str, message: str) -> None:
        """Add a test result to the test results."""
        self.test_results["test_cases"].append({
            "id": test_id,
            "name": test_name,
            "status": status,
            "message": message
        })
        
        if status == "Passed":
            logger.info(f"Test {test_id} - {test_name}: PASSED")
        elif status == "Failed":
            logger.error(f"Test {test_id} - {test_name}: FAILED - {message}")
        else:
            logger.warning(f"Test {test_id} - {test_name}: SKIPPED - {message}")

    def test_rdkit_import(self) -> None:
        """Test that RDKit can be imported."""
        try:
            from rdkit import Chem
            from rdkit.Chem import Descriptors, Lipinski, MolSurf, AllChem, Fragments
            from rdkit.Chem.Scaffolds import MurckoScaffold
            from rdkit.Chem import Draw
            from rdkit.Chem.Draw import rdMolDraw2D
            
            # Try to import IPythonConsole, but don't fail if it's not available
            try:
                from rdkit.Chem.Draw import IPythonConsole
                self.add_test_result("RDKIT-1.1", "IPython Import", "Passed", "IPython is available for RDKit visualization")
            except ImportError:
                self.add_test_result("RDKIT-1.1", "IPython Import", "Skipped", "IPython is not available. Some visualization features may be limited.")
            
            self.add_test_result("RDKIT-1", "RDKit Import", "Passed", "RDKit import successful")
        except ImportError as e:
            self.add_test_result("RDKIT-1", "RDKit Import", "Failed", f"RDKit import failed: {str(e)}")

    def test_molecule_parsing(self) -> None:
        """Test that molecules can be parsed correctly."""
        try:
            from rdkit import Chem
            
            # Test parsing ethanol
            mol = Chem.MolFromSmiles(self.test_molecules["ethanol"])
            if mol is None:
                self.add_test_result("RDKIT-2.1", "Parse Ethanol", "Failed", "Failed to parse ethanol SMILES")
            else:
                self.add_test_result("RDKIT-2.1", "Parse Ethanol", "Passed", "Successfully parsed ethanol SMILES")
            
            # Test parsing glycerol
            mol = Chem.MolFromSmiles(self.test_molecules["glycerol"])
            if mol is None:
                self.add_test_result("RDKIT-2.2", "Parse Glycerol", "Failed", "Failed to parse glycerol SMILES")
            else:
                self.add_test_result("RDKIT-2.2", "Parse Glycerol", "Passed", "Successfully parsed glycerol SMILES")
            
            # Test parsing DMSO
            mol = Chem.MolFromSmiles(self.test_molecules["dmso"])
            if mol is None:
                self.add_test_result("RDKIT-2.3", "Parse DMSO", "Failed", "Failed to parse DMSO SMILES")
            else:
                self.add_test_result("RDKIT-2.3", "Parse DMSO", "Passed", "Successfully parsed DMSO SMILES")
            
            # Test parsing aspirin
            mol = Chem.MolFromSmiles(self.test_molecules["aspirin"])
            if mol is None:
                self.add_test_result("RDKIT-2.4", "Parse Aspirin", "Failed", "Failed to parse aspirin SMILES")
            else:
                self.add_test_result("RDKIT-2.4", "Parse Aspirin", "Passed", "Successfully parsed aspirin SMILES")
        except Exception as e:
            self.add_test_result("RDKIT-2", "Molecule Parsing", "Failed", f"Error testing molecule parsing: {str(e)}")

    def test_property_calculation(self) -> None:
        """Test that molecular properties can be calculated correctly."""
        try:
            from rdkit import Chem
            from rdkit.Chem import Descriptors, Lipinski
            
            # Test calculating properties for ethanol
            mol = Chem.MolFromSmiles(self.test_molecules["ethanol"])
            if mol is None:
                self.add_test_result("RDKIT-3.1", "Ethanol Properties", "Skipped", "Failed to parse ethanol SMILES")
            else:
                # Calculate properties
                mw = Descriptors.MolWt(mol)
                logp = Descriptors.MolLogP(mol)
                tpsa = Descriptors.TPSA(mol)
                h_donors = Lipinski.NumHDonors(mol)
                h_acceptors = Lipinski.NumHAcceptors(mol)
                
                # Check that properties are reasonable
                if abs(mw - 46.07) < 0.1:
                    self.add_test_result("RDKIT-3.1.1", "Ethanol Molecular Weight", "Passed", f"Calculated MW: {mw:.2f}")
                else:
                    self.add_test_result("RDKIT-3.1.1", "Ethanol Molecular Weight", "Failed", f"Calculated MW: {mw:.2f}, expected ~46.07")
                
                if h_donors == 1:
                    self.add_test_result("RDKIT-3.1.2", "Ethanol H-Donors", "Passed", f"Calculated H-Donors: {h_donors}")
                else:
                    self.add_test_result("RDKIT-3.1.2", "Ethanol H-Donors", "Failed", f"Calculated H-Donors: {h_donors}, expected 1")
                
                if h_acceptors == 1:
                    self.add_test_result("RDKIT-3.1.3", "Ethanol H-Acceptors", "Passed", f"Calculated H-Acceptors: {h_acceptors}")
                else:
                    self.add_test_result("RDKIT-3.1.3", "Ethanol H-Acceptors", "Failed", f"Calculated H-Acceptors: {h_acceptors}, expected 1")
            
            # Test calculating properties for glycerol
            mol = Chem.MolFromSmiles(self.test_molecules["glycerol"])
            if mol is None:
                self.add_test_result("RDKIT-3.2", "Glycerol Properties", "Skipped", "Failed to parse glycerol SMILES")
            else:
                # Calculate properties
                mw = Descriptors.MolWt(mol)
                logp = Descriptors.MolLogP(mol)
                tpsa = Descriptors.TPSA(mol)
                h_donors = Lipinski.NumHDonors(mol)
                h_acceptors = Lipinski.NumHAcceptors(mol)
                
                # Check that properties are reasonable
                if abs(mw - 92.09) < 0.1:
                    self.add_test_result("RDKIT-3.2.1", "Glycerol Molecular Weight", "Passed", f"Calculated MW: {mw:.2f}")
                else:
                    self.add_test_result("RDKIT-3.2.1", "Glycerol Molecular Weight", "Failed", f"Calculated MW: {mw:.2f}, expected ~92.09")
                
                if h_donors == 3:
                    self.add_test_result("RDKIT-3.2.2", "Glycerol H-Donors", "Passed", f"Calculated H-Donors: {h_donors}")
                else:
                    self.add_test_result("RDKIT-3.2.2", "Glycerol H-Donors", "Failed", f"Calculated H-Donors: {h_donors}, expected 3")
                
                if h_acceptors == 3:
                    self.add_test_result("RDKIT-3.2.3", "Glycerol H-Acceptors", "Passed", f"Calculated H-Acceptors: {h_acceptors}")
                else:
                    self.add_test_result("RDKIT-3.2.3", "Glycerol H-Acceptors", "Failed", f"Calculated H-Acceptors: {h_acceptors}, expected 3")
        except Exception as e:
            self.add_test_result("RDKIT-3", "Property Calculation", "Failed", f"Error testing property calculation: {str(e)}")

    def test_visualization(self) -> None:
        """Test that molecules can be visualized correctly."""
        try:
            from rdkit import Chem
            from rdkit.Chem import Draw
            from rdkit.Chem.Draw import rdMolDraw2D
            
            # Test visualizing ethanol
            mol = Chem.MolFromSmiles(self.test_molecules["ethanol"])
            if mol is None:
                self.add_test_result("RDKIT-4.1", "Visualize Ethanol", "Skipped", "Failed to parse ethanol SMILES")
            else:
                # Generate SVG
                drawer = rdMolDraw2D.MolDraw2DSVG(300, 200)
                drawer.DrawMolecule(mol)
                drawer.FinishDrawing()
                svg = drawer.GetDrawingText()
                
                # Check that SVG is valid
                if svg and "<svg" in svg and "</svg>" in svg:
                    self.add_test_result("RDKIT-4.1", "Visualize Ethanol", "Passed", "Generated valid SVG")
                else:
                    self.add_test_result("RDKIT-4.1", "Visualize Ethanol", "Failed", "Failed to generate valid SVG")
            
            # Test visualizing glycerol
            mol = Chem.MolFromSmiles(self.test_molecules["glycerol"])
            if mol is None:
                self.add_test_result("RDKIT-4.2", "Visualize Glycerol", "Skipped", "Failed to parse glycerol SMILES")
            else:
                # Generate SVG
                drawer = rdMolDraw2D.MolDraw2DSVG(300, 200)
                drawer.DrawMolecule(mol)
                drawer.FinishDrawing()
                svg = drawer.GetDrawingText()
                
                # Check that SVG is valid
                if svg and "<svg" in svg and "</svg>" in svg:
                    self.add_test_result("RDKIT-4.2", "Visualize Glycerol", "Passed", "Generated valid SVG")
                else:
                    self.add_test_result("RDKIT-4.2", "Visualize Glycerol", "Failed", "Failed to generate valid SVG")
        except Exception as e:
            self.add_test_result("RDKIT-4", "Visualization", "Failed", f"Error testing visualization: {str(e)}")

    def test_substructure_search(self) -> None:
        """Test that substructure searching works correctly."""
        try:
            from rdkit import Chem
            
            # Test searching for hydroxyl group in ethanol
            mol = Chem.MolFromSmiles(self.test_molecules["ethanol"])
            pattern = Chem.MolFromSmarts("[OH]")
            if mol is None or pattern is None:
                self.add_test_result("RDKIT-5.1", "Hydroxyl in Ethanol", "Skipped", "Failed to parse molecules")
            else:
                matches = mol.GetSubstructMatches(pattern)
                if matches and len(matches) == 1:
                    self.add_test_result("RDKIT-5.1", "Hydroxyl in Ethanol", "Passed", f"Found {len(matches)} hydroxyl groups")
                else:
                    self.add_test_result("RDKIT-5.1", "Hydroxyl in Ethanol", "Failed", f"Found {len(matches)} hydroxyl groups, expected 1")
            
            # Test searching for hydroxyl groups in glycerol
            mol = Chem.MolFromSmiles(self.test_molecules["glycerol"])
            pattern = Chem.MolFromSmarts("[OH]")
            if mol is None or pattern is None:
                self.add_test_result("RDKIT-5.2", "Hydroxyl in Glycerol", "Skipped", "Failed to parse molecules")
            else:
                matches = mol.GetSubstructMatches(pattern)
                if matches and len(matches) == 3:
                    self.add_test_result("RDKIT-5.2", "Hydroxyl in Glycerol", "Passed", f"Found {len(matches)} hydroxyl groups")
                else:
                    self.add_test_result("RDKIT-5.2", "Hydroxyl in Glycerol", "Failed", f"Found {len(matches)} hydroxyl groups, expected 3")
            
            # Test searching for sulfoxide group in DMSO
            mol = Chem.MolFromSmiles(self.test_molecules["dmso"])
            pattern = Chem.MolFromSmarts("[#16]=O")
            if mol is None or pattern is None:
                self.add_test_result("RDKIT-5.3", "Sulfoxide in DMSO", "Skipped", "Failed to parse molecules")
            else:
                matches = mol.GetSubstructMatches(pattern)
                if matches and len(matches) == 1:
                    self.add_test_result("RDKIT-5.3", "Sulfoxide in DMSO", "Passed", f"Found {len(matches)} sulfoxide groups")
                else:
                    self.add_test_result("RDKIT-5.3", "Sulfoxide in DMSO", "Failed", f"Found {len(matches)} sulfoxide groups, expected 1")
        except Exception as e:
            self.add_test_result("RDKIT-5", "Substructure Search", "Failed", f"Error testing substructure search: {str(e)}")

    def test_similarity_calculation(self) -> None:
        """Test that molecular similarity can be calculated correctly."""
        try:
            from rdkit import Chem
            from rdkit.Chem import AllChem
            from rdkit import DataStructs
            
            # Test similarity between ethanol and itself (should be 1.0)
            mol1 = Chem.MolFromSmiles(self.test_molecules["ethanol"])
            mol2 = Chem.MolFromSmiles(self.test_molecules["ethanol"])
            if mol1 is None or mol2 is None:
                self.add_test_result("RDKIT-6.1", "Self Similarity", "Skipped", "Failed to parse molecules")
            else:
                fp1 = AllChem.GetMorganFingerprintAsBitVect(mol1, 2)
                fp2 = AllChem.GetMorganFingerprintAsBitVect(mol2, 2)
                similarity = DataStructs.TanimotoSimilarity(fp1, fp2)
                
                if abs(similarity - 1.0) < 0.001:
                    self.add_test_result("RDKIT-6.1", "Self Similarity", "Passed", f"Calculated similarity: {similarity:.2f}")
                else:
                    self.add_test_result("RDKIT-6.1", "Self Similarity", "Failed", f"Calculated similarity: {similarity:.2f}, expected 1.0")
            
            # Test similarity between ethanol and glycerol
            mol1 = Chem.MolFromSmiles(self.test_molecules["ethanol"])
            mol2 = Chem.MolFromSmiles(self.test_molecules["glycerol"])
            if mol1 is None or mol2 is None:
                self.add_test_result("RDKIT-6.2", "Ethanol-Glycerol Similarity", "Skipped", "Failed to parse molecules")
            else:
                fp1 = AllChem.GetMorganFingerprintAsBitVect(mol1, 2)
                fp2 = AllChem.GetMorganFingerprintAsBitVect(mol2, 2)
                similarity = DataStructs.TanimotoSimilarity(fp1, fp2)
                
                # Similarity should be between 0 and 1
                if 0 <= similarity <= 1:
                    self.add_test_result("RDKIT-6.2", "Ethanol-Glycerol Similarity", "Passed", f"Calculated similarity: {similarity:.2f}")
                else:
                    self.add_test_result("RDKIT-6.2", "Ethanol-Glycerol Similarity", "Failed", f"Calculated similarity: {similarity:.2f}, expected value between 0 and 1")
            
            # Test similarity between ethanol and DMSO
            mol1 = Chem.MolFromSmiles(self.test_molecules["ethanol"])
            mol2 = Chem.MolFromSmiles(self.test_molecules["dmso"])
            if mol1 is None or mol2 is None:
                self.add_test_result("RDKIT-6.3", "Ethanol-DMSO Similarity", "Skipped", "Failed to parse molecules")
            else:
                fp1 = AllChem.GetMorganFingerprintAsBitVect(mol1, 2)
                fp2 = AllChem.GetMorganFingerprintAsBitVect(mol2, 2)
                similarity = DataStructs.TanimotoSimilarity(fp1, fp2)
                
                # Similarity should be between 0 and 1
                if 0 <= similarity <= 1:
                    self.add_test_result("RDKIT-6.3", "Ethanol-DMSO Similarity", "Passed", f"Calculated similarity: {similarity:.2f}")
                else:
                    self.add_test_result("RDKIT-6.3", "Ethanol-DMSO Similarity", "Failed", f"Calculated similarity: {similarity:.2f}, expected value between 0 and 1")
        except Exception as e:
            self.add_test_result("RDKIT-6", "Similarity Calculation", "Failed", f"Error testing similarity calculation: {str(e)}")

    def test_scientific_accuracy(self) -> None:
        """Test the scientific accuracy of RDKit calculations."""
        try:
            from rdkit import Chem
            from rdkit.Chem import Descriptors, Lipinski, MolSurf
            
            # Test LogP values
            # Ethanol should have a negative LogP (hydrophilic)
            mol = Chem.MolFromSmiles(self.test_molecules["ethanol"])
            if mol is None:
                self.add_test_result("RDKIT-7.1", "Ethanol LogP", "Skipped", "Failed to parse ethanol SMILES")
            else:
                logp = Descriptors.MolLogP(mol)
                if logp < 0:
                    self.add_test_result("RDKIT-7.1", "Ethanol LogP", "Passed", f"Calculated LogP: {logp:.2f}, correctly hydrophilic")
                else:
                    self.add_test_result("RDKIT-7.1", "Ethanol LogP", "Failed", f"Calculated LogP: {logp:.2f}, expected negative value")
            
            # Aspirin should have a positive LogP (lipophilic)
            mol = Chem.MolFromSmiles(self.test_molecules["aspirin"])
            if mol is None:
                self.add_test_result("RDKIT-7.2", "Aspirin LogP", "Skipped", "Failed to parse aspirin SMILES")
            else:
                logp = Descriptors.MolLogP(mol)
                if logp > 0:
                    self.add_test_result("RDKIT-7.2", "Aspirin LogP", "Passed", f"Calculated LogP: {logp:.2f}, correctly lipophilic")
                else:
                    self.add_test_result("RDKIT-7.2", "Aspirin LogP", "Failed", f"Calculated LogP: {logp:.2f}, expected positive value")
            
            # Test Lipinski's Rule of 5
            # Ethanol should pass all Lipinski rules
            mol = Chem.MolFromSmiles(self.test_molecules["ethanol"])
            if mol is None:
                self.add_test_result("RDKIT-7.3", "Ethanol Lipinski", "Skipped", "Failed to parse ethanol SMILES")
            else:
                mw = Descriptors.MolWt(mol)
                logp = Descriptors.MolLogP(mol)
                h_donors = Lipinski.NumHDonors(mol)
                h_acceptors = Lipinski.NumHAcceptors(mol)
                
                violations = 0
                if mw > 500: violations += 1
                if logp > 5: violations += 1
                if h_donors > 5: violations += 1
                if h_acceptors > 10: violations += 1
                
                if violations == 0:
                    self.add_test_result("RDKIT-7.3", "Ethanol Lipinski", "Passed", "No Lipinski violations")
                else:
                    self.add_test_result("RDKIT-7.3", "Ethanol Lipinski", "Failed", f"{violations} Lipinski violations, expected 0")
        except Exception as e:
            self.add_test_result("RDKIT-7", "Scientific Accuracy", "Failed", f"Error testing scientific accuracy: {str(e)}")

    def save_results(self, filename: str) -> None:
        """Save the test results to a file."""
        try:
            with open(filename, 'w') as f:
                json.dump(self.test_results, f, indent=2)
            logger.info(f"Test results saved to {filename}")
        except Exception as e:
            logger.error(f"Failed to save test results: {str(e)}")

def main():
    """Main function."""
    logger.info("Starting RDKit integration tests")
    
    # Create and run the tests
    test = RDKitIntegrationTest()
    results = test.run_tests()
    
    # Save the results
    test.save_results("rdkit_integration_test_results.json")
    
    # Print summary
    logger.info(f"Test Status: {results['status']}")
    logger.info(f"Total Tests: {results['total_tests']}")
    logger.info(f"Passed Tests: {results['passed_tests']}")
    logger.info(f"Failed Tests: {results['failed_tests']}")
    logger.info(f"Skipped Tests: {results['skipped_tests']}")
    
    # Exit with appropriate status code
    if results["status"] == "Passed":
        logger.info("All tests passed!")
        sys.exit(0)
    else:
        logger.error("Some tests failed. See log for details.")
        sys.exit(1)

if __name__ == "__main__":
    main()