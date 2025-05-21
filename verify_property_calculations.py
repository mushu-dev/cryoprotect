#!/usr/bin/env python3
"""
Verify the accuracy of property calculations by comparing mock_rdkit_formula results 
with actual RDKit calculations for a sample of molecules.
"""

import os
import sys
import json
import psycopg2
from psycopg2.extras import RealDictCursor
from dotenv import load_dotenv

# Try to import rdkit
try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, Lipinski, MolSurf, rdMolDescriptors
    RDKIT_AVAILABLE = True
    print(f"Using RDKit version {Chem.__version__} for calculations")
except ImportError:
    RDKIT_AVAILABLE = False
    print("RDKit not available, attempting installation...")
    try:
        import subprocess
        subprocess.check_call([sys.executable, "-m", "pip", "install", "rdkit", "-q"])
        from rdkit import Chem
        from rdkit.Chem import Descriptors, Lipinski, MolSurf, rdMolDescriptors
        RDKIT_AVAILABLE = True
        print(f"RDKit installed and available: {Chem.__version__}")
    except Exception as e:
        print(f"RDKit installation failed: {e}")
        print("Continuing with mock_rdkit_formula only")
        RDKIT_AVAILABLE = False

# Import our mock RDKit formula module for comparison
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
import mock_rdkit_formula

# Setup logging
import logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Load environment variables
load_dotenv()

def connect_to_db():
    """Connect to the database."""
    db_params = {
        'host': os.getenv('SUPABASE_DB_HOST'),
        'port': os.getenv('SUPABASE_DB_PORT', '5432'),
        'dbname': os.getenv('SUPABASE_DB_NAME', 'postgres'),
        'user': os.getenv('SUPABASE_DB_USER'),
        'password': os.getenv('SUPABASE_DB_PASSWORD'),
        'sslmode': os.getenv('SUPABASE_DB_SSLMODE', 'require')
    }
    
    logger.info(f"Connecting to database")
    return psycopg2.connect(**db_params)

def get_sample_molecules(conn, sample_size=20):
    """Get a sample of molecules with SMILES for verification."""
    query = """
        SELECT id, name, smiles
        FROM molecules
        WHERE smiles IS NOT NULL
        ORDER BY RANDOM()
        LIMIT %s
    """
    
    with conn.cursor(cursor_factory=RealDictCursor) as cursor:
        cursor.execute(query, (sample_size,))
        molecules = cursor.fetchall()
    
    logger.info(f"Got {len(molecules)} sample molecules for verification")
    return molecules

def calculate_rdkit_properties(smiles):
    """Calculate molecular properties using actual RDKit."""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            properties = {
                'Molecular Weight': round(Descriptors.MolWt(mol), 2),
                'LogP': round(Descriptors.MolLogP(mol), 2),
                'TPSA': round(MolSurf.TPSA(mol), 2),
                'Hydrogen Bond Donor Count': Lipinski.NumHDonors(mol),
                'Hydrogen Bond Acceptor Count': Lipinski.NumHAcceptors(mol),
                'Rotatable Bond Count': Descriptors.NumRotatableBonds(mol),
                'Ring Count': Chem.rdMolDescriptors.CalcNumRings(mol),
                'Aromatic Ring Count': Chem.rdMolDescriptors.CalcNumAromaticRings(mol),
                'Heavy Atom Count': mol.GetNumHeavyAtoms(),
                'Molecular Formula': Chem.rdMolDescriptors.CalcMolFormula(mol)
            }
            return properties
        else:
            logger.warning(f"Failed to parse SMILES: {smiles}")
            return None
    except Exception as e:
        logger.error(f"Error calculating RDKit properties: {e}")
        return None

def calculate_mock_properties(smiles):
    """Calculate molecular properties using our mock_rdkit_formula."""
    try:
        properties = {
            'Molecular Weight': mock_rdkit_formula.calculate_molecular_weight(smiles),
            'LogP': mock_rdkit_formula.calculate_logp(smiles),
            'TPSA': mock_rdkit_formula.calculate_tpsa(smiles),
            'Hydrogen Bond Donor Count': mock_rdkit_formula.calculate_h_donors(smiles),
            'Hydrogen Bond Acceptor Count': mock_rdkit_formula.calculate_h_acceptors(smiles),
            'Rotatable Bond Count': mock_rdkit_formula.calculate_rotatable_bonds(smiles),
            'Ring Count': mock_rdkit_formula.calculate_ring_count(smiles),
            'Aromatic Ring Count': mock_rdkit_formula.calculate_aromatic_ring_count(smiles),
            'Heavy Atom Count': mock_rdkit_formula.calculate_heavy_atom_count(smiles),
            'Molecular Formula': mock_rdkit_formula.calculate_molecular_formula(smiles)
        }
        return properties
    except Exception as e:
        logger.error(f"Error calculating mock properties: {e}")
        return None

def evaluate_property_accuracy(rdkit_props, mock_props):
    """Evaluate the accuracy of mock properties compared to RDKit properties."""
    if not rdkit_props or not mock_props:
        return None
    
    accuracy = {}
    
    for prop in rdkit_props:
        if prop in mock_props and rdkit_props[prop] is not None and mock_props[prop] is not None:
            # For numeric properties
            if isinstance(rdkit_props[prop], (int, float)) and isinstance(mock_props[prop], (int, float)):
                # Calculate percent difference
                if rdkit_props[prop] == 0:
                    # Avoid division by zero
                    if mock_props[prop] == 0:
                        accuracy[prop] = 100.0  # Both zero - perfect match
                    else:
                        accuracy[prop] = 0.0  # Mock value not zero - complete mismatch
                else:
                    # Calculate percent error
                    percent_error = abs((mock_props[prop] - rdkit_props[prop]) / rdkit_props[prop]) * 100
                    accuracy[prop] = max(0, 100 - percent_error)
            
            # For string properties like molecular formula
            elif isinstance(rdkit_props[prop], str) and isinstance(mock_props[prop], str):
                if rdkit_props[prop] == mock_props[prop]:
                    accuracy[prop] = 100.0
                else:
                    # Simple string similarity (exact match only for now)
                    accuracy[prop] = 0.0
            
            # For other types, just check equality
            else:
                accuracy[prop] = 100.0 if rdkit_props[prop] == mock_props[prop] else 0.0
    
    # Calculate average accuracy across all properties
    if accuracy:
        overall_accuracy = sum(accuracy.values()) / len(accuracy)
    else:
        overall_accuracy = 0.0
    
    return {
        'property_accuracy': accuracy,
        'overall_accuracy': overall_accuracy
    }

def main():
    """Main function to verify property calculations."""
    import argparse
    parser = argparse.ArgumentParser(description="Verify property calculation accuracy")
    parser.add_argument("--sample-size", type=int, default=20, help="Number of molecules to sample")
    parser.add_argument("--output", type=str, default="property_accuracy_report.json", help="Output file for detailed report")
    args = parser.parse_args()
    
    if not RDKIT_AVAILABLE:
        logger.error("RDKit is required for verification. Run this script in the RDKit container.")
        sys.exit(1)
    
    # Connect to the database
    conn = connect_to_db()
    
    try:
        # Get sample molecules
        molecules = get_sample_molecules(conn, args.sample_size)
        
        if not molecules:
            logger.error("No molecules found for verification. Exiting.")
            return
        
        # Analyze accuracy for each molecule
        all_results = []
        total_accuracy = 0.0
        
        print("\nProperty Calculation Verification Report")
        print("========================================")
        
        for i, molecule in enumerate(molecules, 1):
            molecule_id = molecule['id']
            smiles = molecule['smiles']
            name = molecule['name']
            
            # Calculate properties with both methods
            rdkit_props = calculate_rdkit_properties(smiles)
            mock_props = calculate_mock_properties(smiles)
            
            if rdkit_props and mock_props:
                # Evaluate accuracy
                accuracy_results = evaluate_property_accuracy(rdkit_props, mock_props)
                
                if accuracy_results:
                    molecule_result = {
                        'molecule_id': molecule_id,
                        'name': name,
                        'smiles': smiles,
                        'rdkit_properties': rdkit_props,
                        'mock_properties': mock_props,
                        'accuracy': accuracy_results
                    }
                    
                    all_results.append(molecule_result)
                    total_accuracy += accuracy_results['overall_accuracy']
                    
                    # Print summary for this molecule
                    print(f"\nMolecule {i}: {name}")
                    print(f"SMILES: {smiles}")
                    print(f"Overall Accuracy: {accuracy_results['overall_accuracy']:.2f}%")
                    
                    # Print property comparison for selected important properties
                    key_properties = ['Molecular Weight', 'LogP', 'TPSA', 'Molecular Formula']
                    print("\nProperty Comparison (RDKit vs Mock):")
                    for prop in key_properties:
                        if prop in rdkit_props and prop in mock_props:
                            print(f"  {prop}: {rdkit_props[prop]} vs {mock_props[prop]} " +
                                  f"(Accuracy: {accuracy_results['property_accuracy'].get(prop, 'N/A'):.2f}%)")
        
        # Calculate final overall accuracy
        if all_results:
            average_accuracy = total_accuracy / len(all_results)
            print(f"\nOverall Average Accuracy: {average_accuracy:.2f}%")
            
            # Save detailed report to file
            with open(args.output, 'w') as f:
                json.dump({
                    'summary': {
                        'sample_size': len(all_results),
                        'overall_accuracy': average_accuracy
                    },
                    'detailed_results': all_results
                }, f, indent=2)
            
            print(f"\nDetailed report saved to {args.output}")
        else:
            print("\nNo valid results for comparison.")
    
    except Exception as e:
        logger.error(f"Error verifying property calculations: {e}", exc_info=True)
        sys.exit(1)
    
    finally:
        conn.close()

if __name__ == "__main__":
    main()