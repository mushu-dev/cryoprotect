#\!/usr/bin/env python3
"""
Verify ChEMBL Import Data

This script verifies the completeness and correctness of the ChEMBL data
imported into the database. It checks for:
1. All molecules have complete molecular properties
2. All ChEMBL molecules have PubChem CIDs where possible
3. All properties are properly stored in both tables and JSONB
"""

import os
import sys
import json
import logging
import argparse
import time
from datetime import datetime
from typing import Dict, Any, List, Set
import psycopg2
from psycopg2.extras import RealDictCursor

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler("verify_chembl_import.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

# Database connection parameters (from environment or .env file)
DB_PARAMS = {
    'host': os.getenv('SUPABASE_DB_HOST', 'aws-0-us-east-1.pooler.supabase.com'),
    'port': int(os.getenv('SUPABASE_DB_PORT', '5432')),
    'dbname': os.getenv('SUPABASE_DB_NAME', 'postgres'),
    'user': os.getenv('SUPABASE_DB_USER', 'postgres.tsdlmynydfuypiugmkev'),
    'password': os.getenv('SUPABASE_DB_PASSWORD', 'LDHt$rkaM&Gmf3X@LQ37'),
    'sslmode': 'require'
}

class ChEMBLImportVerifier:
    """
    Verifier for ChEMBL import data.
    
    This class checks the completeness and correctness of the ChEMBL data
    imported into the database.
    """
    
    def __init__(self, db_params: Dict[str, Any]):
        """
        Initialize the verifier.
        
        Args:
            db_params: Database connection parameters
        """
        self.db_params = db_params
        self.conn = None
        
        # Connect to database
        self._connect()
    
    def _connect(self):
        """Connect to the database."""
        try:
            self.conn = psycopg2.connect(**self.db_params)
            logger.info("Connected to database")
        except Exception as e:
            logger.error(f"Database connection error: {str(e)}")
            sys.exit(1)
    
    def get_chembl_molecules(self, limit: int = None) -> List[Dict[str, Any]]:
        """
        Get all molecules with ChEMBL IDs.
        
        Args:
            limit: Maximum number of molecules to fetch
            
        Returns:
            List of molecules with ChEMBL IDs
        """
        limit_clause = f"LIMIT {limit}" if limit else ""
        
        query = f"""
        SELECT id, name, chembl_id, pubchem_cid, smiles, inchi, inchikey,
               properties, updated_at
        FROM molecules
        WHERE chembl_id IS NOT NULL
        ORDER BY updated_at DESC
        {limit_clause};
        """
        
        try:
            with self.conn.cursor(cursor_factory=RealDictCursor) as cursor:
                cursor.execute(query)
                return cursor.fetchall()
        except Exception as e:
            logger.error(f"Error fetching ChEMBL molecules: {str(e)}")
            return []
    
    def check_property_completeness(self, molecules: List[Dict[str, Any]]) -> Dict[str, Any]:
        """
        Check if all molecules have complete molecular properties.
        
        Args:
            molecules: List of molecules to check
            
        Returns:
            Dictionary with property completeness statistics
        """
        # Expected properties
        expected_properties = {
            'LogP', 'TPSA', 'Molecular Weight', 'Heavy Atom Count',
            'Hydrogen Bond Donor Count', 'Hydrogen Bond Acceptor Count',
            'Rotatable Bond Count', 'Ring Count', 'Aromatic Ring Count'
        }
        
        results = {
            "total_molecules": len(molecules),
            "molecules_with_properties": 0,
            "molecules_missing_properties": 0,
            "molecules_with_jsonb_properties": 0,
            "molecules_missing_jsonb_properties": 0,
            "property_coverage": {},
            "molecules_by_property_count": {}
        }
        
        # Initialize property coverage
        for prop in expected_properties:
            results["property_coverage"][prop] = 0
        
        # Initialize molecules by property count
        for i in range(len(expected_properties) + 1):
            results["molecules_by_property_count"][i] = 0
        
        # Process each molecule
        for molecule in molecules:
            # Get properties from the JSONB field
            jsonb_properties = molecule.get('properties', {}) or {}
            
            if isinstance(jsonb_properties, str):
                try:
                    jsonb_properties = json.loads(jsonb_properties)
                except json.JSONDecodeError:
                    jsonb_properties = {}
            
            # Check for properties in database tables
            query = """
            SELECT property_type
            FROM molecular_properties
            WHERE molecule_id = %s;
            """
            
            try:
                with self.conn.cursor() as cursor:
                    cursor.execute(query, (molecule['id'],))
                    db_properties = {row[0] for row in cursor.fetchall()}
                    
                    # Count properties found in database tables
                    db_property_count = sum(1 for prop in expected_properties if prop in db_properties)
                    results["molecules_by_property_count"][db_property_count] += 1
                    
                    if db_property_count == len(expected_properties):
                        results["molecules_with_properties"] += 1
                    else:
                        results["molecules_missing_properties"] += 1
                    
                    # Update property coverage
                    for prop in expected_properties:
                        if prop in db_properties:
                            results["property_coverage"][prop] += 1
                    
                    # Check for properties in JSONB field
                    jsonb_property_count = sum(1 for prop in expected_properties if prop in jsonb_properties)
                    
                    if jsonb_property_count == len(expected_properties):
                        results["molecules_with_jsonb_properties"] += 1
                    else:
                        results["molecules_missing_jsonb_properties"] += 1
            except Exception as e:
                logger.error(f"Error checking properties for molecule {molecule['id']}: {str(e)}")
                results["molecules_missing_properties"] += 1
                results["molecules_missing_jsonb_properties"] += 1
        
        # Calculate percentages
        for prop in expected_properties:
            results["property_coverage"][prop] = {
                "count": results["property_coverage"][prop],
                "percentage": (results["property_coverage"][prop] / results["total_molecules"]) * 100
            }
        
        results["property_completeness_percentage"] = (results["molecules_with_properties"] / results["total_molecules"]) * 100
        results["jsonb_property_completeness_percentage"] = (results["molecules_with_jsonb_properties"] / results["total_molecules"]) * 100
        
        return results
    
    def check_pubchem_cross_references(self, molecules: List[Dict[str, Any]]) -> Dict[str, Any]:
        """
        Check if all ChEMBL molecules have PubChem CIDs.
        
        Args:
            molecules: List of molecules to check
            
        Returns:
            Dictionary with cross-reference statistics
        """
        results = {
            "total_molecules": len(molecules),
            "molecules_with_pubchem_cid": 0,
            "molecules_missing_pubchem_cid": 0,
            "molecules_with_inchikey": 0,
            "molecules_without_inchikey": 0
        }
        
        # Process each molecule
        for molecule in molecules:
            # Check for PubChem CID
            if molecule.get('pubchem_cid'):
                results["molecules_with_pubchem_cid"] += 1
            else:
                results["molecules_missing_pubchem_cid"] += 1
            
            # Check for InChIKey (required for cross-reference resolution)
            if molecule.get('inchikey'):
                results["molecules_with_inchikey"] += 1
            else:
                results["molecules_without_inchikey"] += 1
        
        # Calculate percentages
        results["pubchem_cross_reference_percentage"] = (results["molecules_with_pubchem_cid"] / results["total_molecules"]) * 100
        results["inchikey_coverage_percentage"] = (results["molecules_with_inchikey"] / results["total_molecules"]) * 100
        
        return results
    
    def verify_import(self, limit: int = None) -> Dict[str, Any]:
        """
        Verify the ChEMBL import data.
        
        Args:
            limit: Maximum number of molecules to check
            
        Returns:
            Dictionary with verification results
        """
        # Get ChEMBL molecules
        molecules = self.get_chembl_molecules(limit)
        
        if not molecules:
            logger.error("No ChEMBL molecules found in the database")
            return {
                "error": "No ChEMBL molecules found",
                "timestamp": datetime.now().isoformat()
            }
        
        logger.info(f"Found {len(molecules)} ChEMBL molecules in the database")
        
        # Check property completeness
        property_results = self.check_property_completeness(molecules)
        logger.info(f"Property completeness: {property_results['property_completeness_percentage']:.2f}%")
        logger.info(f"JSONB property completeness: {property_results['jsonb_property_completeness_percentage']:.2f}%")
        
        # Check PubChem cross-references
        cross_ref_results = self.check_pubchem_cross_references(molecules)
        logger.info(f"PubChem cross-reference coverage: {cross_ref_results['pubchem_cross_reference_percentage']:.2f}%")
        logger.info(f"InChIKey coverage: {cross_ref_results['inchikey_coverage_percentage']:.2f}%")
        
        # Combine results
        results = {
            "timestamp": datetime.now().isoformat(),
            "molecule_count": len(molecules),
            "property_completeness": property_results,
            "pubchem_cross_references": cross_ref_results
        }
        
        return results

def main():
    """Main function."""
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Verify ChEMBL import data')
    parser.add_argument('--limit', type=int, default=None, help='Maximum number of molecules to check')
    parser.add_argument('--report-file', type=str, default=None, help='Save report to file')
    args = parser.parse_args()
    
    # Create timestamp for report file
    timestamp = time.strftime("%Y%m%d_%H%M%S", time.localtime())
    if args.report_file is None:
        args.report_file = f"reports/chembl_verification_{timestamp}.json"
    
    # Ensure reports directory exists
    os.makedirs(os.path.dirname(args.report_file), exist_ok=True)
    
    # Initialize verifier
    verifier = ChEMBLImportVerifier(DB_PARAMS)
    
    # Run verification
    logger.info("Running verification...")
    results = verifier.verify_import(args.limit)
    
    # Save report to file
    with open(args.report_file, 'w') as f:
        json.dump(results, f, indent=2)
    
    logger.info(f"Saved verification report to {args.report_file}")
    
    # Print summary
    print("\n=== ChEMBL Import Verification Summary ===")
    print(f"Total molecules checked: {results['molecule_count']}")
    
    prop_results = results["property_completeness"]
    print(f"Property completeness: {prop_results['property_completeness_percentage']:.2f}%")
    print(f"JSONB property completeness: {prop_results['jsonb_property_completeness_percentage']:.2f}%")
    
    xref_results = results["pubchem_cross_references"]
    print(f"PubChem cross-reference coverage: {xref_results['pubchem_cross_reference_percentage']:.2f}%")
    print(f"InChIKey coverage: {xref_results['inchikey_coverage_percentage']:.2f}%")
    
    print(f"\nDetailed report saved to: {args.report_file}")
    print("============================================\n")
    
    # Return code based on verification results
    prop_threshold = 95.0  # 95% property completeness
    xref_threshold = 90.0  # 90% cross-reference coverage
    
    if (prop_results['property_completeness_percentage'] >= prop_threshold and
        prop_results['jsonb_property_completeness_percentage'] >= prop_threshold and
        xref_results['pubchem_cross_reference_percentage'] >= xref_threshold):
        print("Verification PASSED: The import meets or exceeds quality thresholds.")
        return 0
    else:
        print("Verification WARNING: The import does not meet all quality thresholds.")
        print(f"Expected property completeness: >= {prop_threshold}%")
        print(f"Expected JSONB property completeness: >= {prop_threshold}%")
        print(f"Expected PubChem cross-reference coverage: >= {xref_threshold}%")
        print("\nConsider running integrated_chembl_import_fix.py to improve data quality.")
        return 1

if __name__ == "__main__":
    sys.exit(main())
