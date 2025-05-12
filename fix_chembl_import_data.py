#!/usr/bin/env python3
"""
Fix ChEMBL import data by updating missing properties and cross-references.

This script integrates both the enhanced property calculator and the ChEMBL to
PubChem cross-reference resolver to fix all issues identified during the
data verification process.
"""

import os
import sys
import logging
import argparse
import time
from typing import Dict, Any, List, Optional
import json

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler("fix_chembl_import.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

# Import our tools
try:
    from enhanced_property_calculator import PropertyCalculator, DB_PARAMS
    from chembl_pubchem_resolver import ChEMBLPubChemResolver
except ImportError as e:
    logger.error(f"Failed to import required modules: {str(e)}")
    logger.error("Make sure enhanced_property_calculator.py and chembl_pubchem_resolver.py are in the current directory.")
    sys.exit(1)

def run_property_updates(args: argparse.Namespace) -> Dict[str, Any]:
    """
    Run property updates using the enhanced property calculator.
    
    Args:
        args: Command line arguments
        
    Returns:
        Statistics about the property update process
    """
    logger.info("=== Starting Property Updates ===")
    
    try:
        # Initialize property calculator
        calculator = PropertyCalculator(DB_PARAMS)
        
        # Process specific molecule if requested
        if args.molecule_id:
            logger.info(f"Processing properties for molecule ID: {args.molecule_id}")
            
            # Get molecule info
            query = "SELECT name, smiles, inchi FROM molecules WHERE id = %s;"
            with calculator.conn.cursor() as cursor:
                cursor.execute(query, (args.molecule_id,))
                molecule = cursor.fetchone()
                
                if not molecule:
                    logger.error(f"Molecule ID {args.molecule_id} not found")
                    return {"error": "Molecule not found"}
                
                name, smiles, inchi = molecule
                
                logger.info(f"Molecule: {name}")
                
                # Calculate properties
                calculated_properties = calculator.calculate_rdkit_properties(smiles, inchi)
                
                if calculated_properties:
                    # Get existing properties
                    existing_properties = calculator.get_existing_properties(args.molecule_id)
                    
                    # Update properties
                    success = calculator.update_molecule_properties(
                        args.molecule_id, 
                        calculated_properties,
                        existing_properties
                    )
                    
                    if success:
                        logger.info(f"Successfully updated properties for molecule {name}")
                        return {"success": True, "molecule": name, "properties_updated": list(calculated_properties.keys())}
                    else:
                        logger.error(f"Failed to update properties for molecule {name}")
                        return {"success": False, "error": "Failed to update properties"}
                else:
                    logger.error(f"Could not calculate properties for molecule {name}")
                    return {"success": False, "error": "Failed to calculate properties"}
        else:
            # Process all molecules with missing properties
            logger.info("Processing all molecules with missing properties")
            
            stats = calculator.update_all_missing_properties(
                batch_size=args.batch_size,
                limit=args.limit
            )
            
            # Print statistics
            logger.info("=== Property Update Statistics ===")
            logger.info(f"Total molecules processed: {stats['total_molecules']}")
            logger.info(f"Successful updates: {stats['successful_updates']}")
            logger.info(f"Failed updates: {stats['failed_updates']}")
            logger.info(f"Skipped (no structure): {stats['skipped_no_structure']}")
            logger.info(f"Skipped (error): {stats['skipped_error']}")
            logger.info("=================================")
            
            return stats
    except Exception as e:
        logger.error(f"Error during property updates: {str(e)}")
        return {"error": str(e)}

def run_cross_reference_updates(args: argparse.Namespace) -> Dict[str, Any]:
    """
    Run cross-reference updates using the ChEMBL to PubChem resolver.
    
    Args:
        args: Command line arguments
        
    Returns:
        Statistics about the cross-reference update process
    """
    logger.info("=== Starting Cross-Reference Updates ===")
    
    try:
        # Initialize resolver
        resolver = ChEMBLPubChemResolver(DB_PARAMS)
        
        if args.chembl_id:
            # Look up a specific ChEMBL ID
            logger.info(f"Looking up ChEMBL ID: {args.chembl_id}")
            
            pubchem_cid = resolver.find_pubchem_by_chembl_id(args.chembl_id)
            
            if pubchem_cid:
                logger.info(f"Found PubChem CID: {pubchem_cid}")
                
                # Get molecule with this ChEMBL ID
                query = "SELECT id, name FROM molecules WHERE chembl_id = %s;"
                with resolver.conn.cursor() as cursor:
                    cursor.execute(query, (args.chembl_id,))
                    molecule = cursor.fetchone()
                    
                    if molecule:
                        molecule_id, name = molecule
                        
                        # Update PubChem CID
                        success = resolver.update_pubchem_cid(molecule_id, pubchem_cid, "chembl_id")
                        
                        if success:
                            logger.info(f"Successfully updated PubChem CID for molecule {name}")
                            return {"success": True, "molecule": name, "pubchem_cid": pubchem_cid}
                        else:
                            logger.error(f"Failed to update PubChem CID for molecule {name}")
                            return {"success": False, "error": "Failed to update PubChem CID"}
                    else:
                        logger.error(f"No molecule found with ChEMBL ID {args.chembl_id}")
                        return {"success": False, "error": "No molecule found with this ChEMBL ID"}
            else:
                logger.error(f"Could not find PubChem CID for ChEMBL ID {args.chembl_id}")
                return {"success": False, "error": "PubChem CID not found"}
        elif args.molecule_id:
            # Process specific molecule
            logger.info(f"Processing cross-references for molecule ID: {args.molecule_id}")
            
            # Get molecule info
            query = "SELECT name, chembl_id, inchikey, smiles FROM molecules WHERE id = %s;"
            with resolver.conn.cursor() as cursor:
                cursor.execute(query, (args.molecule_id,))
                molecule = cursor.fetchone()
                
                if not molecule:
                    logger.error(f"Molecule ID {args.molecule_id} not found")
                    return {"error": "Molecule not found"}
                
                name, chembl_id, inchikey, smiles = molecule
                
                logger.info(f"Molecule: {name}")
                logger.info(f"ChEMBL ID: {chembl_id}")
                
                if not chembl_id:
                    logger.error(f"Molecule {name} does not have a ChEMBL ID")
                    return {"success": False, "error": "No ChEMBL ID"}
                
                # Try to find PubChem CID using different methods
                pubchem_cid = None
                method = ""
                
                # Method 1: Direct ChEMBL ID lookup
                pubchem_cid = resolver.find_pubchem_by_chembl_id(chembl_id)
                if pubchem_cid:
                    method = "chembl_id"
                
                # Method 2: InChIKey lookup
                if not pubchem_cid and inchikey:
                    pubchem_cid = resolver.find_pubchem_by_inchikey(inchikey)
                    if pubchem_cid:
                        method = "inchikey"
                
                # Method 3: SMILES lookup
                if not pubchem_cid and smiles:
                    pubchem_cid = resolver.find_pubchem_by_smiles(smiles)
                    if pubchem_cid:
                        method = "smiles"
                
                # Update PubChem CID if found
                if pubchem_cid:
                    logger.info(f"Found PubChem CID: {pubchem_cid} (method: {method})")
                    
                    success = resolver.update_pubchem_cid(args.molecule_id, pubchem_cid, method)
                    
                    if success:
                        logger.info(f"Successfully updated PubChem CID for molecule {name}")
                        return {"success": True, "molecule": name, "pubchem_cid": pubchem_cid, "method": method}
                    else:
                        logger.error(f"Failed to update PubChem CID for molecule {name}")
                        return {"success": False, "error": "Failed to update PubChem CID"}
                else:
                    logger.error(f"Could not find PubChem CID for molecule {name}")
                    return {"success": False, "error": "PubChem CID not found"}
        else:
            # Process all molecules without PubChem CIDs
            logger.info("Processing all ChEMBL molecules without PubChem CIDs")
            
            stats = resolver.resolve_all_missing_pubchem_cids(
                batch_size=args.batch_size,
                limit=args.limit
            )
            
            # Print statistics
            logger.info("=== PubChem CID Resolution Statistics ===")
            logger.info(f"Total molecules processed: {stats['total_molecules']}")
            logger.info(f"Resolved by ChEMBL ID: {stats['resolved_by_chembl']}")
            logger.info(f"Resolved by InChIKey: {stats['resolved_by_inchikey']}")
            logger.info(f"Resolved by SMILES: {stats['resolved_by_smiles']}")
            logger.info(f"Failed to resolve: {stats['failed_to_resolve']}")
            logger.info("========================================")
            
            return stats
    except Exception as e:
        logger.error(f"Error during cross-reference updates: {str(e)}")
        return {"error": str(e)}

def create_verification_report(property_stats: Dict[str, Any], cross_ref_stats: Dict[str, Any]) -> Dict[str, Any]:
    """
    Create a verification report after running the fixes.
    
    Args:
        property_stats: Statistics from property updates
        cross_ref_stats: Statistics from cross-reference updates
        
    Returns:
        Comprehensive verification report
    """
    # Generate timestamp
    timestamp = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
    
    # Calculate success rates
    property_success_rate = 0
    if property_stats.get('total_molecules', 0) > 0:
        property_success_rate = property_stats.get('successful_updates', 0) / property_stats.get('total_molecules', 1) * 100
    
    cross_ref_success_rate = 0
    if cross_ref_stats.get('total_molecules', 0) > 0:
        cross_ref_success_rate = (
            cross_ref_stats.get('resolved_by_chembl', 0) + 
            cross_ref_stats.get('resolved_by_inchikey', 0) + 
            cross_ref_stats.get('resolved_by_smiles', 0)
        ) / cross_ref_stats.get('total_molecules', 1) * 100
    
    # Generate report
    report = {
        "timestamp": timestamp,
        "property_updates": {
            "total_molecules": property_stats.get('total_molecules', 0),
            "successful_updates": property_stats.get('successful_updates', 0),
            "failed_updates": property_stats.get('failed_updates', 0),
            "skipped_no_structure": property_stats.get('skipped_no_structure', 0),
            "skipped_error": property_stats.get('skipped_error', 0),
            "success_rate": property_success_rate
        },
        "cross_reference_updates": {
            "total_molecules": cross_ref_stats.get('total_molecules', 0),
            "resolved_by_chembl": cross_ref_stats.get('resolved_by_chembl', 0),
            "resolved_by_inchikey": cross_ref_stats.get('resolved_by_inchikey', 0),
            "resolved_by_smiles": cross_ref_stats.get('resolved_by_smiles', 0),
            "failed_to_resolve": cross_ref_stats.get('failed_to_resolve', 0),
            "success_rate": cross_ref_success_rate
        },
        "overall": {
            "total_molecules_processed": property_stats.get('total_molecules', 0) + cross_ref_stats.get('total_molecules', 0),
            "overall_success_rate": (property_success_rate + cross_ref_success_rate) / 2 if property_stats.get('total_molecules', 0) > 0 and cross_ref_stats.get('total_molecules', 0) > 0 else 0
        }
    }
    
    return report

def main():
    """Main function."""
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Fix ChEMBL import data issues')
    parser.add_argument('--batch-size', type=int, default=50, help='Batch size')
    parser.add_argument('--limit', type=int, default=None, help='Maximum number of molecules to process')
    parser.add_argument('--molecule-id', type=str, default=None, help='Process a specific molecule ID')
    parser.add_argument('--chembl-id', type=str, default=None, help='Look up a specific ChEMBL ID')
    parser.add_argument('--properties-only', action='store_true', help='Only update properties')
    parser.add_argument('--cross-refs-only', action='store_true', help='Only update cross-references')
    parser.add_argument('--report-file', type=str, default=None, help='Save report to file')
    args = parser.parse_args()
    
    # Create timestamp for report file
    timestamp = time.strftime("%Y%m%d_%H%M%S", time.localtime())
    if args.report_file is None:
        args.report_file = f"reports/chembl_fixes_{timestamp}.json"
    
    # Ensure reports directory exists
    os.makedirs(os.path.dirname(args.report_file), exist_ok=True)
    
    # Run updates
    property_stats = {"total_molecules": 0}
    cross_ref_stats = {"total_molecules": 0}
    
    # Update properties if requested
    if not args.cross_refs_only:
        property_stats = run_property_updates(args)
    
    # Update cross-references if requested
    if not args.properties_only:
        cross_ref_stats = run_cross_reference_updates(args)
    
    # Create verification report
    report = create_verification_report(property_stats, cross_ref_stats)
    
    # Save report to file
    with open(args.report_file, 'w') as f:
        json.dump(report, f, indent=2)
    
    logger.info(f"Saved report to {args.report_file}")
    
    # Print summary
    print("\n=== ChEMBL Import Fix Summary ===")
    print(f"Property updates: {property_stats.get('successful_updates', 0)} successful, {property_stats.get('failed_updates', 0)} failed")
    print(f"Cross-reference updates: {cross_ref_stats.get('resolved_by_chembl', 0) + cross_ref_stats.get('resolved_by_inchikey', 0) + cross_ref_stats.get('resolved_by_smiles', 0)} successful, {cross_ref_stats.get('failed_to_resolve', 0)} failed")
    print(f"Overall success rate: {report['overall']['overall_success_rate']:.2f}%")
    print(f"Report saved to: {args.report_file}")
    print("================================\n")

if __name__ == "__main__":
    main()