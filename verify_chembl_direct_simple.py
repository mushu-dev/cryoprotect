#!/usr/bin/env python3
"""
Verify ChEMBL import data directly in Supabase database - Simplified Version.

This script connects directly to the Supabase database using hard-coded credentials
from the .env file and verifies the ChEMBL molecule import.
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
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

try:
    import psycopg2
    from psycopg2.extras import RealDictCursor
except ImportError:
    logger.error("psycopg2 is required. Install it with: pip install psycopg2-binary")
    sys.exit(1)

# Direct connection parameters from .env file
DB_PARAMS = {
    'host': 'aws-0-us-east-1.pooler.supabase.com',
    'port': 5432,
    'dbname': 'postgres',
    'user': 'postgres.tsdlmynydfuypiugmkev',
    'password': 'LDHt$rkaM&Gmf3X@LQ37',
    'sslmode': 'require'
}

def connect_to_database() -> Tuple[Optional[Any], str]:
    """Connect to the database and return connection object."""
    try:
        # Try to connect
        conn = psycopg2.connect(**DB_PARAMS)
        return conn, "Connected successfully"
    except Exception as e:
        return None, f"Failed to connect: {str(e)}"

def get_tables(conn) -> List[Dict[str, Any]]:
    """Get list of tables in the public schema."""
    query = """
    SELECT table_name
    FROM information_schema.tables
    WHERE table_schema = 'public'
    ORDER BY table_name;
    """
    
    with conn.cursor(cursor_factory=RealDictCursor) as cursor:
        cursor.execute(query)
        return cursor.fetchall()

def get_molecule_count(conn) -> int:
    """Get total count of molecules in the database."""
    query = """
    SELECT COUNT(*) as count
    FROM molecules;
    """
    
    with conn.cursor(cursor_factory=RealDictCursor) as cursor:
        cursor.execute(query)
        result = cursor.fetchone()
        return result['count'] if result else 0

def get_chembl_molecules(conn, limit: int) -> List[Dict[str, Any]]:
    """Get ChEMBL molecules from the database."""
    query = f"""
    SELECT id, name, smiles, inchi, inchikey, formula, molecular_weight, 
           data_source, pubchem_cid, chembl_id, properties
    FROM molecules
    WHERE chembl_id IS NOT NULL
    ORDER BY created_at DESC
    LIMIT {limit};
    """
    
    with conn.cursor(cursor_factory=RealDictCursor) as cursor:
        cursor.execute(query)
        return cursor.fetchall()

def get_molecular_properties(conn, molecule_ids: List[str]) -> Dict[str, List[Dict[str, Any]]]:
    """Get molecular properties for the given molecule IDs."""
    if not molecule_ids:
        return {}
        
    # Convert list of UUIDs to a string for SQL IN clause
    ids_str = "', '".join(molecule_ids)
    
    query = f"""
    SELECT molecule_id, property_type, property_value, method
    FROM molecular_properties
    WHERE molecule_id IN ('{ids_str}')
    ORDER BY molecule_id, property_type;
    """
    
    with conn.cursor(cursor_factory=RealDictCursor) as cursor:
        cursor.execute(query)
        results = cursor.fetchall()
    
    # Group properties by molecule_id
    grouped = {}
    for prop in results:
        molecule_id = prop['molecule_id']
        if molecule_id not in grouped:
            grouped[molecule_id] = []
        grouped[molecule_id].append(prop)
        
    return grouped

def check_property_completeness(molecules: List[Dict[str, Any]], properties: Dict[str, List[Dict[str, Any]]]) -> Dict[str, Any]:
    """Check if all expected properties are calculated for each molecule."""
    results = {
        "total_molecules": len(molecules),
        "molecules_with_properties": 0,
        "molecules_missing_properties": 0,
        "molecules_with_jsonb_properties": 0,
        "property_counts": {},
        "molecules_with_issues": []
    }
    
    # Define expected properties
    expected_properties = {
        "logp", "tpsa", "molecular_weight", "heavy_atom_count", "h_bond_donor_count", 
        "h_bond_acceptor_count", "rotatable_bond_count", "ring_count"
    }
    
    for molecule in molecules:
        molecule_id = molecule['id']
        has_issues = False
        missing_properties = []
        
        # Check if molecule has properties in the molecular_properties table
        molecule_properties = properties.get(molecule_id, [])
        
        # Check if molecule has properties in the JSONB properties field
        jsonb_properties = molecule.get('properties', {})
        if jsonb_properties:
            results["molecules_with_jsonb_properties"] += 1
        
        # Convert property types to a set for easy checking
        property_types = {prop['property_type'] for prop in molecule_properties}
        
        # Check for missing properties
        for prop in expected_properties:
            if prop not in property_types and (not jsonb_properties or prop not in jsonb_properties):
                missing_properties.append(prop)
                has_issues = True
        
        # Update property counts
        for prop_type in property_types:
            if prop_type not in results["property_counts"]:
                results["property_counts"][prop_type] = 0
            results["property_counts"][prop_type] += 1
        
        # Update molecule counters
        if len(molecule_properties) > 0 or jsonb_properties:
            results["molecules_with_properties"] += 1
        else:
            results["molecules_missing_properties"] += 1
            has_issues = True
        
        # Add molecule to issues list if needed
        if has_issues:
            results["molecules_with_issues"].append({
                "id": molecule_id,
                "chembl_id": molecule.get('chembl_id'),
                "missing_properties": missing_properties,
                "has_molecular_properties": len(molecule_properties) > 0,
                "has_jsonb_properties": bool(jsonb_properties)
            })
    
    return results

def check_cross_references(molecules: List[Dict[str, Any]]) -> Dict[str, Any]:
    """Check if molecules have cross-references to PubChem."""
    results = {
        "total_molecules": len(molecules),
        "molecules_with_pubchem_cid": 0,
        "molecules_with_chembl_id": 0,
        "molecules_with_both": 0,
        "molecules_missing_pubchem": []
    }
    
    for molecule in molecules:
        has_pubchem = molecule.get('pubchem_cid') is not None
        has_chembl = molecule.get('chembl_id') is not None
        
        if has_pubchem:
            results["molecules_with_pubchem_cid"] += 1
        
        if has_chembl:
            results["molecules_with_chembl_id"] += 1
        
        if has_pubchem and has_chembl:
            results["molecules_with_both"] += 1
        
        if has_chembl and not has_pubchem:
            results["molecules_missing_pubchem"].append({
                "id": molecule['id'],
                "chembl_id": molecule['chembl_id']
            })
    
    return results

def generate_report(
    tables: List[Dict[str, Any]], 
    total_molecules: int,
    chembl_molecules: List[Dict[str, Any]], 
    property_results: Dict[str, Any],
    cross_ref_results: Dict[str, Any]
) -> Dict[str, Any]:
    """Generate a comprehensive verification report."""
    report = {
        "timestamp": datetime.now().isoformat(),
        "database_info": {
            "tables": [table['table_name'] for table in tables],
            "total_molecules": total_molecules,
            "chembl_molecules_checked": len(chembl_molecules)
        },
        "property_verification": property_results,
        "cross_reference_verification": cross_ref_results,
        "overall_health": {
            "percent_with_properties": 0,
            "percent_with_cross_refs": 0,
            "status": "unknown"
        }
    }
    
    # Calculate overall health metrics
    if len(chembl_molecules) > 0:
        prop_percent = (property_results["molecules_with_properties"] / len(chembl_molecules)) * 100
        xref_percent = (cross_ref_results["molecules_with_both"] / len(chembl_molecules)) * 100
        
        report["overall_health"]["percent_with_properties"] = prop_percent
        report["overall_health"]["percent_with_cross_refs"] = xref_percent
        
        # Determine overall status
        if prop_percent >= 95 and xref_percent >= 90:
            report["overall_health"]["status"] = "excellent"
        elif prop_percent >= 80 and xref_percent >= 75:
            report["overall_health"]["status"] = "good"
        elif prop_percent >= 60 and xref_percent >= 50:
            report["overall_health"]["status"] = "fair"
        else:
            report["overall_health"]["status"] = "poor"
    
    return report

def main():
    """Main function."""
    limit = 100  # Default limit
    
    logger.info("Verifying ChEMBL import data in Supabase with direct connection")
    
    # Connect to the database
    logger.info("Connecting to the database...")
    conn, message = connect_to_database()
    
    if not conn:
        logger.error(f"Database connection failed: {message}")
        sys.exit(1)
    
    logger.info(message)
    
    try:
        # Get list of tables
        logger.info("Getting list of tables...")
        tables = get_tables(conn)
        logger.info(f"Found {len(tables)} tables in the database")
        
        # Get total molecule count
        logger.info("Getting total molecule count...")
        total_molecules = get_molecule_count(conn)
        logger.info(f"Total molecules in database: {total_molecules}")
        
        # Get ChEMBL molecules
        logger.info(f"Getting up to {limit} ChEMBL molecules...")
        chembl_molecules = get_chembl_molecules(conn, limit)
        logger.info(f"Found {len(chembl_molecules)} ChEMBL molecules")
        
        # Print sample data for verification
        if chembl_molecules:
            logger.info("Sample molecule data:")
            sample = chembl_molecules[0]
            for key, value in sample.items():
                if key not in ('smiles', 'inchi', 'properties'):
                    logger.info(f"  {key}: {value}")
        
        # Get molecular properties
        if chembl_molecules:
            logger.info("Getting molecular properties...")
            molecule_ids = [molecule['id'] for molecule in chembl_molecules]
            molecular_properties = get_molecular_properties(conn, molecule_ids)
            logger.info(f"Found properties for {len(molecular_properties)} molecules")
            
            # Check property completeness
            logger.info("Checking property completeness...")
            property_results = check_property_completeness(chembl_molecules, molecular_properties)
            
            # Check cross-references
            logger.info("Checking cross-references...")
            cross_ref_results = check_cross_references(chembl_molecules)
            
            # Generate report
            logger.info("Generating verification report...")
            report = generate_report(
                tables, 
                total_molecules, 
                chembl_molecules, 
                property_results, 
                cross_ref_results
            )
            
            # Save report to file
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            report_file = f"reports/chembl_data_verification_{timestamp}.json"
            os.makedirs("reports", exist_ok=True)
            
            with open(report_file, 'w') as f:
                json.dump(report, f, indent=2)
            
            logger.info(f"Report saved to {report_file}")
            
            # Print summary
            print("\n===== ChEMBL Import Verification Summary =====")
            print(f"Total molecules in database: {total_molecules}")
            print(f"ChEMBL molecules checked: {len(chembl_molecules)}")
            if len(chembl_molecules) > 0:
                print(f"Molecules with properties: {property_results['molecules_with_properties']} ({property_results['molecules_with_properties']/len(chembl_molecules)*100:.1f}%)")
                print(f"Molecules with both ChEMBL and PubChem IDs: {cross_ref_results['molecules_with_both']} ({cross_ref_results['molecules_with_both']/len(chembl_molecules)*100:.1f}%)")
                print(f"Overall health status: {report['overall_health']['status'].upper()}")
            
                # Print issues summary if any
                if property_results["molecules_with_issues"]:
                    print(f"Found {len(property_results['molecules_with_issues'])} molecules with property issues")
                    
                if cross_ref_results["molecules_missing_pubchem"]:
                    print(f"Found {len(cross_ref_results['molecules_missing_pubchem'])} ChEMBL molecules without PubChem CIDs")
            print("=============================================\n")
                
        else:
            logger.warning("No ChEMBL molecules found in the database!")
    
    finally:
        # Close database connection
        if conn:
            conn.close()
            logger.info("Database connection closed")

if __name__ == "__main__":
    main()