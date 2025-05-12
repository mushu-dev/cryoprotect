#!/usr/bin/env python3
"""
Standardize property types in the CryoProtect database.

This script:
1. Identifies duplicate property types with different capitalizations
2. Creates a mapping between variant property names
3. Standardizes property type names
4. Updates all molecular_properties records to use the standardized types
5. Removes deprecated property types

Usage:
    python standardize_property_types.py [--dry-run]
"""

import os
import sys
import json
import logging
import argparse
from datetime import datetime
from typing import Dict, List, Set, Tuple, Any, Optional
import psycopg2
from psycopg2.extras import RealDictCursor
from dotenv import load_dotenv

# Load environment variables
load_dotenv()

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler(f'property_types_standardization_{datetime.now().strftime("%Y%m%d_%H%M%S")}.log'),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

# Property type standardization mapping
# Format: "current_name": "standardized_name"
PROPERTY_TYPE_MAPPING = {
    # Case standardization
    "logp": "LogP",
    "logP": "LogP",
    "xlogp3": "XLogP3",
    "XLogP3": "XLogP3",
    "molecular_weight": "Molecular Weight",
    "mw_freebase": "Molecular Weight",
    "mw_monoisotopic": "Monoisotopic Mass",
    "full_mwt": "Molecular Weight",
    
    # Property name standardization
    "h_bond_donors": "Hydrogen Bond Donor Count",
    "hbd": "Hydrogen Bond Donor Count",
    "hbd_lipinski": "Hydrogen Bond Donor Count",
    
    "h_bond_acceptors": "Hydrogen Bond Acceptor Count",
    "hba": "Hydrogen Bond Acceptor Count",
    "hba_lipinski": "Hydrogen Bond Acceptor Count",
    
    "heavy_atoms": "Heavy Atom Count",
    
    "rotatable_bonds": "Rotatable Bond Count",
    "rtb": "Rotatable Bond Count",
    
    "psa": "TPSA",
    "tpsa": "TPSA",
    "topological polar surface area": "TPSA",
    "Topological Polar Surface Area": "TPSA",
    
    "aromatic_rings": "Aromatic Ring Count",
}

# Properties to mark as deprecated (will be kept but marked as not in use)
DEPRECATED_PROPERTIES = [
    "test_unknown_prop_2",
    "test_property_TEST_PUBLIC_203357",
    "Unknown Property Type"
]

def connect_to_db() -> psycopg2.extensions.connection:
    """Connect to the database."""
    db_params = {
        'host': os.getenv('SUPABASE_DB_HOST'),
        'port': os.getenv('SUPABASE_DB_PORT', '5432'),
        'dbname': os.getenv('SUPABASE_DB_NAME', 'postgres'),
        'user': os.getenv('SUPABASE_DB_USER'),
        'password': os.getenv('SUPABASE_DB_PASSWORD'),
        'sslmode': os.getenv('SUPABASE_DB_SSLMODE', 'require')
    }
    
    # Log connection attempt without password
    log_params = db_params.copy()
    if 'password' in log_params:
        log_params['password'] = '********'
    logger.info(f"Connecting to database with parameters: {log_params}")
    
    return psycopg2.connect(**db_params)

def get_property_types(conn: psycopg2.extensions.connection) -> List[Dict[str, Any]]:
    """Get all property types from the database."""
    with conn.cursor(cursor_factory=RealDictCursor) as cursor:
        cursor.execute("""
            SELECT 
                pt.id, 
                pt.name, 
                pt.description, 
                pt.data_type, 
                pt.units,
                COUNT(mp.id) as usage_count
            FROM 
                property_types pt
            LEFT JOIN 
                molecular_properties mp ON pt.id = mp.property_type_id
            GROUP BY 
                pt.id, pt.name, pt.description, pt.data_type, pt.units
            ORDER BY 
                LOWER(pt.name)
        """)
        return cursor.fetchall()

def find_duplicate_property_types(property_types: List[Dict[str, Any]]) -> Dict[str, List[Dict[str, Any]]]:
    """
    Find duplicate property types based on case-insensitive name comparison.
    
    Returns a dictionary where keys are lowercase names and values are lists of
    property types with that name (case-insensitive).
    """
    duplicates = {}
    
    # Group by lowercase name
    for prop_type in property_types:
        name_lower = prop_type['name'].lower()
        if name_lower not in duplicates:
            duplicates[name_lower] = []
        duplicates[name_lower].append(prop_type)
    
    # Filter to only include names with multiple entries
    return {name: props for name, props in duplicates.items() if len(props) > 1}

def create_standardization_plan(
    property_types: List[Dict[str, Any]],
    mapping: Dict[str, str]
) -> Tuple[Dict[str, Dict[str, Any]], Dict[str, str]]:
    """
    Create a plan for standardizing property types.
    
    Returns:
        - Dictionary of standardized property types (keyed by standardized name)
        - Mapping from old property type IDs to standardized property type IDs
    """
    standardized_types = {}
    id_mapping = {}
    
    # First pass: identify the best property type for each standardized name
    for prop_type in property_types:
        current_name = prop_type['name']
        std_name = mapping.get(current_name.lower(), current_name)
        
        # Skip deprecated properties
        if current_name in DEPRECATED_PROPERTIES:
            logger.info(f"Skipping deprecated property: {current_name}")
            continue
        
        # If we haven't seen this standardized name yet, or this one has more usage,
        # consider it as the canonical version
        if (std_name not in standardized_types or 
            prop_type['usage_count'] > standardized_types[std_name]['usage_count']):
            standardized_types[std_name] = prop_type
    
    # Second pass: create mapping from old IDs to standardized IDs
    for prop_type in property_types:
        current_name = prop_type['name']
        
        # Skip deprecated properties
        if current_name in DEPRECATED_PROPERTIES:
            continue
            
        std_name = mapping.get(current_name.lower(), current_name)
        
        # Map this property type ID to the standardized version's ID
        if std_name in standardized_types and prop_type['id'] != standardized_types[std_name]['id']:
            id_mapping[prop_type['id']] = standardized_types[std_name]['id']
    
    return standardized_types, id_mapping

def update_property_types(
    conn: psycopg2.extensions.connection,
    standardized_types: Dict[str, Dict[str, Any]],
    id_mapping: Dict[str, str],
    dry_run: bool = False
) -> Dict[str, Any]:
    """
    Update property types in the database.
    
    Args:
        conn: Database connection
        standardized_types: Dictionary of standardized property types
        id_mapping: Mapping from old property type IDs to standardized property type IDs
        dry_run: If True, don't actually update the database
        
    Returns:
        Dictionary with statistics about the update
    """
    stats = {
        "properties_updated": 0,
        "names_standardized": 0,
        "records_remapped": 0,
        "deprecated_properties": len(DEPRECATED_PROPERTIES),
    }
    
    with conn.cursor() as cursor:
        # 1. Update property type names to standardized versions
        for std_name, prop_type in standardized_types.items():
            current_name = prop_type['name']
            if current_name != std_name:
                stats["names_standardized"] += 1
                if not dry_run:
                    cursor.execute("""
                        UPDATE property_types
                        SET name = %s, updated_at = NOW()
                        WHERE id = %s
                    """, (std_name, prop_type['id']))
                logger.info(f"Updated property type name: {current_name} -> {std_name}")
        
        # 2. Update molecular_properties to use standardized property type IDs
        for old_id, new_id in id_mapping.items():
            if not dry_run:
                cursor.execute("""
                    UPDATE molecular_properties
                    SET property_type_id = %s, updated_at = NOW()
                    WHERE property_type_id = %s
                """, (new_id, old_id))
                
                # Get count of updated rows
                count = cursor.rowcount
                stats["records_remapped"] += count
                logger.info(f"Remapped {count} properties from type ID {old_id} to {new_id}")
        
        # 3. Mark deprecated properties
        if not dry_run:
            placeholders = ', '.join(['%s'] * len(DEPRECATED_PROPERTIES))
            if DEPRECATED_PROPERTIES:
                cursor.execute(f"""
                    UPDATE property_types
                    SET description = CONCAT(description, ' [DEPRECATED]')
                    WHERE name IN ({placeholders})
                """, DEPRECATED_PROPERTIES)
                logger.info(f"Marked {len(DEPRECATED_PROPERTIES)} properties as deprecated")
    
    if not dry_run:
        conn.commit()
        logger.info("Changes committed to database")
    else:
        conn.rollback()
        logger.info("Dry run: rolled back all changes")
    
    return stats

def main():
    parser = argparse.ArgumentParser(description='Standardize property types in the database')
    parser.add_argument('--dry-run', action='store_true', help='Perform a dry run without making changes')
    args = parser.parse_args()
    
    dry_run = args.dry_run
    if dry_run:
        logger.info("Running in dry-run mode - no changes will be made")
    
    try:
        # Connect to database
        conn = connect_to_db()
        
        # Get all property types
        property_types = get_property_types(conn)
        logger.info(f"Found {len(property_types)} property types in the database")
        
        # Find duplicates
        duplicates = find_duplicate_property_types(property_types)
        
        # Log duplicate property types
        for name_lower, props in duplicates.items():
            logger.info(f"Found duplicate property type (case-insensitive): {name_lower}")
            for prop in props:
                logger.info(f"  - {prop['name']} (ID: {prop['id']}, Usage: {prop['usage_count']})")
        
        # Create standardization plan
        standardized_types, id_mapping = create_standardization_plan(
            property_types, PROPERTY_TYPE_MAPPING
        )
        
        # Log standardization plan
        logger.info(f"Standardization plan created for {len(standardized_types)} property types")
        logger.info(f"Will remap {len(id_mapping)} property types to standardized versions")
        
        # Update property types
        stats = update_property_types(conn, standardized_types, id_mapping, dry_run)
        
        # Log results
        if dry_run:
            logger.info("DRY RUN SUMMARY:")
        else:
            logger.info("UPDATE SUMMARY:")
            
        logger.info(f"Names standardized: {stats['names_standardized']}")
        logger.info(f"Records remapped: {stats['records_remapped']}")
        logger.info(f"Deprecated properties: {stats['deprecated_properties']}")
        
        # Create report file
        report = {
            "timestamp": datetime.now().isoformat(),
            "dry_run": dry_run,
            "property_counts": {
                "total": len(property_types),
                "duplicates": sum(len(props) for props in duplicates.values()),
                "standardized": len(standardized_types),
                "deprecated": len(DEPRECATED_PROPERTIES)
            },
            "updates": stats,
            "duplicate_groups": {name: [p['name'] for p in props] for name, props in duplicates.items()},
            "id_mapping": {old_id: new_id for old_id, new_id in id_mapping.items()},
            "name_mapping": {old: new for old, new in PROPERTY_TYPE_MAPPING.items()}
        }
        
        report_file = f"property_types_standardization_report_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json"
        with open(report_file, 'w') as f:
            json.dump(report, f, indent=2)
        logger.info(f"Report saved to {report_file}")
        
        # Return success
        return 0
        
    except Exception as e:
        logger.exception(f"Error standardizing property types: {e}")
        return 1
    
if __name__ == "__main__":
    sys.exit(main())