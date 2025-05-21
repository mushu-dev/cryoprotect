#!/usr/bin/env python3
"""
Fix Property Standardization Issues

This script addresses the remaining property standardization issues identified
in the verification report, ensuring consistent naming conventions and
eliminating duplicate property types.
"""

import os
import sys
import psycopg2
from psycopg2.extras import RealDictCursor
from dotenv import load_dotenv
import logging
from datetime import datetime
import json

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler(f"property_standardization_fix_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger()

# Load environment variables
load_dotenv()

# Standardization mapping - use canonical names
PROPERTY_STANDARDIZATION = {
    "logp": "LogP",
    "log_p": "LogP",
    "molecular_weight": "Molecular Weight",
    "mol_weight": "Molecular Weight",
    "mw": "Molecular Weight",
    "tpsa": "TPSA",
    "topological_polar_surface_area": "TPSA",
    "h_donors": "Hydrogen Bond Donor Count",
    "hbd": "Hydrogen Bond Donor Count",
    "donor_count": "Hydrogen Bond Donor Count",
    "h_acceptors": "Hydrogen Bond Acceptor Count",
    "hba": "Hydrogen Bond Acceptor Count",
    "acceptor_count": "Hydrogen Bond Acceptor Count",
    "rotatable_bonds": "Rotatable Bond Count",
    "ring_count": "Ring Count",
    "aromatic_rings": "Aromatic Ring Count",
    "heavy_atoms": "Heavy Atom Count",
    "formula": "Molecular Formula",
    "molecular_formula": "Molecular Formula"
}

def connect_to_db():
    """Connect to the database."""
    db_params = {
        'host': os.getenv('SUPABASE_DB_HOST'),
        'port': os.getenv('SUPABASE_DB_PORT', '5432'),
        'dbname': os.getenv('SUPABASE_DB_NAME', 'postgres'),
        'user': os.getenv('SUPABASE_DB_USER'),
        'password': os.getenv('SUPABASE_DB_PASSWORD')
    }
    
    try:
        conn = psycopg2.connect(**db_params)
        logger.info("Connected to database")
        return conn
    except Exception as e:
        logger.error(f"Database connection error: {e}")
        return None

def find_standardization_candidates(conn):
    """Find property types that need standardization."""
    try:
        standardization_candidates = []
        
        with conn.cursor(cursor_factory=RealDictCursor) as cursor:
            # Get all property types
            cursor.execute("""
                SELECT id, name, description, data_type, units
                FROM property_types
                ORDER BY name
            """)
            
            property_types = cursor.fetchall()
            logger.info(f"Found {len(property_types)} property types")
            
            # Identify candidates for standardization
            for prop in property_types:
                name_lower = prop['name'].lower()
                
                # Check if this property should be standardized
                for pattern, standard_name in PROPERTY_STANDARDIZATION.items():
                    if pattern == name_lower or pattern in name_lower:
                        if prop['name'] != standard_name:
                            standardization_candidates.append({
                                "id": prop['id'],
                                "current_name": prop['name'],
                                "standard_name": standard_name,
                                "description": prop['description'],
                                "data_type": prop['data_type'],
                                "units": prop['units']
                            })
                            break
            
            logger.info(f"Found {len(standardization_candidates)} properties that need standardization")
            return standardization_candidates
    except Exception as e:
        logger.error(f"Error finding standardization candidates: {e}")
        raise

def find_property_duplicates(conn):
    """Find duplicate property types (by standard name)."""
    try:
        duplicates = []
        
        with conn.cursor(cursor_factory=RealDictCursor) as cursor:
            # Find property types with the same standardized name
            std_name_counts = {}
            std_name_props = {}
            
            # Get all property types
            cursor.execute("""
                SELECT id, name, description, data_type, units,
                       (SELECT COUNT(*) FROM molecular_properties WHERE property_type_id = property_types.id) as usage_count
                FROM property_types
                ORDER BY name
            """)
            
            for prop in cursor.fetchall():
                # Get standardized name for this property
                name_lower = prop['name'].lower()
                std_name = None
                
                for pattern, standard_name in PROPERTY_STANDARDIZATION.items():
                    if pattern == name_lower or pattern in name_lower:
                        std_name = standard_name
                        break
                
                if not std_name:
                    std_name = prop['name']  # Use current name if no standard exists
                
                # Track by standardized name
                if std_name not in std_name_counts:
                    std_name_counts[std_name] = 0
                    std_name_props[std_name] = []
                
                std_name_counts[std_name] += 1
                std_name_props[std_name].append({
                    "id": prop['id'],
                    "name": prop['name'],
                    "description": prop['description'],
                    "data_type": prop['data_type'],
                    "units": prop['units'],
                    "usage_count": prop['usage_count']
                })
            
            # Identify duplicates
            for std_name, count in std_name_counts.items():
                if count > 1:
                    duplicates.append({
                        "standard_name": std_name,
                        "count": count,
                        "properties": std_name_props[std_name]
                    })
            
            logger.info(f"Found {len(duplicates)} duplicate property type groups")
            return duplicates
    except Exception as e:
        logger.error(f"Error finding duplicate properties: {e}")
        raise

def select_primary_property(properties):
    """Select the primary property from a group of duplicates."""
    if not properties:
        return None
    
    # Sort by usage count (descending) and then by name
    sorted_props = sorted(properties, key=lambda p: (-p['usage_count'], p['name']))
    
    # Select the most-used property as primary
    primary = sorted_props[0]
    logger.info(f"Selected primary property: {primary['name']} (ID: {primary['id']})")
    
    return primary

def standardize_property_types(conn, standardization_candidates):
    """Apply standardization to property types."""
    if not standardization_candidates:
        logger.info("No properties to standardize")
        return {}
    
    try:
        results = {"updated": 0, "errors": 0, "skipped_duplicates": 0}
        
        # First, check for existing properties with target names
        existing_names = {}
        with conn.cursor(cursor_factory=RealDictCursor) as cursor:
            cursor.execute("SELECT id, name FROM property_types")
            for row in cursor.fetchall():
                existing_names[row['name']] = row['id']
        
        # Group candidates by target name
        name_groups = {}
        for candidate in standardization_candidates:
            if candidate['standard_name'] not in name_groups:
                name_groups[candidate['standard_name']] = []
            name_groups[candidate['standard_name']].append(candidate)
        
        with conn.cursor() as cursor:
            # Process each target name group
            for std_name, candidates in name_groups.items():
                # Check if the standard name already exists
                if std_name in existing_names:
                    # Target property already exists, so these will need merging instead
                    logger.info(f"Standard name '{std_name}' already exists, will merge instead of rename")
                    for candidate in candidates:
                        results["skipped_duplicates"] += 1
                    continue
                
                # If there are multiple candidates for this standard name,
                # rename just the one with the highest usage count
                if len(candidates) > 1:
                    # Get usage count for each candidate
                    for cand in candidates:
                        cursor.execute("""
                            SELECT COUNT(*) FROM molecular_properties 
                            WHERE property_type_id = %s
                        """, (cand['id'],))
                        cand['usage_count'] = cursor.fetchone()[0]
                    
                    # Sort by usage count (descending)
                    candidates.sort(key=lambda c: c.get('usage_count', 0), reverse=True)
                    
                    # Only rename the most used one
                    candidate = candidates[0]
                    cursor.execute("""
                        UPDATE property_types
                        SET name = %s, updated_at = NOW()
                        WHERE id = %s
                    """, (std_name, candidate['id']))
                    
                    results["updated"] += 1
                    logger.info(f"Standardized property: {candidate['current_name']} -> {std_name}")
                    
                    # The rest will be merged later
                    for cand in candidates[1:]:
                        results["skipped_duplicates"] += 1
                else:
                    # Just one candidate, safe to rename
                    candidate = candidates[0]
                    cursor.execute("""
                        UPDATE property_types
                        SET name = %s, updated_at = NOW()
                        WHERE id = %s
                    """, (std_name, candidate['id']))
                    
                    results["updated"] += 1
                    logger.info(f"Standardized property: {candidate['current_name']} -> {std_name}")
        
        logger.info(f"Standardized {results['updated']} properties, skipped {results['skipped_duplicates']} duplicates, with {results['errors']} errors")
        return results
    except Exception as e:
        logger.error(f"Error standardizing properties: {e}")
        conn.rollback()
        raise

def merge_duplicate_properties(conn, duplicate_groups):
    """Merge duplicate property types."""
    if not duplicate_groups:
        logger.info("No duplicate properties to merge")
        return {}
    
    try:
        results = {"groups_merged": 0, "properties_affected": 0, "records_updated": 0, "errors": 0}
        
        for group in duplicate_groups:
            try:
                # Select primary property
                primary = select_primary_property(group['properties'])
                if not primary:
                    logger.warning(f"Could not select primary property for group: {group['standard_name']}")
                    results["errors"] += 1
                    continue
                
                # Get duplicates (all except primary)
                duplicates = [p for p in group['properties'] if p['id'] != primary['id']]
                logger.info(f"Merging {len(duplicates)} duplicates into {primary['name']} (ID: {primary['id']})")
                
                # Update molecular_properties to use the primary property type
                with conn.cursor() as cursor:
                    for dup in duplicates:
                        cursor.execute("""
                            UPDATE molecular_properties
                            SET property_type_id = %s
                            WHERE property_type_id = %s
                        """, (primary['id'], dup['id']))
                        
                        updated = cursor.rowcount
                        results["records_updated"] += updated
                        logger.info(f"Updated {updated} records from {dup['name']} to {primary['name']}")
                        
                        # Delete the duplicate property type
                        cursor.execute("""
                            DELETE FROM property_types
                            WHERE id = %s
                        """, (dup['id'],))
                        
                        logger.info(f"Deleted duplicate property: {dup['name']} (ID: {dup['id']})")
                
                results["properties_affected"] += len(duplicates)
                results["groups_merged"] += 1
                
            except Exception as e:
                logger.error(f"Error merging duplicate group {group['standard_name']}: {e}")
                results["errors"] += 1
        
        logger.info(f"Merged {results['groups_merged']} groups with {results['properties_affected']} duplicate properties")
        logger.info(f"Updated {results['records_updated']} property records")
        logger.info(f"Encountered {results['errors']} errors")
        
        return results
    except Exception as e:
        logger.error(f"Error merging duplicate properties: {e}")
        conn.rollback()
        raise

def main():
    """Main function."""
    conn = connect_to_db()
    if not conn:
        logger.error("Failed to connect to database")
        return 1
    
    try:
        # Start transaction
        conn.autocommit = False
        
        # Find properties needing standardization
        logger.info("Finding properties that need standardization...")
        standardization_candidates = find_standardization_candidates(conn)
        
        # Standardize property types
        logger.info("Standardizing property types...")
        standardization_results = standardize_property_types(conn, standardization_candidates)
        
        # Find duplicate properties after standardization
        logger.info("Finding duplicate properties...")
        duplicate_groups = find_property_duplicates(conn)
        
        # Merge duplicate properties
        logger.info("Merging duplicate properties...")
        merge_results = merge_duplicate_properties(conn, duplicate_groups)
        
        # Commit changes
        conn.commit()
        logger.info("Changes committed to database")
        
        # Generate report
        report = {
            "timestamp": datetime.now().isoformat(),
            "standardization_candidates": len(standardization_candidates),
            "standardization_results": standardization_results,
            "duplicate_groups": len(duplicate_groups),
            "merge_results": merge_results
        }
        
        # Save report
        report_file = f"property_standardization_fix_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json"
        with open(report_file, 'w') as f:
            json.dump(report, f, indent=2)
        
        logger.info(f"Report saved to {report_file}")
        logger.info("Property standardization fixed successfully")
        
        return 0
    except Exception as e:
        logger.exception(f"Error fixing property standardization: {e}")
        conn.rollback()
        return 1
    finally:
        conn.close()

if __name__ == "__main__":
    sys.exit(main())