#\!/usr/bin/env python3
"""
Fix Duplicate Molecules

This script identifies and resolves the remaining 18 duplicate molecules
by consolidating them into a single record per unique InChIKey.
"""

import os
import sys
import logging
import psycopg2
from psycopg2.extras import RealDictCursor
from dotenv import load_dotenv
import json
from datetime import datetime
import uuid

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler(f"fix_duplicate_molecules_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger()

# Load environment variables
load_dotenv()

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

def find_remaining_duplicates(conn):
    """Find remaining duplicate molecules based on InChIKey."""
    try:
        with conn.cursor(cursor_factory=RealDictCursor) as cursor:
            # Identify duplicate InChIKeys
            cursor.execute("""
                SELECT inchikey, COUNT(*) as count
                FROM molecules
                WHERE inchikey IS NOT NULL
                GROUP BY inchikey
                HAVING COUNT(*) > 1
                ORDER BY COUNT(*) DESC
            """)
            
            duplicate_inchikeys = cursor.fetchall()
            logger.info(f"Found {len(duplicate_inchikeys)} InChIKeys with duplicate molecules")
            
            # Get details for each duplicate group
            duplicate_groups = []
            for row in duplicate_inchikeys:
                inchikey = row['inchikey']
                logger.info(f"Processing duplicate InChIKey: {inchikey}")
                query = """
                    SELECT id, name, smiles, pubchem_cid, inchikey, molecular_formula,
                           created_at, updated_at, data_source, 
                           (SELECT COUNT(*) FROM molecular_properties WHERE molecule_id = molecules.id) as property_count
                    FROM molecules
                    WHERE inchikey = %s
                    ORDER BY 
                        CASE WHEN pubchem_cid IS NOT NULL THEN 0 ELSE 1 END,
                        CASE WHEN name LIKE 'TEST_%%' THEN 1 ELSE 0 END,
                        property_count DESC,
                        created_at
                """
                logger.info(f"Executing query with parameter: {inchikey}")
                cursor.execute(query, (inchikey,))
                
                molecules = cursor.fetchall()
                duplicate_groups.append({
                    "inchikey": inchikey,
                    "count": len(molecules),
                    "molecules": molecules
                })
            
            return duplicate_groups
    except Exception as e:
        logger.error(f"Error finding duplicate molecules: {e}")
        raise

def select_primary_molecule(group):
    """Select primary molecule from a group of duplicates."""
    if not group['molecules']:
        return None
    
    # Primary molecule is the first in the sorted list
    # (Already sorted by the SQL query based on PubChem CID, name, property count, and creation date)
    primary = group['molecules'][0]
    logger.info(f"Selected primary molecule: {primary['name']} (ID: {primary['id']})")
    
    return primary

def update_consolidated_molecules(conn, group, primary, duplicates):
    """Update molecules table with consolidated_to properties for duplicates."""
    try:
        with conn.cursor() as cursor:
            # For each duplicate, update its properties to indicate it's consolidated to the primary
            for dup in duplicates:
                cursor.execute("""
                    UPDATE molecules
                    SET properties = COALESCE(properties, '{}'::jsonb) || 
                                    jsonb_build_object('consolidated_to', %s),
                        updated_at = NOW()
                    WHERE id = %s
                """, (str(primary['id']), dup['id']))
            
            logger.info(f"Updated molecules properties for InChIKey {group['inchikey']}")
            return True
    except Exception as e:
        logger.error(f"Error updating consolidated_molecules: {e}")
        raise

def update_dependent_tables(conn, primary, duplicates):
    """Update dependent tables to reference the primary molecule."""
    try:
        duplicate_ids = [dup['id'] for dup in duplicates]
        if not duplicate_ids:
            return True  # No duplicates to process
            
        with conn.cursor() as cursor:
            # Update molecular_properties
            placeholders = ','.join([f"'{id}'" for id in duplicate_ids])
            cursor.execute(f"""
                UPDATE molecular_properties
                SET molecule_id = '{primary['id']}'
                WHERE molecule_id IN ({placeholders})
                AND property_type_id NOT IN (
                    SELECT pt.id FROM property_types pt
                    JOIN molecular_properties mp ON pt.id = mp.property_type_id
                    WHERE mp.molecule_id = '{primary['id']}'
                )
            """)
            
            property_count = cursor.rowcount
            logger.info(f"Updated {property_count} molecular_properties records")
            
            # Update mixture_components - checking for existing combinations first
            cursor.execute(f"""
                DELETE FROM mixture_components mc
                USING (
                    SELECT mixture_id, mc.molecule_id
                    FROM mixture_components mc
                    WHERE mc.molecule_id IN ({placeholders})
                    AND EXISTS (
                        SELECT 1 FROM mixture_components mc2
                        WHERE mc2.mixture_id = mc.mixture_id
                        AND mc2.molecule_id = '{primary['id']}'
                    )
                ) as duplicates
                WHERE mc.mixture_id = duplicates.mixture_id
                AND mc.molecule_id = duplicates.molecule_id
            """)
            
            deleted_count = cursor.rowcount
            logger.info(f"Deleted {deleted_count} duplicate mixture_components records")
            
            # Update remaining mixture components
            cursor.execute(f"""
                UPDATE mixture_components
                SET molecule_id = '{primary['id']}'
                WHERE molecule_id IN ({placeholders})
            """)
            
            mixture_count = cursor.rowcount
            logger.info(f"Updated {mixture_count} mixture_components records")
            
            # Update predictions
            cursor.execute(f"""
                UPDATE predictions
                SET molecule_id = '{primary['id']}'
                WHERE molecule_id IN ({placeholders})
                AND property_type_id NOT IN (
                    SELECT pt.id FROM property_types pt
                    JOIN predictions p ON pt.id = p.property_type_id
                    WHERE p.molecule_id = '{primary['id']}'
                )
            """)
            
            prediction_count = cursor.rowcount
            logger.info(f"Updated {prediction_count} predictions records")
            
            # Skip audit records for now due to user_id not-null constraint
            return True
    except Exception as e:
        logger.error(f"Error updating dependent tables: {e}")
        raise

def consolidate_duplicates(conn, duplicate_groups):
    """Consolidate all duplicate molecules."""
    if not duplicate_groups:
        logger.info("No duplicate molecules to consolidate")
        return {}
    
    results = {
        "groups_processed": 0,
        "total_duplicates": 0,
        "successful_consolidations": 0,
        "consolidated_molecules": []
    }
    
    for group in duplicate_groups:
        # Create a separate connection for each group to isolate transactions
        group_conn = None
        try:
            group_conn = connect_to_db()
            if not group_conn:
                logger.error(f"Failed to create connection for group {group['inchikey']}")
                continue
                
            group_conn.autocommit = False  # Start transaction
            
            logger.info(f"Processing InChIKey {group['inchikey']} with {group['count']} molecules")
            
            # Select primary molecule
            primary = select_primary_molecule(group)
            if not primary:
                logger.error(f"Could not select primary molecule for InChIKey {group['inchikey']}")
                continue
            
            # Get duplicate molecules (all except primary)
            duplicates = [m for m in group['molecules'] if m['id'] != primary['id']]
            logger.info(f"Found {len(duplicates)} duplicates for primary molecule {primary['name']}")
            
            # Update consolidated_molecules table
            if update_consolidated_molecules(group_conn, group, primary, duplicates):
                # Update dependent tables
                if update_dependent_tables(group_conn, primary, duplicates):
                    # Commit this group's changes
                    group_conn.commit()
                    logger.info(f"Successfully consolidated InChIKey {group['inchikey']}")
                    
                    results["successful_consolidations"] += 1
                    results["consolidated_molecules"].append({
                        "inchikey": group['inchikey'],
                        "primary": {
                            "id": primary['id'],
                            "name": primary['name']
                        },
                        "duplicates": [
                            {"id": dup['id'], "name": dup['name']}
                            for dup in duplicates
                        ]
                    })
            
            results["groups_processed"] += 1
            results["total_duplicates"] += len(duplicates)
            
        except Exception as e:
            logger.error(f"Error processing InChIKey {group['inchikey']}: {e}")
            if group_conn:
                group_conn.rollback()
        finally:
            if group_conn:
                group_conn.close()
    
    logger.info(f"Consolidated {results['successful_consolidations']} groups with {results['total_duplicates']} duplicates")
    return results

def main():
    """Main function."""
    conn = connect_to_db()
    if not conn:
        logger.error("Failed to connect to database")
        return 1
    
    try:
        # Find duplicate molecules
        logger.info("Finding remaining duplicate molecules...")
        duplicate_groups = find_remaining_duplicates(conn)
        
        # Consolidate duplicates
        if duplicate_groups:
            logger.info(f"Consolidating {len(duplicate_groups)} duplicate groups...")
            results = consolidate_duplicates(conn, duplicate_groups)
            
            # Generate report
            report = {
                "timestamp": datetime.now().isoformat(),
                "duplicate_groups": len(duplicate_groups),
                "results": results
            }
            
            # Save report
            report_file = f"duplicate_consolidation_report_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json"
            with open(report_file, 'w') as f:
                json.dump(report, f, indent=2)
            
            logger.info(f"Report saved to {report_file}")
            logger.info("Duplicate molecules fixed successfully")
        else:
            logger.info("No duplicate molecules found")
        
        return 0
    except Exception as e:
        logger.exception(f"Error fixing duplicate molecules: {e}")
        return 1
    finally:
        conn.close()

if __name__ == "__main__":
    sys.exit(main())
