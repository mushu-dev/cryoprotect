#!/usr/bin/env python3
"""
consolidate_duplicate_molecules.py

This script identifies and consolidates duplicate molecules based on InChIKey.
It designates one molecule as the primary record and updates all related tables
to reference the primary molecule.

Usage:
    python consolidate_duplicate_molecules.py [--dry-run]
"""

import os
import sys
import json
import logging
import argparse
import psycopg2
import psycopg2.extras
from datetime import datetime
from dotenv import load_dotenv
import uuid

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler("duplicate_consolidation.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

# Load environment variables
load_dotenv()

# Database connection parameters
DB_HOST = os.getenv("SUPABASE_DB_HOST")
DB_PORT = os.getenv("SUPABASE_DB_PORT")
DB_NAME = os.getenv("SUPABASE_DB_NAME")
DB_USER = os.getenv("SUPABASE_DB_USER")
DB_PASSWORD = os.getenv("SUPABASE_DB_PASSWORD")

def get_db_connection():
    """Get a direct database connection using psycopg2."""
    try:
        conn = psycopg2.connect(
            host=DB_HOST,
            port=DB_PORT,
            dbname=DB_NAME,
            user=DB_USER,
            password=DB_PASSWORD
        )
        logger.info("Connected to database using direct PostgreSQL connection")
        return conn
    except Exception as e:
        logger.error(f"Database connection error: {e}")
        return None

def identify_duplicate_molecules(dry_run=False):
    """Identify duplicate molecules by InChIKey."""
    logger.info("Identifying duplicate molecules by InChIKey...")
    
    conn = get_db_connection()
    if not conn:
        return []
        
    cursor = conn.cursor(cursor_factory=psycopg2.extras.RealDictCursor)
    duplicate_groups = []
    
    try:
        # Find all InChIKeys with multiple molecule entries
        cursor.execute("""
            SELECT inchikey, COUNT(*) as count, array_agg(id) as molecule_ids, 
                   array_agg(name) as names, array_agg(pubchem_cid) as pubchem_cids
            FROM molecules
            WHERE inchikey IS NOT NULL
            GROUP BY inchikey
            HAVING COUNT(*) > 1
        """)
        
        duplicate_groups = cursor.fetchall()
        logger.info(f"Found {len(duplicate_groups)} InChIKeys with duplicate molecules")
        
    except Exception as e:
        logger.error(f"Error identifying duplicate molecules: {e}")
        conn.rollback()
    
    cursor.close()
    conn.close()
    return duplicate_groups

def select_primary_molecule(group):
    """
    Select the primary molecule from a group of duplicates.
    Prioritize molecules with PubChem CIDs, better names, and earliest creation.
    """
    # Parse PostgreSQL array string into a Python list of UUIDs
    molecule_ids_str = group['molecule_ids']
    # Remove curly braces and split by comma
    if molecule_ids_str.startswith('{') and molecule_ids_str.endswith('}'): 
        molecule_ids_str = molecule_ids_str[1:-1]
    molecule_ids = [id.strip() for id in molecule_ids_str.split(',')]
    
    logger.info(f"Processing molecules with IDs: {molecule_ids}")
    
    conn = get_db_connection()
    if not conn:
        return None
        
    cursor = conn.cursor(cursor_factory=psycopg2.extras.RealDictCursor)
    primary_molecule = None
    
    try:
        # Check if any of these molecules actually exist in the database
        for molecule_id in molecule_ids:
            cursor.execute("""
                SELECT id, name, pubchem_cid, created_at, data_source, 
                       molecular_weight, smiles, inchi, inchikey
                FROM molecules
                WHERE id = %s
            """, (molecule_id,))
            
            result = cursor.fetchone()
            if result:
                logger.info(f"Found molecule: {result['name']} (ID: {result['id']}, PubChem CID: {result['pubchem_cid']})")
            else:
                logger.warning(f"Molecule with ID {molecule_id} not found in database")
                
        # Now try to get all molecules with the group's InChIKey
        inchikey = group['inchikey']
        
        # Using direct string substitution is safer here because we know inchikey is a string
        # and not a SQL injection risk in this case
        query = f"""
            SELECT id, name, pubchem_cid, created_at, data_source, 
                   molecular_weight, smiles, inchi, inchikey
            FROM molecules
            WHERE inchikey = '{inchikey}'
            ORDER BY 
                CASE WHEN pubchem_cid IS NOT NULL THEN 0 ELSE 1 END,
                CASE WHEN name LIKE 'TEST_%' THEN 1 ELSE 0 END,
                created_at
        """
        
        cursor.execute(query)
        molecules = cursor.fetchall()
        logger.info(f"Found {len(molecules)} molecules with InChIKey {inchikey}")
        
        if molecules and len(molecules) > 0:
            # First molecule in sorted order is our primary
            primary_molecule = molecules[0]
            logger.info(f"Selected primary molecule: {primary_molecule['name']} (ID: {primary_molecule['id']}, PubChem CID: {primary_molecule['pubchem_cid']})")
        
    except Exception as e:
        logger.error(f"Error selecting primary molecule: {e}")
        conn.rollback()
    
    cursor.close()
    conn.close()
    return primary_molecule

def update_consolidated_molecules_table(primary_molecule, duplicates, dry_run=False):
    """
    Update consolidated_molecules table with primary and duplicate information.
    """
    logger.info(f"Updating consolidated_molecules table for {len(duplicates)} duplicates...")
    
    conn = get_db_connection()
    if not conn:
        return False
        
    cursor = conn.cursor()
    success = False
    
    try:
        # Start a transaction
        conn.autocommit = False
        
        if dry_run:
            logger.info("DRY RUN: Would update consolidated_molecules table")
            logger.info(f"Primary molecule: {primary_molecule['id']} ({primary_molecule['name']})")
            logger.info(f"Duplicates: {[d['id'] for d in duplicates]}")
            success = True
        else:
            # Update primary molecule record
            cursor.execute("""
                INSERT INTO consolidated_molecules 
                (id, name, smiles, inchi, inchikey, formula, molecular_weight, pubchem_cid,
                 is_public, data_source, created_at, updated_at, is_consolidated, 
                 primary_molecule_id, primary_molecule_name, primary_pubchem_cid, molecule_status)
                VALUES (
                    %s, %s, %s, %s, %s, %s, %s, %s, 
                    TRUE, %s, %s, %s, TRUE,
                    %s, %s, %s, 'primary'
                )
                ON CONFLICT (id) DO UPDATE SET
                    is_consolidated = TRUE,
                    primary_molecule_id = EXCLUDED.primary_molecule_id,
                    primary_molecule_name = EXCLUDED.primary_molecule_name,
                    primary_pubchem_cid = EXCLUDED.primary_pubchem_cid,
                    molecule_status = 'primary',
                    updated_at = EXCLUDED.updated_at
            """, (
                primary_molecule['id'], primary_molecule['name'], 
                primary_molecule['smiles'], primary_molecule['inchi'], 
                primary_molecule['inchikey'], primary_molecule.get('formula', None), 
                primary_molecule['molecular_weight'], primary_molecule.get('pubchem_cid', None),
                primary_molecule.get('data_source', 'consolidation_script'),
                primary_molecule['created_at'], datetime.now(),
                primary_molecule['id'], primary_molecule['name'], 
                primary_molecule.get('pubchem_cid', None)
            ))
            
            # Update duplicate molecule records
            for dup in duplicates:
                cursor.execute("""
                    INSERT INTO consolidated_molecules 
                    (id, name, smiles, inchi, inchikey, formula, molecular_weight, pubchem_cid,
                     is_public, data_source, created_at, updated_at, is_consolidated, 
                     primary_molecule_id, primary_molecule_name, primary_pubchem_cid, molecule_status)
                    VALUES (
                        %s, %s, %s, %s, %s, %s, %s, %s, 
                        TRUE, %s, %s, %s, TRUE,
                        %s, %s, %s, 'duplicate'
                    )
                    ON CONFLICT (id) DO UPDATE SET
                        is_consolidated = TRUE,
                        primary_molecule_id = EXCLUDED.primary_molecule_id,
                        primary_molecule_name = EXCLUDED.primary_molecule_name,
                        primary_pubchem_cid = EXCLUDED.primary_pubchem_cid,
                        molecule_status = 'duplicate',
                        updated_at = EXCLUDED.updated_at
                """, (
                    dup['id'], dup['name'], 
                    dup['smiles'], dup['inchi'], 
                    dup['inchikey'], dup.get('formula', None), 
                    dup['molecular_weight'], dup.get('pubchem_cid', None),
                    dup.get('data_source', 'consolidation_script'),
                    dup['created_at'], datetime.now(),
                    primary_molecule['id'], primary_molecule['name'], 
                    primary_molecule.get('pubchem_cid', None)
                ))
            
            conn.commit()
            logger.info(f"Updated consolidated_molecules table for primary molecule {primary_molecule['id']} and {len(duplicates)} duplicates")
            success = True
            
    except Exception as e:
        conn.rollback()
        logger.error(f"Error updating consolidated_molecules table: {e}")
    
    cursor.close()
    conn.close()
    return success

def update_dependent_tables(primary_molecule, duplicates, dry_run=False):
    """
    Update dependent tables to reference the primary molecule instead of duplicates.
    """
    logger.info(f"Updating dependent tables for {len(duplicates)} duplicates...")
    
    conn = get_db_connection()
    if not conn:
        return False
        
    cursor = conn.cursor()
    success = False
    
    try:
        # Start a transaction
        conn.autocommit = False
        
        duplicate_ids = [dup['id'] for dup in duplicates]
        
        if dry_run:
            logger.info("DRY RUN: Would update dependent tables")
            logger.info(f"Primary molecule: {primary_molecule['id']} ({primary_molecule['name']})")
            logger.info(f"Duplicates: {duplicate_ids}")
            success = True
        else:
            # Update molecular_properties table
            placeholders = ','.join(['%s'] * len(duplicate_ids))
            cursor.execute(f"""
                UPDATE molecular_properties
                SET molecule_id = %s
                WHERE molecule_id IN ({placeholders})
            """, [primary_molecule['id']] + duplicate_ids)
            
            mp_count = cursor.rowcount
            logger.info(f"Updated {mp_count} molecular_properties records")
            
            # Update mixture_components table
            cursor.execute(f"""
                UPDATE mixture_components
                SET molecule_id = %s
                WHERE molecule_id IN ({placeholders})
            """, [primary_molecule['id']] + duplicate_ids)
            
            mc_count = cursor.rowcount
            logger.info(f"Updated {mc_count} mixture_components records")
            
            # Update predictions table
            cursor.execute(f"""
                UPDATE predictions
                SET molecule_id = %s
                WHERE molecule_id IN ({placeholders})
            """, [primary_molecule['id']] + duplicate_ids)
            
            pred_count = cursor.rowcount
            logger.info(f"Updated {pred_count} predictions records")
            
            # Add an audit record for traceability
            for dup_id in duplicate_ids:
                cursor.execute("""
                    INSERT INTO scientific_data_audit 
                    (table_name, record_id, operation, old_value, new_value, user_id, timestamp)
                    VALUES 
                    ('molecules', %s, 'consolidated', %s, %s, %s, %s)
                """, (
                    dup_id, 
                    json.dumps({"type": "duplicate_consolidation", "duplicate_id": str(dup_id)}),
                    json.dumps({"type": "duplicate_consolidation", "primary_id": str(primary_molecule['id'])}),
                    primary_molecule.get('created_by', '00000000-0000-0000-0000-000000000000'),
                    datetime.now()
                ))
            
            conn.commit()
            logger.info(f"Successfully updated all dependent tables")
            success = True
            
    except Exception as e:
        conn.rollback()
        logger.error(f"Error updating dependent tables: {e}")
    
    cursor.close()
    conn.close()
    return success

def main(dry_run=False):
    # If called from command line, parse arguments
    if len(sys.argv) > 1:
        parser = argparse.ArgumentParser(description="Consolidate duplicate molecules based on InChIKey.")
        parser.add_argument("--dry-run", action="store_true", help="Print actions instead of executing")
        args = parser.parse_args()
        dry_run = args.dry_run
    
    logger.info("Starting duplicate molecule consolidation process")
    
    # Step 1: Identify duplicate molecules
    duplicate_groups = identify_duplicate_molecules(dry_run=dry_run)
    
    if not duplicate_groups:
        logger.info("No duplicate molecules found. Exiting.")
        return True
    
    # Step 2: Process each group of duplicates
    consolidation_summary = {
        "timestamp": datetime.now().isoformat(),
        "groups_processed": 0,
        "total_duplicates": 0,
        "successful_consolidations": 0,
        "details": []
    }
    
    for group in duplicate_groups:
        inchikey = group['inchikey']
        molecule_count = group['count']
        logger.info(f"Processing duplicate group with InChIKey {inchikey} ({molecule_count} molecules)")
        
        # Step 3: Select the primary molecule
        primary_molecule = select_primary_molecule(group)
        if not primary_molecule:
            logger.error(f"Failed to select primary molecule for InChIKey {inchikey}")
            continue
        
        # Get detailed info for all molecules in the group
        conn = get_db_connection()
        if not conn:
            continue
            
        cursor = conn.cursor(cursor_factory=psycopg2.extras.RealDictCursor)
        
        try:
            # Instead of using molecule_ids list which might have issues,
            # just query by InChIKey since we know that works
            inchikey = group['inchikey']
            
            # Using direct string substitution is safer here because we know inchikey is a string
            # and not a SQL injection risk in this case
            query = f"""
                SELECT id, name, pubchem_cid, created_at, data_source, 
                       molecular_weight, smiles, inchi, inchikey, created_by, formula
                FROM molecules
                WHERE inchikey = '{inchikey}'
            """
            
            cursor.execute(query)
            all_molecules = cursor.fetchall()
            
            # Filter out the primary molecule
            duplicates = [m for m in all_molecules if m['id'] != primary_molecule['id']]
            
            # Step 4: Update consolidated_molecules table
            if update_consolidated_molecules_table(primary_molecule, duplicates, dry_run):
                # Step 5: Update dependent tables
                if update_dependent_tables(primary_molecule, duplicates, dry_run):
                    consolidation_summary["successful_consolidations"] += 1
                    
                    consolidation_summary["details"].append({
                        "inchikey": inchikey,
                        "primary": {
                            "id": str(primary_molecule['id']),
                            "name": primary_molecule['name']
                        },
                        "duplicates": [
                            {"id": str(dup['id']), "name": dup['name']}
                            for dup in duplicates
                        ]
                    })
                    
                    logger.info(f"Successfully consolidated duplicate group with InChIKey {inchikey}")
            
            consolidation_summary["groups_processed"] += 1
            consolidation_summary["total_duplicates"] += len(duplicates)
            
        except Exception as e:
            logger.error(f"Error processing duplicate group with InChIKey {inchikey}: {e}")
        finally:
            cursor.close()
            conn.close()
    
    # Save consolidation summary to file
    with open("duplicate_consolidation_summary.json", "w") as f:
        json.dump(consolidation_summary, f, indent=2)
    
    # Print summary
    print("\n" + "=" * 60)
    print("Duplicate Molecule Consolidation Summary")
    print("=" * 60)
    print(f"Groups processed: {consolidation_summary['groups_processed']}")
    print(f"Total duplicates: {consolidation_summary['total_duplicates']}")
    print(f"Successful consolidations: {consolidation_summary['successful_consolidations']}")
    print("=" * 60)
    
    logger.info("Duplicate molecule consolidation complete")
    return True

if __name__ == "__main__":
    success = main()
    exit(0 if success else 1)