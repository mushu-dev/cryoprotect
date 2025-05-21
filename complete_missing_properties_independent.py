#!/usr/bin/env python3
"""
Complete missing molecular properties in the CryoProtect database.

This script:
1. Identifies molecules with missing essential properties but valid SMILES
2. Uses mock_rdkit_formula to calculate properties from SMILES without requiring RDKit
3. Updates the database with the calculated properties

This script is designed to work without requiring RDKit or containers.
"""

import os
import sys
import time
import json
import logging
import psycopg2
from psycopg2.extras import RealDictCursor
from dotenv import load_dotenv

# Import our custom mock_rdkit_formula module
import mock_rdkit_formula

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler("property_completion.log"),
        logging.StreamHandler()
    ]
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

def get_molecules_needing_properties(conn, limit=None):
    """Get molecules that need property calculations."""
    query = """
        SELECT id, name, smiles
        FROM molecules
        WHERE smiles IS NOT NULL
        AND data_source != 'MCP_Verification_Atomic'
        AND id NOT IN (
            SELECT DISTINCT molecule_id 
            FROM molecular_properties 
            WHERE property_type = 'Molecular Weight'
        )
    """
    
    if limit:
        query += f" LIMIT {limit}"
    
    with conn.cursor(cursor_factory=RealDictCursor) as cursor:
        cursor.execute(query)
        molecules = cursor.fetchall()
    
    logger.info(f"Found {len(molecules)} molecules needing property calculations")
    return molecules

def get_property_type_id_map(conn):
    """Get a mapping of property type names to their IDs."""
    property_map = {}
    
    with conn.cursor() as cursor:
        cursor.execute("SELECT id, name FROM property_types")
        for row in cursor.fetchall():
            property_id, property_name = row
            property_map[property_name] = property_id
    
    return property_map

def calculate_properties(smiles):
    """Calculate molecular properties from SMILES."""
    if not smiles:
        return None
    
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
            'Heavy Atom Count': mock_rdkit_formula.calculate_heavy_atom_count(smiles)
        }
        
        return properties
    except Exception as e:
        logger.error(f"Error calculating properties: {e}")
        return None

def property_exists(conn, molecule_id, property_type_id):
    """Check if a property already exists for a molecule."""
    with conn.cursor() as cursor:
        cursor.execute("""
            SELECT COUNT(*) FROM molecular_properties
            WHERE molecule_id = %s AND property_type_id = %s
        """, (molecule_id, property_type_id))
        count = cursor.fetchone()[0]
        return count > 0

def add_property(conn, molecule_id, property_type, property_type_id, value, source="mock_rdkit_formula"):
    """Add a property to the molecular_properties table."""
    # Determine value type (numeric or text)
    if isinstance(value, (int, float)):
        numeric_value = value
        text_value = None
    else:
        numeric_value = None
        text_value = str(value)
    
    # Check if property already exists
    exists = property_exists(conn, molecule_id, property_type_id)
    
    with conn.cursor() as cursor:
        try:
            if exists:
                # Update existing property
                cursor.execute("""
                    UPDATE molecular_properties
                    SET 
                        numeric_value = %s,
                        text_value = %s,
                        source = %s,
                        updated_at = NOW()
                    WHERE molecule_id = %s AND property_type_id = %s
                """, (
                    numeric_value,
                    text_value,
                    source,
                    molecule_id,
                    property_type_id
                ))
            else:
                # Insert new property
                cursor.execute("""
                    INSERT INTO molecular_properties
                    (id, molecule_id, property_type_id, property_type, numeric_value, text_value, source, created_at, updated_at)
                    VALUES (uuid_generate_v4(), %s, %s, %s, %s, %s, %s, NOW(), NOW())
                """, (
                    molecule_id,
                    property_type_id,
                    property_type,
                    numeric_value,
                    text_value,
                    source
                ))
            
            # Commit immediately for each property
            conn.commit()
            return True
        except Exception as e:
            conn.rollback()
            logger.error(f"Error adding property {property_type} to molecule {molecule_id}: {e}")
            return False

def update_molecule_metadata(conn, molecule_id):
    """Update molecule's metadata to indicate properties have been added."""
    try:
        with conn.cursor() as cursor:
            # Check if molecules table has a metadata column
            cursor.execute("""
                SELECT column_name 
                FROM information_schema.columns 
                WHERE table_name = 'molecules' AND column_name = 'metadata'
            """)
            has_metadata = cursor.fetchone() is not None
            
            if has_metadata:
                cursor.execute("""
                    UPDATE molecules
                    SET 
                        updated_at = NOW(),
                        metadata = COALESCE(metadata, '{}'::jsonb) || %s::jsonb
                    WHERE id = %s
                """, (
                    json.dumps({
                        "properties_completed": True,
                        "properties_completed_at": time.strftime("%Y-%m-%d %H:%M:%S")
                    }),
                    molecule_id
                ))
            else:
                # Just update the timestamp if no metadata column
                cursor.execute("""
                    UPDATE molecules
                    SET updated_at = NOW()
                    WHERE id = %s
                """, (molecule_id,))
            
            conn.commit()
            return cursor.rowcount > 0
    except Exception as e:
        conn.rollback()
        logger.error(f"Error updating molecule metadata for {molecule_id}: {e}")
        return False

def main():
    """Main function to complete missing properties."""
    import argparse
    parser = argparse.ArgumentParser(description="Complete missing molecular properties")
    parser.add_argument("--limit", type=int, default=None, help="Limit the number of molecules to process")
    parser.add_argument("--batch-size", type=int, default=50, help="Number of molecules to process in each batch")
    parser.add_argument("--dry-run", action="store_true", help="Don't actually update the database")
    args = parser.parse_args()
    
    if args.dry_run:
        logger.info("Running in dry run mode - no data will be written to the database")
    
    # Connect to the database
    conn = connect_to_db()
    
    try:
        # Get property type ID mapping
        property_type_id_map = get_property_type_id_map(conn)
        logger.info(f"Loaded {len(property_type_id_map)} property types from database")
        
        # Check that we have all needed property type IDs
        needed_properties = [
            'Molecular Weight', 'LogP', 'TPSA', 
            'Hydrogen Bond Donor Count', 'Hydrogen Bond Acceptor Count', 
            'Rotatable Bond Count', 'Ring Count', 'Aromatic Ring Count', 'Heavy Atom Count'
        ]
        missing_property_types = [prop for prop in needed_properties if prop not in property_type_id_map]
        
        if missing_property_types:
            logger.error(f"Missing property types in database: {missing_property_types}")
            return
        
        # Get molecules needing properties
        molecules = get_molecules_needing_properties(conn, args.limit)
        
        if not molecules:
            logger.info("No molecules needing property calculations found. Exiting.")
            return
        
        # Process molecules in batches
        total_updated = 0
        batch_count = (len(molecules) + args.batch_size - 1) // args.batch_size
        
        for batch_idx in range(batch_count):
            start_idx = batch_idx * args.batch_size
            end_idx = min(start_idx + args.batch_size, len(molecules))
            batch = molecules[start_idx:end_idx]
            
            logger.info(f"Processing batch {batch_idx+1}/{batch_count} ({len(batch)} molecules)")
            
            batch_updated = 0
            
            for molecule in batch:
                molecule_id = molecule['id']
                smiles = molecule['smiles']
                name = molecule['name']
                
                # Calculate properties
                properties = calculate_properties(smiles)
                
                if properties:
                    if not args.dry_run:
                        # Add each property to the database
                        property_added_count = 0
                        try:
                            for prop_type, value in properties.items():
                                if value is not None and prop_type in property_type_id_map:
                                    property_type_id = property_type_id_map[prop_type]
                                    success = add_property(conn, molecule_id, prop_type, property_type_id, value)
                                    if success:
                                        property_added_count += 1
                            
                            # Update molecule metadata if any properties were added
                            if property_added_count > 0:
                                update_molecule_metadata(conn, molecule_id)
                                batch_updated += 1
                                
                                # Log every 10th update to avoid too much output
                                if batch_updated % 10 == 0 or batch_updated == 1:
                                    logger.info(f"Added {property_added_count} properties to {name}")
                        except Exception as e:
                            logger.error(f"Error processing molecule {name} ({molecule_id}): {e}")
                            # Continue with the next molecule even if this one fails
                    else:
                        batch_updated += 1
                        
                        # Log every 10th update to avoid too much output
                        if batch_updated % 10 == 0 or batch_updated == 1:
                            logger.info(f"Would add properties to {name}: {properties}")
            
            # No need to commit batch - we're committing per property
            total_updated += batch_updated
            if not args.dry_run:
                logger.info(f"Added property updates for {batch_updated} molecules in batch {batch_idx+1}")
            else:
                logger.info(f"Would have updated {batch_updated} molecules in batch {batch_idx+1}")
            
            # Brief pause between batches
            if batch_idx < batch_count - 1:
                time.sleep(1)
        
        # Final summary
        logger.info(f"Property update complete. {'Would have updated' if args.dry_run else 'Updated'} {total_updated} molecules.")
        
    except Exception as e:
        conn.rollback()
        logger.error(f"Error completing properties: {e}", exc_info=True)
        sys.exit(1)
    
    finally:
        conn.close()

if __name__ == "__main__":
    main()