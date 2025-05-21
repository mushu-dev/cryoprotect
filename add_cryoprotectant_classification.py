#!/usr/bin/env python3
"""
Add classification fields for cryoprotectants.

This script adds is_known_cryoprotectant and cryoprotectant_type fields
to the molecules table and populates them for known cryoprotectants.
"""

import os
import sys
import argparse
import psycopg2
from psycopg2.extras import RealDictCursor
import json
from datetime import datetime

# Known cryoprotectant molecules and their classifications
KNOWN_CRYOPROTECTANTS = {
    # Penetrating cryoprotectants
    'dimethyl sulfoxide': {'type': 'PENETRATING', 'alias': ['dmso']},
    'glycerol': {'type': 'PENETRATING', 'alias': []},
    'propylene glycol': {'type': 'PENETRATING', 'alias': ['1,2-propanediol']},
    'ethylene glycol': {'type': 'PENETRATING', 'alias': ['1,2-ethanediol']},
    'methanol': {'type': 'PENETRATING', 'alias': []},
    'formamide': {'type': 'PENETRATING', 'alias': []},
    'n-methylformamide': {'type': 'PENETRATING', 'alias': []},
    '3-methoxy-1,2-propanediol': {'type': 'PENETRATING', 'alias': []},
    'acetamide': {'type': 'PENETRATING', 'alias': []},
    'propanediol': {'type': 'PENETRATING', 'alias': ['1,3-propanediol']},
    'butanediol': {'type': 'PENETRATING', 'alias': ['1,4-butanediol']},
    
    # Non-penetrating cryoprotectants
    'trehalose': {'type': 'NON_PENETRATING', 'alias': []},
    'sucrose': {'type': 'NON_PENETRATING', 'alias': []},
    'glucose': {'type': 'NON_PENETRATING', 'alias': []},
    'maltose': {'type': 'NON_PENETRATING', 'alias': []},
    'raffinose': {'type': 'NON_PENETRATING', 'alias': []},
    'dextran': {'type': 'NON_PENETRATING', 'alias': []},
    'polyvinylpyrrolidone': {'type': 'NON_PENETRATING', 'alias': ['pvp']},
    'polyethylene glycol': {'type': 'NON_PENETRATING', 'alias': ['peg']},
    'mannitol': {'type': 'NON_PENETRATING', 'alias': []},
    'sorbitol': {'type': 'NON_PENETRATING', 'alias': []},
    'xylitol': {'type': 'NON_PENETRATING', 'alias': []},
    'inositol': {'type': 'NON_PENETRATING', 'alias': []},
    
    # Amino acid cryoprotectants
    'proline': {'type': 'AMINO_ACID', 'alias': []},
    'glycine': {'type': 'AMINO_ACID', 'alias': []},
    'alanine': {'type': 'AMINO_ACID', 'alias': []},
    'serine': {'type': 'AMINO_ACID', 'alias': []},
    'glutamate': {'type': 'AMINO_ACID', 'alias': ['glutamic acid']},
    'glutamine': {'type': 'AMINO_ACID', 'alias': []},
    'valine': {'type': 'AMINO_ACID', 'alias': []},
    'lysine': {'type': 'AMINO_ACID', 'alias': []},
    'arginine': {'type': 'AMINO_ACID', 'alias': []},
    'histidine': {'type': 'AMINO_ACID', 'alias': []},
    'tryptophan': {'type': 'AMINO_ACID', 'alias': []},
    
    # Osmolytes and other cryoprotectants
    'trimethylamine n-oxide': {'type': 'OSMOLYTE', 'alias': ['tmao']},
    'betaine': {'type': 'OSMOLYTE', 'alias': []},
    'ectoine': {'type': 'OSMOLYTE', 'alias': []},
    'hydroxyectoine': {'type': 'OSMOLYTE', 'alias': []},
    'meso-erythritol': {'type': 'POLYOL', 'alias': ['erythritol']},
    'hydroxyethyl starch': {'type': 'POLYMER', 'alias': ['hes']},
}

def connect_to_db(db_params):
    """Connect to the database and return a connection."""
    try:
        conn = psycopg2.connect(**db_params)
        return conn
    except psycopg2.Error as e:
        print(f"Database connection error: {e}")
        sys.exit(1)

def add_classification_fields(conn, dry_run=False):
    """Add classification fields to molecules table if they don't exist."""
    with conn.cursor() as cursor:
        # Check if fields already exist
        cursor.execute("""
            SELECT 
                COUNT(*) 
            FROM 
                information_schema.columns 
            WHERE 
                table_schema = 'public' 
                AND table_name = 'molecules' 
                AND column_name = 'is_known_cryoprotectant'
        """)
        
        if cursor.fetchone()[0] == 0:
            # Add the fields
            if not dry_run:
                print("Adding classification fields to molecules table...")
                cursor.execute("""
                    ALTER TABLE molecules 
                    ADD COLUMN is_known_cryoprotectant BOOLEAN DEFAULT false,
                    ADD COLUMN cryoprotectant_type VARCHAR(50)
                """)
                conn.commit()
                print("Fields added successfully.")
            else:
                print("[DRY RUN] Would add classification fields to molecules table")
        else:
            print("Classification fields already exist in molecules table.")

def classify_known_cryoprotectants(conn, dry_run=False):
    """Classify known cryoprotectant molecules."""
    # First, get all molecules
    with conn.cursor(cursor_factory=RealDictCursor) as cursor:
        cursor.execute("""
            SELECT id, name, LOWER(name) as lower_name
            FROM molecules
        """)
        
        all_molecules = cursor.fetchall()
        
        # Create a dictionary mapping lowercase names to molecule IDs
        name_to_id = {}
        for molecule in all_molecules:
            name_to_id[molecule['lower_name']] = molecule['id']
            
            # Check if name contains any alias
            for cp_name, cp_info in KNOWN_CRYOPROTECTANTS.items():
                for alias in cp_info['alias']:
                    if alias.lower() in molecule['lower_name']:
                        name_to_id[alias.lower()] = molecule['id']
        
        # Track matches
        matches = []
        
        # Classify known cryoprotectants
        for cp_name, cp_info in KNOWN_CRYOPROTECTANTS.items():
            found = False
            
            # Check exact match
            if cp_name.lower() in name_to_id:
                molecule_id = name_to_id[cp_name.lower()]
                matches.append({
                    'id': molecule_id,
                    'name': cp_name,
                    'type': cp_info['type'],
                    'match_type': 'exact'
                })
                found = True
            
            # Check aliases
            if not found:
                for alias in cp_info['alias']:
                    if alias.lower() in name_to_id:
                        molecule_id = name_to_id[alias.lower()]
                        matches.append({
                            'id': molecule_id,
                            'name': alias,
                            'type': cp_info['type'],
                            'match_type': 'alias'
                        })
                        found = True
                        break
            
            # Check partial match (if no exact or alias match)
            if not found:
                for mol_name, mol_id in name_to_id.items():
                    if cp_name.lower() in mol_name or any(alias.lower() in mol_name for alias in cp_info['alias']):
                        matches.append({
                            'id': mol_id,
                            'name': cp_name,
                            'type': cp_info['type'],
                            'match_type': 'partial'
                        })
                        break
        
        # Update the database
        if not dry_run:
            print(f"\nUpdating {len(matches)} molecules with classification information...")
            for match in matches:
                cursor.execute("""
                    UPDATE molecules
                    SET 
                        is_known_cryoprotectant = true,
                        cryoprotectant_type = %s,
                        updated_at = NOW()
                    WHERE id = %s
                """, (match['type'], match['id']))
            
            conn.commit()
            print("Classification complete.")
        else:
            print(f"\n[DRY RUN] Would update {len(matches)} molecules with classification information")
        
        # Print matches by type
        match_types = {}
        for match in matches:
            if match['type'] not in match_types:
                match_types[match['type']] = []
            match_types[match['type']].append(match)
        
        print("\nMatched cryoprotectants by type:")
        for type_name, type_matches in match_types.items():
            print(f"  {type_name}: {len(type_matches)} molecules")
            for match in type_matches[:5]:  # Show only first 5
                print(f"    - {match['name']} (match type: {match['match_type']})")
            if len(type_matches) > 5:
                print(f"    - ... and {len(type_matches) - 5} more")
        
        return matches

def create_audit_log(conn, matches, dry_run=False):
    """Create an audit log entry for the classification."""
    if dry_run:
        return
    
    with conn.cursor() as cursor:
        cursor.execute("""
            INSERT INTO scientific_data_audit (
                operation_type,
                table_name,
                operation_details,
                performed_by,
                performed_at,
                affected_rows,
                operation_status,
                operation_metadata
            ) VALUES (
                'UPDATE',
                'molecules',
                'Add cryoprotectant classification',
                'database_remediation_script',
                NOW(),
                %s,
                'COMPLETED',
                %s
            )
            RETURNING id
        """, (
            len(matches),
            json.dumps({
                'operation': 'add_cryoprotectant_classification',
                'timestamp': datetime.now().isoformat(),
                'details': f"Classified {len(matches)} molecules as known cryoprotectants",
                'types': {
                    type_name: len([m for m in matches if m['type'] == type_name])
                    for type_name in set(m['type'] for m in matches)
                }
            })
        ))
        
        audit_id = cursor.fetchone()[0]
        conn.commit()
        
        print(f"\nCreated audit log entry with ID: {audit_id}")

def main():
    parser = argparse.ArgumentParser(description="Add cryoprotectant classification fields")
    parser.add_argument("--host", default=os.environ.get("DB_HOST", "aws-0-us-east-1.pooler.supabase.com"), 
                        help="Database host")
    parser.add_argument("--port", default=os.environ.get("DB_PORT", "5432"), 
                        help="Database port")
    parser.add_argument("--dbname", default=os.environ.get("DB_NAME", "postgres"), 
                        help="Database name")
    parser.add_argument("--user", default=os.environ.get("DB_USER", "postgres.tsdlmynydfuypiugmkev"), 
                        help="Database user")
    parser.add_argument("--password", default=os.environ.get("DB_PASSWORD"), 
                        help="Database password")
    parser.add_argument("--dry-run", action="store_true", 
                        help="Show what would be done without making changes")
    
    args = parser.parse_args()
    
    # Prepare database connection parameters
    db_params = {
        'host': args.host,
        'port': args.port,
        'dbname': args.dbname,
        'user': args.user,
        'password': args.password,
        'sslmode': 'require'
    }
    
    # Ensure password is provided
    if not db_params['password']:
        print("Error: Database password is required. Set DB_PASSWORD environment variable or use --password.")
        sys.exit(1)
    
    print(f"Starting cryoprotectant classification process at {datetime.now().isoformat()}")
    if args.dry_run:
        print("DRY RUN MODE: No changes will be made to the database")
    
    # Connect to the database
    conn = connect_to_db(db_params)
    
    try:
        # Add classification fields
        add_classification_fields(conn, args.dry_run)
        
        # Classify known cryoprotectants
        matches = classify_known_cryoprotectants(conn, args.dry_run)
        
        # Create audit log
        create_audit_log(conn, matches, args.dry_run)
        
    except Exception as e:
        print(f"Error during cryoprotectant classification: {e}")
        conn.rollback()
        sys.exit(1)
    
    finally:
        # Close the connection
        if conn:
            conn.close()
    
    print(f"Cryoprotectant classification process completed at {datetime.now().isoformat()}")

if __name__ == "__main__":
    main()