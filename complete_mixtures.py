#!/usr/bin/env python3
"""
Complete or remove placeholder mixtures in the database.

This script analyzes existing mixtures without components and either
completes them with proper components or removes them if appropriate.
"""

import os
import sys
import argparse
import psycopg2
from psycopg2.extras import RealDictCursor
import json
import uuid
from datetime import datetime

# Known cryoprotectant mixture compositions
MIXTURE_COMPOSITIONS = {
    'VS55 Vitrification Solution': [
        {'molecule': 'dimethyl sulfoxide', 'concentration': 8.4, 'unit': '%w/v'},
        {'molecule': 'formamide', 'concentration': 22.3, 'unit': '%w/v'},
        {'molecule': 'propylene glycol', 'concentration': 16.8, 'unit': '%w/v'}
    ],
    'DP6 Vitrification Solution': [
        {'molecule': 'dimethyl sulfoxide', 'concentration': 23.0, 'unit': '%w/v'},
        {'molecule': 'propylene glycol', 'concentration': 16.0, 'unit': '%w/v'}
    ],
    'M22 Vitrification Solution': [
        {'molecule': 'dimethyl sulfoxide', 'concentration': 22.3, 'unit': '%w/v'},
        {'molecule': 'formamide', 'concentration': 12.8, 'unit': '%w/v'},
        {'molecule': 'propylene glycol', 'concentration': 8.4, 'unit': '%w/v'},
        {'molecule': 'ethylene glycol', 'concentration': 9.9, 'unit': '%w/v'},
        {'molecule': 'n-methylformamide', 'concentration': 1.7, 'unit': '%w/v'},
        {'molecule': '3-methoxy-1,2-propanediol', 'concentration': 4.9, 'unit': '%w/v'}
    ],
    'EAFS10/10 Solution': [
        {'molecule': 'ethylene glycol', 'concentration': 10.0, 'unit': '%w/v'},
        {'molecule': 'acetamide', 'concentration': 10.0, 'unit': '%w/v'},
        {'molecule': 'formamide', 'concentration': 10.0, 'unit': '%w/v'},
        {'molecule': 'dimethyl sulfoxide', 'concentration': 10.0, 'unit': '%w/v'}
    ],
    'Glycerol/Trehalose Solution': [
        {'molecule': 'glycerol', 'concentration': 15.0, 'unit': '%w/v'},
        {'molecule': 'trehalose', 'concentration': 15.0, 'unit': '%w/v'}
    ],
    'Betaine/Glycerol Solution': [
        {'molecule': 'glycerol', 'concentration': 10.0, 'unit': '%w/v'},
        {'molecule': 'betaine', 'concentration': 5.0, 'unit': '%w/v'}
    ],
    'DMSO/Proline Solution': [
        {'molecule': 'dimethyl sulfoxide', 'concentration': 10.0, 'unit': '%w/v'},
        {'molecule': 'proline', 'concentration': 5.0, 'unit': '%w/v'}
    ],
    'Hydroxyectoine/DMSO Solution': [
        {'molecule': 'dimethyl sulfoxide', 'concentration': 7.5, 'unit': '%w/v'},
        {'molecule': 'hydroxyectoine', 'concentration': 5.0, 'unit': '%w/v'}
    ]
}

# Should keep even if no known composition
KEEP_MIXTURES = [
    'VS55 Vitrification Solution',  # Already has components
    'DP6 Vitrification Solution',
    'M22 Vitrification Solution',
    'EAFS10/10 Solution',
    'Glycerol/Trehalose Solution',
    'Betaine/Glycerol Solution',
    'DMSO/Proline Solution',
    'Hydroxyectoine/DMSO Solution'
]

def connect_to_db():
    """Connect to the database using environment variables."""
    from dotenv import load_dotenv

    # Load environment variables from .env file
    load_dotenv()

    db_params = {
        'host': os.getenv('SUPABASE_DB_HOST'),
        'port': os.getenv('SUPABASE_DB_PORT', '5432'),
        'dbname': os.getenv('SUPABASE_DB_NAME', 'postgres'),
        'user': os.getenv('SUPABASE_DB_USER'),
        'password': os.getenv('SUPABASE_DB_PASSWORD'),
        'sslmode': os.getenv('SUPABASE_DB_SSLMODE', 'require')
    }

    # Ensure required parameters are present
    if not all([db_params['host'], db_params['user'], db_params['password']]):
        print("Error: Missing required database connection parameters in environment variables.")
        print("Make sure you have a valid .env file with SUPABASE_DB_* variables.")
        sys.exit(1)

    try:
        conn = psycopg2.connect(**db_params)
        return conn
    except psycopg2.Error as e:
        print(f"Database connection error: {e}")
        sys.exit(1)

def get_mixture_status(conn):
    """Get status of all mixtures including component counts."""
    with conn.cursor(cursor_factory=RealDictCursor) as cursor:
        cursor.execute("""
            SELECT 
                m.id, 
                m.name, 
                m.description,
                COUNT(mc.id) as component_count
            FROM 
                mixtures m
            LEFT JOIN 
                mixture_components mc ON m.id = mc.mixture_id
            GROUP BY 
                m.id, m.name, m.description
            ORDER BY 
                component_count DESC, m.name
        """)
        
        return cursor.fetchall()

def find_molecule_by_name(conn, name):
    """Find a molecule by name or similar name."""
    with conn.cursor(cursor_factory=RealDictCursor) as cursor:
        # Try exact match first
        cursor.execute("""
            SELECT id, name, smiles
            FROM molecules
            WHERE LOWER(name) = LOWER(%s)
            LIMIT 1
        """, (name,))
        
        result = cursor.fetchone()
        if result:
            return result
        
        # Try partial match
        cursor.execute("""
            SELECT id, name, smiles
            FROM molecules
            WHERE LOWER(name) LIKE LOWER(%s)
            LIMIT 1
        """, (f"%{name}%",))
        
        return cursor.fetchone()

def complete_mixture_components(conn, mixture, dry_run=False):
    """Add components to a mixture based on known compositions."""
    print(f"\nProcessing mixture: {mixture['name']}")
    
    # Check if we have a known composition for this mixture
    if mixture['name'] not in MIXTURE_COMPOSITIONS:
        print(f"  No known composition for '{mixture['name']}'")
        return False
    
    composition = MIXTURE_COMPOSITIONS[mixture['name']]
    print(f"  Found composition with {len(composition)} components")
    
    # Track molecules that were found and added
    found_molecules = []
    not_found_molecules = []
    
    # For each component, find the molecule and add it
    for component in composition:
        molecule_name = component['molecule']
        molecule = find_molecule_by_name(conn, molecule_name)
        
        if molecule:
            found_molecules.append({
                'molecule': molecule,
                'concentration': component['concentration'],
                'unit': component['unit']
            })
            print(f"  Found molecule '{molecule_name}' (ID: {molecule['id']})")
        else:
            not_found_molecules.append(molecule_name)
            print(f"  Could not find molecule '{molecule_name}'")
    
    # If we found all molecules, add them as components
    if len(not_found_molecules) == 0:
        if not dry_run:
            with conn.cursor() as cursor:
                for component in found_molecules:
                    # Create a new component
                    component_id = str(uuid.uuid4())
                    cursor.execute("""
                        INSERT INTO mixture_components (
                            id, mixture_id, molecule_id, concentration, unit, 
                            created_at, updated_at
                        ) VALUES (
                            %s, %s, %s, %s, %s, NOW(), NOW()
                        )
                    """, (
                        component_id,
                        mixture['id'],
                        component['molecule']['id'],
                        component['concentration'],
                        component['unit']
                    ))
                
                conn.commit()
                
            print(f"  Added {len(found_molecules)} components to mixture '{mixture['name']}'")
            return True
        else:
            print(f"  [DRY RUN] Would add {len(found_molecules)} components to mixture '{mixture['name']}'")
            return True
    else:
        print(f"  Could not complete mixture '{mixture['name']}' - missing {len(not_found_molecules)} molecules")
        return False

def remove_placeholder_mixture(conn, mixture, dry_run=False):
    """Remove a placeholder mixture with no components."""
    print(f"\nRemoving placeholder mixture: {mixture['name']}")
    
    # Check for references to this mixture
    with conn.cursor(cursor_factory=RealDictCursor) as cursor:
        cursor.execute("""
            SELECT 
                'experiments' as table_name, COUNT(*) as ref_count 
            FROM 
                experiments 
            WHERE 
                mixture_id = %s
            UNION ALL
            SELECT 
                'predictions' as table_name, COUNT(*) as ref_count 
            FROM 
                predictions 
            WHERE 
                mixture_id = %s
        """, (mixture['id'], mixture['id']))
        
        references = cursor.fetchall()
        
        has_references = False
        for ref in references:
            if ref['ref_count'] > 0:
                has_references = True
                print(f"  Cannot remove: Referenced by {ref['ref_count']} rows in {ref['table_name']}")
        
        if has_references:
            return False
        
        if not dry_run:
            cursor.execute("""
                DELETE FROM mixtures
                WHERE id = %s
            """, (mixture['id'],))
            
            conn.commit()
            
            print(f"  Removed mixture '{mixture['name']}'")
            return True
        else:
            print(f"  [DRY RUN] Would remove mixture '{mixture['name']}'")
            return True

def create_audit_log(conn, completed, removed, dry_run=False):
    """Create an audit log entry for the mixture updates."""
    if dry_run:
        print("\n[DRY RUN] Would create audit log entry")
        return

    try:
        with conn.cursor() as cursor:
            cursor.execute("""
                INSERT INTO scientific_data_audit (
                    table_name,
                    operation,
                    user_id,
                    timestamp,
                    old_data,
                    new_data,
                    application_context
                ) VALUES (
                    'mixtures',
                    'UPDATE',
                    NULL,
                    NOW(),
                    NULL,
                    %s::jsonb,
                    'database_remediation_script'
                )
                RETURNING id
            """, (
                json.dumps({
                    'operation': 'complete_mixtures',
                    'timestamp': datetime.now().isoformat(),
                    'details': f"Completed {len(completed)} mixtures, removed {len(removed)} placeholders",
                    'completed': [m['name'] for m in completed],
                    'removed': [m['name'] for m in removed]
                }),
            ))

            audit_id = cursor.fetchone()[0]
            conn.commit()

            print(f"\nCreated audit log entry with ID: {audit_id}")
    except Exception as e:
        print(f"Error creating audit log: {e}")
        print("The data changes were applied, but the audit log could not be created.")
        print("This might be due to incompatibility with the existing audit table structure.")

def process_mixtures(conn, dry_run=False):
    """Process all mixtures, completing or removing as appropriate."""
    # Get current mixture status
    mixtures = get_mixture_status(conn)
    
    print(f"\nFound {len(mixtures)} mixtures in the database")
    for mixture in mixtures:
        print(f"  {mixture['name']}: {mixture['component_count']} components")
    
    # Track completed and removed mixtures
    completed_mixtures = []
    removed_mixtures = []
    
    # Process each mixture
    for mixture in mixtures:
        # Skip mixtures that already have components
        if mixture['component_count'] > 0:
            print(f"\nSkipping {mixture['name']} - already has {mixture['component_count']} components")
            continue
        
        # Try to complete the mixture
        if mixture['name'] in MIXTURE_COMPOSITIONS:
            if complete_mixture_components(conn, mixture, dry_run):
                completed_mixtures.append(mixture)
        elif mixture['name'] not in KEEP_MIXTURES:
            # Remove placeholder mixtures that we don't need to keep
            if remove_placeholder_mixture(conn, mixture, dry_run):
                removed_mixtures.append(mixture)
    
    # Create audit log
    if not dry_run and (completed_mixtures or removed_mixtures):
        create_audit_log(conn, completed_mixtures, removed_mixtures)
    
    # Print summary
    print("\nMixture Processing Summary:")
    print(f"  Completed {len(completed_mixtures)} mixtures")
    print(f"  Removed {len(removed_mixtures)} placeholder mixtures")
    print(f"  Skipped {len(mixtures) - len(completed_mixtures) - len(removed_mixtures)} mixtures")

def main():
    parser = argparse.ArgumentParser(description="Complete or remove placeholder mixtures")
    parser.add_argument("--dry-run", action="store_true",
                        help="Show what would be done without making changes")
    parser.add_argument("--yes", action="store_true",
                        help="Automatically confirm changes without prompting")

    args = parser.parse_args()

    print(f"Starting mixture completion process at {datetime.now().isoformat()}")
    if args.dry_run:
        print("DRY RUN MODE: No changes will be made to the database")
    if args.yes:
        print("AUTO-CONFIRM MODE: Changes will be applied without confirmation")

    # Connect to the database using environment variables
    conn = connect_to_db()
    
    try:
        # Process mixtures
        process_mixtures(conn, args.dry_run)
        
    except Exception as e:
        print(f"Error during mixture processing: {e}")
        conn.rollback()
        sys.exit(1)
    
    finally:
        # Close the connection
        if conn:
            conn.close()
    
    print(f"Mixture completion process finished at {datetime.now().isoformat()}")

if __name__ == "__main__":
    main()