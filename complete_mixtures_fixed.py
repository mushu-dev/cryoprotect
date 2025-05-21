#!/usr/bin/env python3
"""
This script completes mixture data by either:
1. Adding proper components and composition for known mixtures
2. Removing placeholder mixtures with no scientific value

Part of the CryoProtect database remediation plan.
"""

import os
import sys
import argparse
import psycopg2
from psycopg2.extras import RealDictCursor
import json
import uuid
from datetime import datetime
from dotenv import load_dotenv

# Load environment variables
load_dotenv()

# Known mixture compositions for common cryoprotectant mixtures
KNOWN_MIXTURES = {
    "VS55 Vitrification Solution": {
        "components": [
            {"name": "dimethyl sulfoxide", "formula": "C2H6OS", "concentration": 8.4, "unit": "%w/v"},
            {"name": "propylene glycol", "formula": "C3H8O2", "concentration": 16.84, "unit": "%w/v"},
            {"name": "formamide", "formula": "CH3NO", "concentration": 16.84, "unit": "%w/v"}
        ],
        "description": "VS55 is a well-known vitrification solution used for organ cryopreservation.",
        "reference": "Brockbank, K.G., et al. (2010). Comparison of VS55 and DP6 for successful rabbit kidney vitrification. Cryobiology, 61, 327-332."
    },
    "M22 Vitrification Solution": {
        "components": [
            {"name": "dimethyl sulfoxide", "formula": "C2H6OS", "concentration": 22.3, "unit": "%w/v"},
            {"name": "formamide", "formula": "CH3NO", "concentration": 12.8, "unit": "%w/v"},
            {"name": "ethylene glycol", "formula": "C2H6O2", "concentration": 16.8, "unit": "%w/v"}
        ],
        "description": "M22 is an advanced vitrification solution designed for minimal toxicity.",
        "reference": "Fahy, G.M., et al. (2009). Physical and biological aspects of renal vitrification. Organogenesis, 5, 167-175."
    },
    "DP6 Vitrification Solution": {
        "components": [
            {"name": "dimethyl sulfoxide", "formula": "C2H6OS", "concentration": 6, "unit": "%w/v"},
            {"name": "propylene glycol", "formula": "C3H8O2", "concentration": 6, "unit": "%w/v"}
        ],
        "description": "DP6 is a dimethyl sulfoxide and propylene glycol mixture used in cryopreservation.",
        "reference": "Taylor, M.J., et al. (1995). A new solution for life without blood. Asanguineous low flow perfusion of a whole-body perfusate during 3 hours of cardiac arrest and profound hypothermia. Circulation, 91, 431-444."
    },
    "EAFS10/10 Solution": {
        "components": [
            {"name": "ethylene glycol", "formula": "C2H6O2", "concentration": 10, "unit": "%w/v"},
            {"name": "acetamide", "formula": "C2H5NO", "concentration": 10, "unit": "%w/v"},
            {"name": "formamide", "formula": "CH3NO", "concentration": 10, "unit": "%w/v"},
            {"name": "dimethyl sulfoxide", "formula": "C2H6OS", "concentration": 10, "unit": "%w/v"}
        ],
        "description": "EAFS10/10 is a vitrification solution for embryo cryopreservation.",
        "reference": "Pedro, P.B., et al. (2005). Permeability of mouse oocytes and embryos at various developmental stages to five cryoprotectants. Journal of Reproduction and Development, 51, 235-246."
    },
    "Glycerol/Trehalose Solution": {
        "components": [
            {"name": "glycerol", "formula": "C3H8O3", "concentration": 15, "unit": "%w/v"},
            {"name": "trehalose", "formula": "C12H22O11", "concentration": 15, "unit": "%w/v"}
        ],
        "description": "Glycerol/Trehalose combination for cell cryopreservation.",
        "reference": "Hubálek, Z. (2003). Protectants used in the cryopreservation of microorganisms. Cryobiology, 46, 205-229."
    },
    "Betaine/Glycerol Solution": {
        "components": [
            {"name": "glycerol", "formula": "C3H8O3", "concentration": 10, "unit": "%w/v"},
            {"name": "betaine", "formula": "C5H11NO2", "concentration": 5, "unit": "%w/v"}
        ],
        "description": "Betaine/Glycerol combination for improved cryopreservation.",
        "reference": "Dong, Q., et al. (2009). Improvement of post-thaw sperm motility using betaine. Cryobiology, 59, 369-372."
    },
    "DMSO/Proline Solution": {
        "components": [
            {"name": "dimethyl sulfoxide", "formula": "C2H6OS", "concentration": 10, "unit": "%w/v"},
            {"name": "proline", "formula": "C5H9NO2", "concentration": 5, "unit": "%w/v"}
        ],
        "description": "DMSO combined with proline for improved cryopreservation outcomes.",
        "reference": "Withers, L.A., King, P.J. (1979). Proline: A novel cryoprotectant for the freeze preservation of cultured cells of Zea mays L. Plant Physiology, 64, 675-678."
    },
    "Hydroxyectoine/DMSO Solution": {
        "components": [
            {"name": "dimethyl sulfoxide", "formula": "C2H6OS", "concentration": 7.5, "unit": "%w/v"},
            {"name": "hydroxyectoine", "formula": "C7H12N2O3", "concentration": 5, "unit": "%w/v"}
        ],
        "description": "Combination of DMSO with hydroxyectoine for improved cryopreservation.",
        "reference": "Bursy, J., et al. (2008). Synthesis and uptake of the compatible solutes ectoine and hydroxyectoine by Streptomyces coelicolor A3(2) in response to salt and heat stresses. Applied and Environmental Microbiology, 74, 7286-7296."
    }
}

# Mixtures to keep even if we can't add components yet
KEEP_MIXTURES = list(KNOWN_MIXTURES.keys())

def connect_to_db():
    """Connect to the database using environment variables."""
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

def get_placeholder_mixtures(conn):
    """Get mixtures that are potentially placeholders (no components)."""
    with conn.cursor(cursor_factory=RealDictCursor) as cursor:
        cursor.execute("""
            SELECT m.id, m.name, m.description
            FROM mixtures m
            LEFT JOIN mixture_components mc ON m.id = mc.mixture_id
            WHERE mc.id IS NULL
            ORDER BY m.name
        """)
        
        return cursor.fetchall()

def get_known_mixtures(conn):
    """Get mixtures that match our known cryoprotectant mixtures by name."""
    with conn.cursor(cursor_factory=RealDictCursor) as cursor:
        # Create an SQL IN clause with all known mixture names
        known_names = list(KNOWN_MIXTURES.keys())
        placeholders = ', '.join(['%s'] * len(known_names))
        
        query = f"""
            SELECT m.id, m.name, m.description, 
                   EXISTS (
                       SELECT 1 FROM mixture_components mc 
                       WHERE mc.mixture_id = m.id
                   ) AS has_components
            FROM mixtures m
            WHERE m.name IN ({placeholders})
            ORDER BY m.name
        """
        
        cursor.execute(query, known_names)
        
        return cursor.fetchall()

def get_molecule_by_name(conn, name):
    """Find a molecule by name (fuzzy match)."""
    with conn.cursor(cursor_factory=RealDictCursor) as cursor:
        cursor.execute("""
            SELECT id, name, formula 
            FROM molecules
            WHERE LOWER(name) LIKE LOWER(%s)
            ORDER BY LENGTH(name) ASC
            LIMIT 1
        """, [f"%{name}%"])
        
        result = cursor.fetchone()
        if not result:
            print(f"  Warning: Component '{name}' not found in database, will be created")
        return result

def create_molecule(conn, name, formula=None):
    """Create a new molecule record."""
    with conn.cursor(cursor_factory=RealDictCursor) as cursor:
        molecule_id = str(uuid.uuid4())
        cursor.execute("""
            INSERT INTO molecules (
                id, name, formula, data_source, created_at, updated_at
            ) VALUES (
                %s, %s, %s, 'mixture_component', NOW(), NOW()
            )
            RETURNING id, name, formula
        """, [molecule_id, name, formula])
        
        result = cursor.fetchone()
        print(f"  Created new molecule: {name} (ID: {result['id']})")
        return result

def add_mixture_components(conn, mixture_id, components, dry_run=False):
    """Add components to a mixture."""
    for component_data in components:
        # Try to find the molecule in the database
        molecule = get_molecule_by_name(conn, component_data["name"])
        
        # If not found, create it
        if not molecule and not dry_run:
            molecule = create_molecule(conn, component_data["name"], component_data.get("formula"))
        elif not molecule and dry_run:
            print(f"  [DRY RUN] Would create molecule for {component_data['name']}")
            molecule = {"id": "new_molecule_placeholder", "name": component_data["name"]}
        
        # Add the component to the mixture
        if not dry_run:
            with conn.cursor() as cursor:
                cursor.execute("""
                    INSERT INTO mixture_components (
                        id, mixture_id, molecule_id, concentration, concentration_unit, 
                        role, created_at, updated_at, properties
                    ) VALUES (
                        %s, %s, %s, %s, %s, 
                        %s, NOW(), NOW(), %s::jsonb
                    )
                """, [
                    str(uuid.uuid4()),
                    mixture_id,
                    molecule["id"],
                    component_data["concentration"],
                    component_data["unit"],
                    "cryoprotectant",
                    json.dumps({})
                ])
                print(f"  Added component: {component_data['name']} at {component_data['concentration']} {component_data['unit']}")
        else:
            print(f"  [DRY RUN] Would add component: {component_data['name']} at {component_data['concentration']} {component_data['unit']}")

def update_mixture_description(conn, mixture_id, description, reference=None, dry_run=False):
    """Update mixture description and reference."""
    if not dry_run:
        with conn.cursor() as cursor:
            # Update the description and add reference to the properties
            properties = {"reference": reference} if reference else {}
            cursor.execute("""
                UPDATE mixtures
                SET description = %s,
                    properties = %s::jsonb,
                    updated_at = NOW()
                WHERE id = %s
            """, [description, json.dumps(properties), mixture_id])
            print(f"  Updated mixture description and properties")
    else:
        print(f"  [DRY RUN] Would update mixture description and properties")

def remove_placeholder_mixture(conn, mixture_id, dry_run=False):
    """Remove a placeholder mixture with no scientific value."""
    # First check if this mixture is referenced elsewhere
    with conn.cursor(cursor_factory=RealDictCursor) as cursor:
        cursor.execute("""
            SELECT 'experiments' as table_name, COUNT(*) as ref_count 
            FROM experiments 
            WHERE mixture_id = %s
            UNION ALL
            SELECT 'predictions' as table_name, COUNT(*) as ref_count 
            FROM predictions 
            WHERE mixture_id = %s
        """, [mixture_id, mixture_id])
        
        references = cursor.fetchall()
        
        for ref in references:
            if ref['ref_count'] > 0:
                print(f"  Cannot remove: Referenced by {ref['ref_count']} rows in {ref['table_name']}")
                return False
    
    # If no references, we can remove it
    if not dry_run:
        with conn.cursor() as cursor:
            cursor.execute("""
                DELETE FROM mixtures
                WHERE id = %s
            """, [mixture_id])
            print(f"  Removed placeholder mixture")
    else:
        print(f"  [DRY RUN] Would remove placeholder mixture")
    
    return True

def completion_needed(mixture):
    """Determine if a mixture needs completion."""
    return not mixture.get('has_components', False)

def create_audit_log(conn, action, details, dry_run=False):
    """Create an audit log entry."""
    if dry_run:
        print(f"[DRY RUN] Would create audit log entry for {action}")
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
                    %s,
                    NULL,
                    NOW(),
                    NULL,
                    %s::jsonb,
                    'database_remediation_script'
                )
            """, [action, json.dumps(details)])
        
        conn.commit()
        print(f"Created audit log entry for {action}")
    except Exception as e:
        print(f"Error creating audit log: {e}")
        print("This might be due to incompatibility with the existing audit table structure.")

def complete_mixtures(conn, dry_run=False, auto_confirm=False):
    """Complete mixture data in the database."""
    # Get known mixtures
    known_mixtures = get_known_mixtures(conn)
    if known_mixtures:
        print(f"\nFound {len(known_mixtures)} known mixtures:")
        for mixture in known_mixtures:
            status = "needs completion" if completion_needed(mixture) else "already has components"
            print(f"  {mixture['name']} ({status})")
            
        # Process each known mixture
        if not dry_run and not auto_confirm:
            confirmation = input("\nProceed with completing known mixtures? (y/n): ")
            if confirmation.lower() != 'y':
                print("Skipping known mixtures.")
                return
        
        completed_count = 0
        for mixture in known_mixtures:
            if completion_needed(mixture):
                mixture_data = KNOWN_MIXTURES.get(mixture['name'])
                if not mixture_data:
                    print(f"Warning: No data found for known mixture {mixture['name']}")
                    continue
                
                print(f"\nCompleting mixture: {mixture['name']} (ID: {mixture['id']})")
                # Add components
                add_mixture_components(conn, mixture['id'], mixture_data['components'], dry_run)
                
                # Update description
                update_mixture_description(
                    conn, 
                    mixture['id'], 
                    mixture_data['description'], 
                    mixture_data.get('reference'),
                    dry_run
                )
                
                completed_count += 1
                
        if completed_count > 0:
            action = "completed" if not dry_run else "would complete"
            print(f"\n{action.capitalize()} {completed_count} known mixtures")
            
            # Log to audit table
            if not dry_run:
                create_audit_log(conn, 'UPDATE', {
                    'action': 'complete_known_mixtures',
                    'count': completed_count,
                    'mixtures': [m['name'] for m in known_mixtures if completion_needed(m)]
                }, dry_run)
    else:
        print("No known mixtures found in the database.")
    
    # Get placeholder mixtures
    placeholder_mixtures = get_placeholder_mixtures(conn)
    if placeholder_mixtures:
        print(f"\nFound {len(placeholder_mixtures)} placeholder mixtures with no components:")
        for mixture in placeholder_mixtures:
            is_known = mixture['name'] in KNOWN_MIXTURES
            status = "will be completed" if is_known else "will be removed"
            print(f"  {mixture['name']} ({status})")
            
        # Ask for confirmation
        if not dry_run and not auto_confirm:
            confirmation = input("\nProceed with handling placeholder mixtures? (y/n): ")
            if confirmation.lower() != 'y':
                print("Skipping placeholder mixtures.")
                return
        
        removed_count = 0
        # Remove placeholder mixtures (except those we know)
        for mixture in placeholder_mixtures:
            is_known = mixture['name'] in KEEP_MIXTURES  # Use KEEP_MIXTURES, not just KNOWN_MIXTURES
            
            if not is_known:
                print(f"\nRemoving placeholder mixture: {mixture['name']} (ID: {mixture['id']})")
                if remove_placeholder_mixture(conn, mixture['id'], dry_run):
                    removed_count += 1
        
        if removed_count > 0:
            action = "removed" if not dry_run else "would remove"
            print(f"\n{action.capitalize()} {removed_count} placeholder mixtures")
            
            # Log to audit table
            if not dry_run:
                create_audit_log(conn, 'DELETE', {
                    'action': 'remove_placeholder_mixtures',
                    'count': removed_count
                }, dry_run)
    else:
        print("No placeholder mixtures found in the database.")

def main():
    parser = argparse.ArgumentParser(description="Complete and clean up mixture data")
    parser.add_argument("--dry-run", action="store_true", 
                        help="Show what would be done without making changes")
    parser.add_argument("--yes", action="store_true",
                        help="Automatically confirm changes without prompting")
    
    args = parser.parse_args()
    
    print(f"Starting mixture data completion process at {datetime.now().isoformat()}")
    if args.dry_run:
        print("DRY RUN MODE: No changes will be made to the database")
    if args.yes:
        print("AUTO-CONFIRM MODE: Changes will be applied without confirmation")
    
    # Connect to the database
    conn = connect_to_db()
    
    try:
        # Begin transaction
        conn.autocommit = False
        
        # Complete mixtures
        complete_mixtures(conn, args.dry_run, args.yes)
        
        # Commit transaction if not dry run
        if not args.dry_run:
            conn.commit()
            print("\nChanges committed successfully")
        else:
            conn.rollback()
            print("\nDry run complete, changes rolled back")
            
    except Exception as e:
        conn.rollback()
        print(f"Error during mixture completion: {e}")
        return 1
    finally:
        conn.close()
    
    print(f"Mixture data completion process finished at {datetime.now().isoformat()}")
    return 0

if __name__ == "__main__":
    sys.exit(main())