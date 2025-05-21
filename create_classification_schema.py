#!/usr/bin/env python3
"""
Create classification schema for cryoprotectant molecules.

This script adds the necessary property types for classifying molecules by:
1. Cryoprotectant type (penetrating, non-penetrating, etc.)
2. Mechanism of action
3. Chemical class

These classifications enable better scientific organization and analysis.
"""

import os
import sys
import uuid
import psycopg2
from psycopg2.extras import RealDictCursor
from dotenv import load_dotenv

# Load environment variables
load_dotenv()

# Classification definitions
CRYOPROTECTANT_TYPES = [
    {
        "name": "cryoprotectant_type",
        "data_type": "text",
        "description": "Type of cryoprotectant based on cellular penetration ability",
        "units": None,
        "values": [
            "PENETRATING", 
            "NON_PENETRATING", 
            "OSMOLYTE", 
            "ANTIFREEZE", 
            "OTHER"
        ]
    },
    {
        "name": "mechanism_of_action",
        "data_type": "text",
        "description": "Primary mechanism by which the agent provides cryoprotection",
        "units": None,
        "values": [
            "VITRIFICATION", 
            "OSMOTIC_BUFFER", 
            "ICE_BLOCKER", 
            "MEMBRANE_STABILIZER", 
            "DEHYDRATION", 
            "ANTIOXIDANT", 
            "UNKNOWN"
        ]
    },
    {
        "name": "chemical_class",
        "data_type": "text",
        "description": "Chemical classification of the cryoprotectant",
        "units": None,
        "values": [
            "ALCOHOL", 
            "SUGAR", 
            "AMINO_ACID", 
            "AMIDE", 
            "POLYOL", 
            "SULFOXIDE", 
            "POLYMER", 
            "PEPTIDE", 
            "PROTEIN", 
            "INORGANIC", 
            "OTHER"
        ]
    },
    {
        "name": "toxicity_level",
        "data_type": "text",
        "description": "Relative toxicity level of the cryoprotectant at effective concentrations",
        "units": None,
        "values": [
            "LOW", 
            "MEDIUM", 
            "HIGH", 
            "UNKNOWN"
        ]
    }
]

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

def property_type_exists(conn, name):
    """Check if a property type already exists."""
    with conn.cursor(cursor_factory=RealDictCursor) as cursor:
        cursor.execute("""
            SELECT EXISTS(
                SELECT 1 FROM property_types WHERE name = %s
            ) AS exists
        """, (name,))
        result = cursor.fetchone()
        return result['exists']

def create_property_type(conn, classification):
    """Create a new property type for classification."""
    if property_type_exists(conn, classification['name']):
        print(f"Property type '{classification['name']}' already exists.")
        return None
        
    with conn.cursor(cursor_factory=RealDictCursor) as cursor:
        prop_type_id = str(uuid.uuid4())
        cursor.execute("""
            INSERT INTO property_types (
                id, name, data_type, description, units, created_at, updated_at
            ) VALUES (
                %s, %s, %s, %s, %s, NOW(), NOW()
            ) RETURNING id
        """, (
            prop_type_id,
            classification['name'],
            classification['data_type'],
            classification['description'],
            classification['units']
        ))
        
        result = cursor.fetchone()
        if result:
            print(f"Created property type: {classification['name']} (ID: {result['id']})")
            return result['id']
        return None

def create_audit_log(conn, details):
    """Create an audit log entry."""
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
                    'property_types',
                    'INSERT',
                    NULL,
                    NOW(),
                    NULL,
                    %s::jsonb,
                    'database_remediation_script'
                )
            """, (details,))
        
        print("Created audit log entry")
    except Exception as e:
        print(f"Error creating audit log: {e}")
        print("This might be due to incompatibility with the existing audit table structure.")

def create_enum_documentation(conn, classifications):
    """Create documentation for enum values in project properties."""
    try:
        # Check if projects table exists and has properties column
        with conn.cursor(cursor_factory=RealDictCursor) as cursor:
            cursor.execute("""
                SELECT EXISTS (
                    SELECT FROM information_schema.tables 
                    WHERE table_schema = 'public' AND table_name = 'projects'
                ) AS table_exists
            """)
            
            if not cursor.fetchone()['table_exists']:
                print("Projects table does not exist, skipping enum documentation.")
                return

            # Check if properties column exists
            cursor.execute("""
                SELECT EXISTS (
                    SELECT FROM information_schema.columns
                    WHERE table_schema = 'public' AND table_name = 'projects' AND column_name = 'properties'
                ) AS column_exists
            """)
            
            if not cursor.fetchone()['column_exists']:
                print("Properties column does not exist in projects table, skipping enum documentation.")
                return
            
            # Get current properties from first project
            cursor.execute("SELECT id, properties FROM projects LIMIT 1")
            project = cursor.fetchone()
            
            if not project:
                print("No projects found, skipping enum documentation.")
                return
            
            # Create enum documentation
            enum_docs = {}
            for classification in classifications:
                enum_docs[classification['name']] = {
                    "description": classification['description'],
                    "values": classification['values'],
                    "data_type": classification['data_type']
                }
            
            # Update properties
            properties = project['properties'] or {}
            properties['enum_documentation'] = enum_docs
            
            # Update project
            cursor.execute("""
                UPDATE projects
                SET properties = %s::jsonb
                WHERE id = %s
            """, (properties, project['id']))
            
            print(f"Added enum documentation to project {project['id']}")
    except Exception as e:
        print(f"Error creating enum documentation: {e}")

def main():
    """Create classification schema for molecules."""
    print("Creating classification schema for cryoprotectant molecules...")
    
    # Connect to the database
    conn = connect_to_db()
    conn.autocommit = False  # Start a transaction
    
    try:
        created_types = []
        
        # Create property types for classification
        for classification in CRYOPROTECTANT_TYPES:
            print(f"\nProcessing classification: {classification['name']}")
            prop_type_id = create_property_type(conn, classification)
            if prop_type_id:
                created_types.append({
                    "id": prop_type_id,
                    "name": classification['name'],
                    "values": classification['values']
                })
        
        # Create enum documentation
        print("\nCreating enum documentation...")
        create_enum_documentation(conn, CRYOPROTECTANT_TYPES)
        
        # Create audit log
        if created_types:
            create_audit_log(conn, {
                "action": "create_classification_schema",
                "created_types": created_types
            })
        
        # Commit transaction
        conn.commit()
        print("\nClassification schema created successfully")
        
    except Exception as e:
        conn.rollback()
        print(f"Error creating classification schema: {e}")
        sys.exit(1)
    
    finally:
        conn.close()

if __name__ == "__main__":
    main()