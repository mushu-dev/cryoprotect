#!/usr/bin/env python3
"""
Create indexes for common scientific query patterns in CryoProtect database.

This script identifies critical indexes needed for efficient scientific queries:
1. Indexes for molecular property filtering
2. Indexes for searching molecules by name or identifier
3. Indexes for efficient mixture component lookups
4. Indexes for experiment data access patterns
"""

import os
import sys
import json
import psycopg2
from psycopg2.extras import RealDictCursor
from dotenv import load_dotenv

# Load environment variables
load_dotenv()

# Define indexes to create
SCIENTIFIC_INDEXES = [
    {
        "name": "idx_molecules_name_trgm",
        "table": "molecules",
        "columns": "name",
        "method": "gin",
        "options": "gin_trgm_ops",
        "description": "Text search index for molecule names using trigram matching"
    },
    {
        "name": "idx_molecules_formula",
        "table": "molecules",
        "columns": "formula",
        "method": None,
        "options": None,
        "description": "Index for searching by chemical formula"
    },
    {
        "name": "idx_molecules_pubchem_cid",
        "table": "molecules", 
        "columns": "pubchem_cid",
        "method": None,
        "options": None,
        "description": "Index for looking up molecules by PubChem CID"
    },
    {
        "name": "idx_molecules_smiles_md5",
        "table": "molecules",
        "columns": "md5(smiles::text)",
        "method": None,
        "options": None,
        "description": "Index for efficient SMILES-based lookups (using MD5 hash)"
    },
    {
        "name": "idx_molecular_properties_molecule_property",
        "table": "molecular_properties",
        "columns": "molecule_id, property_type_id",
        "method": None,
        "options": None,
        "description": "Index for filtering molecules by property types"
    },
    {
        "name": "idx_molecular_properties_numeric",
        "table": "molecular_properties",
        "columns": "property_type_id, numeric_value",
        "method": None,
        "options": None,
        "description": "Index for numerical property range queries"
    },
    {
        "name": "idx_molecular_properties_text",
        "table": "molecular_properties",
        "columns": "property_type_id, text_value",
        "method": None,
        "options": None,
        "description": "Index for categorical property filtering"
    },
    {
        "name": "idx_mixture_components_mixture",
        "table": "mixture_components",
        "columns": "mixture_id",
        "method": None,
        "options": None,
        "description": "Index for looking up components of a mixture"
    },
    {
        "name": "idx_mixture_components_molecule",
        "table": "mixture_components",
        "columns": "molecule_id",
        "method": None,
        "options": None,
        "description": "Index for finding mixtures containing a specific molecule"
    },
    {
        "name": "idx_mixtures_name_trgm",
        "table": "mixtures",
        "columns": "name",
        "method": "gin",
        "options": "gin_trgm_ops",
        "description": "Text search index for mixture names using trigram matching"
    }
]

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
    return psycopg2.connect(**db_params)

def check_extensions(conn):
    """Check if required PostgreSQL extensions are installed."""
    required_extensions = ["pg_trgm"]
    missing_extensions = []
    
    with conn.cursor() as cursor:
        for ext in required_extensions:
            cursor.execute(f"SELECT 1 FROM pg_extension WHERE extname = '{ext}'")
            if not cursor.fetchone():
                missing_extensions.append(ext)
    
    return missing_extensions

def enable_extensions(conn, extensions):
    """Enable required PostgreSQL extensions."""
    with conn.cursor() as cursor:
        for ext in extensions:
            try:
                cursor.execute(f"CREATE EXTENSION IF NOT EXISTS {ext}")
                print(f"Enabled extension: {ext}")
            except Exception as e:
                print(f"Failed to enable extension {ext}: {e}")
                return False
    
    return True

def check_existing_indexes(conn):
    """Check which indexes already exist."""
    existing_indexes = set()
    
    with conn.cursor() as cursor:
        cursor.execute("""
            SELECT indexname 
            FROM pg_indexes 
            WHERE schemaname = 'public'
        """)
        for row in cursor.fetchall():
            existing_indexes.add(row[0])
    
    return existing_indexes

def create_index(conn, index_def, force_recreate=False):
    """Create an index according to definition."""
    index_name = index_def["name"]
    table = index_def["table"]
    columns = index_def["columns"]
    method = index_def.get("method")
    options = index_def.get("options")
    
    method_clause = f"USING {method}" if method else ""
    options_clause = f"({columns} {options})" if options else f"({columns})"
    
    # Build the create index statement
    create_stmt = f"CREATE INDEX {index_name} ON {table} {method_clause} {options_clause}"
    
    # Build a drop index statement if recreating
    drop_stmt = f"DROP INDEX IF EXISTS {index_name}" if force_recreate else None
    
    try:
        with conn.cursor() as cursor:
            # Drop index if recreating
            if drop_stmt:
                cursor.execute(drop_stmt)
                print(f"Dropped existing index: {index_name}")
            
            # Create index
            cursor.execute(create_stmt)
            print(f"Created index: {index_name}")
            
            return True
    except Exception as e:
        print(f"Failed to create index {index_name}: {e}")
        return False

def main():
    """Create indexes for common scientific query patterns."""
    print("Creating indexes for common scientific query patterns...")
    
    # Connect to database
    conn = connect_to_db()
    conn.autocommit = False
    
    try:
        # Check required extensions
        missing_extensions = check_extensions(conn)
        if missing_extensions:
            print(f"Enabling required extensions: {', '.join(missing_extensions)}")
            if not enable_extensions(conn, missing_extensions):
                print("Failed to enable required extensions. Aborting.")
                return
        
        # Check existing indexes
        existing_indexes = check_existing_indexes(conn)
        print(f"Found {len(existing_indexes)} existing indexes")
        
        # Create indexes
        created_indexes = []
        skipped_indexes = []
        failed_indexes = []
        
        for index_def in SCIENTIFIC_INDEXES:
            index_name = index_def["name"]
            
            # Skip if index already exists
            if index_name in existing_indexes:
                print(f"Index already exists: {index_name}")
                skipped_indexes.append(index_def)
                continue
            
            # Create the index
            if create_index(conn, index_def):
                created_indexes.append(index_def)
            else:
                failed_indexes.append(index_def)
        
        # Analyze tables for better query planning
        analyze_tables = set(index_def["table"] for index_def in SCIENTIFIC_INDEXES)
        with conn.cursor() as cursor:
            for table in analyze_tables:
                cursor.execute(f"ANALYZE {table}")
                print(f"Analyzed table: {table}")
        
        # Commit changes
        conn.commit()
        
        # Print summary
        print("\nIndex Creation Summary:")
        print(f"Created {len(created_indexes)} indexes")
        print(f"Skipped {len(skipped_indexes)} existing indexes")
        print(f"Failed to create {len(failed_indexes)} indexes")
        
        # Save results to file
        results = {
            "created_indexes": created_indexes,
            "skipped_indexes": skipped_indexes,
            "failed_indexes": failed_indexes,
            "analyzed_tables": list(analyze_tables)
        }
        with open("index_creation_report.json", "w") as f:
            json.dump(results, f, indent=2)
        
        print("\nScientific index creation completed successfully")
        print("Report saved to 'index_creation_report.json'")
    
    except Exception as e:
        conn.rollback()
        print(f"Error creating indexes: {e}")
        sys.exit(1)
    
    finally:
        conn.close()

if __name__ == "__main__":
    main()