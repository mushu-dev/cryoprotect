#!/usr/bin/env python3
"""
Script to check molecule and property counts in the database.
"""

import os
import sys
import psycopg2
from psycopg2.extras import RealDictCursor

# Database connection parameters
DB_PARAMS = {
    'host': 'aws-0-us-east-1.pooler.supabase.com',
    'port': 5432,
    'dbname': 'postgres',
    'user': 'postgres.tsdlmynydfuypiugmkev',
    'password': 'LDHt$rkaM&Gmf3X@LQ37',
    'sslmode': 'require'
}

def get_db_connection():
    """Get database connection using psycopg2."""
    return psycopg2.connect(**DB_PARAMS)

def get_molecule_count():
    """Get the total count of molecules in the database."""
    conn = get_db_connection()
    cursor = conn.cursor()
    
    cursor.execute("SELECT COUNT(*) FROM molecules")
    count = cursor.fetchone()[0]
    
    print(f"Total molecules: {count}")
    
    cursor.close()
    conn.close()
    return count

def get_property_count():
    """Get the total count of molecular properties in the database."""
    conn = get_db_connection()
    cursor = conn.cursor()
    
    cursor.execute("SELECT COUNT(*) FROM molecular_properties")
    count = cursor.fetchone()[0]
    
    print(f"Total properties: {count}")
    
    cursor.close()
    conn.close()
    return count

def get_counts_by_source():
    """Get molecule counts by data source."""
    conn = get_db_connection()
    cursor = conn.cursor(cursor_factory=RealDictCursor)
    
    cursor.execute("""
        SELECT data_source, COUNT(*) as count
        FROM molecules
        WHERE data_source IS NOT NULL
        GROUP BY data_source
        ORDER BY count DESC
    """)
    sources = cursor.fetchall()
    
    print("\nMolecules by source:")
    for src in sources:
        print(f"  {src['data_source']}: {src['count']}")
    
    cursor.close()
    conn.close()

def get_pubchem_molecules():
    """Get the count of molecules with PubChem CIDs."""
    conn = get_db_connection()
    cursor = conn.cursor()
    
    cursor.execute("SELECT COUNT(*) FROM molecules WHERE pubchem_cid IS NOT NULL")
    count = cursor.fetchone()[0]
    
    print(f"\nMolecules with PubChem CID: {count}")
    
    cursor.close()
    conn.close()
    return count

def get_chembl_molecules():
    """Get the count of molecules with ChEMBL IDs."""
    conn = get_db_connection()
    cursor = conn.cursor()
    
    cursor.execute("SELECT COUNT(*) FROM molecules WHERE chembl_id IS NOT NULL")
    count = cursor.fetchone()[0]
    
    print(f"Molecules with ChEMBL ID: {count}")
    
    cursor.close()
    conn.close()
    return count

if __name__ == "__main__":
    print("Checking database counts...")
    get_molecule_count()
    get_property_count()
    get_counts_by_source()
    get_pubchem_molecules()
    get_chembl_molecules()
    print("\nDone.")