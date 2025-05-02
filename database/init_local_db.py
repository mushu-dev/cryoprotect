#!/usr/bin/env python3
"""
Initialize local PostgreSQL database for development.

This script:
1. Creates required database and schema
2. Applies all migrations
3. Creates test data for development
"""

import os
import sys
import logging
import argparse
import psycopg2
from psycopg2.extensions import ISOLATION_LEVEL_AUTOCOMMIT
from dotenv import load_dotenv
import glob
import re

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='Initialize local PostgreSQL database')
    parser.add_argument('--reset', action='store_true', help='Reset database if it exists')
    parser.add_argument('--skip-migrations', action='store_true', help='Skip applying migrations')
    parser.add_argument('--skip-test-data', action='store_true', help='Skip creating test data')
    return parser.parse_args()

def get_connection_params():
    """Get database connection parameters from environment."""
    load_dotenv()
    
    # Get connection parameters
    params = {
        'host': os.getenv('LOCAL_DB_HOST', 'localhost'),
        'port': os.getenv('LOCAL_DB_PORT', '5432'),
        'user': os.getenv('LOCAL_DB_USER', 'postgres'),
        'password': os.getenv('LOCAL_DB_PASSWORD', '')
    }
    
    # Database name
    db_name = os.getenv('LOCAL_DB_NAME', 'cryoprotect')
    
    return params, db_name

def create_database(params, db_name, reset=False):
    """
    Create database if it doesn't exist.
    
    Args:
        params: Connection parameters
        db_name: Database name
        reset: Whether to reset database if it exists
        
    Returns:
        bool: True if database created, False if it already exists
    """
    conn = None
    try:
        # Connect to PostgreSQL server
        conn = psycopg2.connect(**params)
        conn.set_isolation_level(ISOLATION_LEVEL_AUTOCOMMIT)
        
        with conn.cursor() as cursor:
            # Check if database exists
            cursor.execute("SELECT 1 FROM pg_database WHERE datname = %s", (db_name,))
            exists = cursor.fetchone() is not None
            
            if exists:
                if reset:
                    logger.info(f"Resetting database '{db_name}'...")
                    # Terminate all connections
                    cursor.execute(f"SELECT pg_terminate_backend(pid) FROM pg_stat_activity WHERE datname = '{db_name}'")
                    # Drop database
                    cursor.execute(f"DROP DATABASE {db_name}")
                    logger.info(f"Database '{db_name}' dropped")
                else:
                    logger.info(f"Database '{db_name}' already exists")
                    return False
            
            # Create database
            cursor.execute(f"CREATE DATABASE {db_name}")
            logger.info(f"Database '{db_name}' created")
            
            return True
    except Exception as e:
        logger.error(f"Error creating database: {str(e)}")
        return False
    finally:
        if conn:
            conn.close()

def apply_migrations(params, db_name):
    """
    Apply all migrations in order.
    
    Args:
        params: Connection parameters
        db_name: Database name
        
    Returns:
        bool: True if migrations applied successfully, False otherwise
    """
    conn = None
    try:
        # Connect to database
        conn_params = params.copy()
        conn_params['dbname'] = db_name
        conn = psycopg2.connect(**conn_params)
        
        # Get migration files
        migrations_dir = os.path.join(os.getcwd(), 'migrations')
        migration_files = glob.glob(os.path.join(migrations_dir, '*.sql'))
        
        # Sort migration files by numeric prefix
        def get_migration_number(file_path):
            match = re.search(r'^(\d+)', os.path.basename(file_path))
            return int(match.group(1)) if match else 0
            
        migration_files.sort(key=get_migration_number)
        
        if not migration_files:
            logger.warning("No migration files found")
            return True
            
        # Apply migrations
        for migration_file in migration_files:
            logger.info(f"Applying migration: {os.path.basename(migration_file)}")
            
            with conn.cursor() as cursor:
                # Create migrations table if it doesn't exist
                cursor.execute("""
                    CREATE TABLE IF NOT EXISTS migrations (
                        id SERIAL PRIMARY KEY,
                        filename VARCHAR(255) NOT NULL,
                        applied_at TIMESTAMP WITH TIME ZONE DEFAULT CURRENT_TIMESTAMP
                    )
                """)
                
                # Check if migration already applied
                cursor.execute("SELECT 1 FROM migrations WHERE filename = %s", (os.path.basename(migration_file),))
                if cursor.fetchone():
                    logger.info(f"Migration {os.path.basename(migration_file)} already applied")
                    continue
                    
                # Read and apply migration
                with open(migration_file, 'r') as f:
                    sql = f.read()
                    cursor.execute(sql)
                    
                # Record migration
                cursor.execute("INSERT INTO migrations (filename) VALUES (%s)", (os.path.basename(migration_file),))
                
                # Commit transaction
                conn.commit()
                
        logger.info("All migrations applied successfully")
        return True
    except Exception as e:
        if conn:
            conn.rollback()
        logger.error(f"Error applying migrations: {str(e)}")
        return False
    finally:
        if conn:
            conn.close()

def create_test_data(params, db_name):
    """
    Create test data for development.
    
    Args:
        params: Connection parameters
        db_name: Database name
        
    Returns:
        bool: True if test data created successfully, False otherwise
    """
    conn = None
    try:
        # Connect to database
        conn_params = params.copy()
        conn_params['dbname'] = db_name
        conn = psycopg2.connect(**conn_params)
        
        logger.info("Creating test data...")
        
        with conn.cursor() as cursor:
            # Create test user
            cursor.execute("""
                INSERT INTO auth.users (id, email, encrypted_password, email_confirmed_at, created_at, updated_at)
                VALUES (
                    '00000000-0000-0000-0000-000000000000',
                    'test@example.com',
                    '$2a$10$abcdefghijklmnopqrstuvwxyz',
                    NOW(),
                    NOW(),
                    NOW()
                )
                ON CONFLICT (id) DO NOTHING
            """)
            
            # Create test user profile
            cursor.execute("""
                INSERT INTO user_profile (id, auth_user_id, display_name, email, affiliation, created_at, updated_at)
                VALUES (
                    '00000000-0000-0000-0000-000000000001',
                    '00000000-0000-0000-0000-000000000000',
                    'Test User',
                    'test@example.com',
                    'Test Organization',
                    NOW(),
                    NOW()
                )
                ON CONFLICT (id) DO NOTHING
            """)
            
            # Create test molecules
            cursor.execute("""
                INSERT INTO molecules (
                    id, name, formula, molecular_weight, smiles, inchi, inchi_key,
                    chembl_id, pubchem_cid, data_source, created_at, updated_at
                )
                VALUES
                    (
                        '00000000-0000-0000-0000-000000000010',
                        'Dimethyl sulfoxide',
                        'C2H6OS',
                        78.13,
                        'CS(=O)C',
                        'InChI=1S/C2H6OS/c1-4(2)3/h1-2H3',
                        'IAZDPXIOMUYVGZ-UHFFFAOYSA-N',
                        'CHEMBL422',
                        '679',
                        'test_data',
                        NOW(),
                        NOW()
                    ),
                    (
                        '00000000-0000-0000-0000-000000000011',
                        'Glycerol',
                        'C3H8O3',
                        92.09,
                        'C(C(CO)O)O',
                        'InChI=1S/C3H8O3/c4-1-3(6)2-5/h3-6H,1-2H2',
                        'PEDCQBHIVMGVHV-UHFFFAOYSA-N',
                        'CHEMBL692',
                        '753',
                        'test_data',
                        NOW(),
                        NOW()
                    ),
                    (
                        '00000000-0000-0000-0000-000000000012',
                        'Trehalose',
                        'C12H22O11',
                        342.30,
                        'C(C1C(C(C(C(O1)OC2C(C(C(C(O2)CO)O)O)O)O)O)O)O',
                        'InChI=1S/C12H22O11/c13-1-4-7(16)8(17)9(18)11(21-4)23-12-10(19)6(15)5(14)3(2-20)22-12/h3-19H,1-2H2/t3-,4+,5-,6-,7+,8+,9-,10-,11-,12+/m1/s1',
                        'HBVRQNNGHBPKND-DZGCQCFKSA-N',
                        'CHEMBL15132',
                        '7427',
                        'test_data',
                        NOW(),
                        NOW()
                    )
                ON CONFLICT (id) DO NOTHING
            """)
            
            # Create test property types
            cursor.execute("""
                INSERT INTO property_types (id, name, description, data_type, unit, created_at, updated_at)
                VALUES
                    (
                        '00000000-0000-0000-0000-000000000020',
                        'logP',
                        'Octanol-water partition coefficient',
                        'numeric',
                        '',
                        NOW(),
                        NOW()
                    ),
                    (
                        '00000000-0000-0000-0000-000000000021',
                        'molecular_weight',
                        'Molecular weight',
                        'numeric',
                        'g/mol',
                        NOW(),
                        NOW()
                    ),
                    (
                        '00000000-0000-0000-0000-000000000022',
                        'h_bond_donors',
                        'Number of hydrogen bond donors',
                        'numeric',
                        '',
                        NOW(),
                        NOW()
                    ),
                    (
                        '00000000-0000-0000-0000-000000000023',
                        'h_bond_acceptors',
                        'Number of hydrogen bond acceptors',
                        'numeric',
                        '',
                        NOW(),
                        NOW()
                    )
                ON CONFLICT (id) DO NOTHING
            """)
            
            # Create test molecular properties
            cursor.execute("""
                INSERT INTO molecular_properties (
                    id, molecule_id, property_type_id, value, source, confidence, created_at, updated_at
                )
                VALUES
                    (
                        '00000000-0000-0000-0000-000000000030',
                        '00000000-0000-0000-0000-000000000010',
                        '00000000-0000-0000-0000-000000000020',
                        -1.35,
                        'test_data',
                        0.95,
                        NOW(),
                        NOW()
                    ),
                    (
                        '00000000-0000-0000-0000-000000000031',
                        '00000000-0000-0000-0000-000000000010',
                        '00000000-0000-0000-0000-000000000021',
                        78.13,
                        'test_data',
                        1.0,
                        NOW(),
                        NOW()
                    ),
                    (
                        '00000000-0000-0000-0000-000000000032',
                        '00000000-0000-0000-0000-000000000010',
                        '00000000-0000-0000-0000-000000000022',
                        0,
                        'test_data',
                        1.0,
                        NOW(),
                        NOW()
                    ),
                    (
                        '00000000-0000-0000-0000-000000000033',
                        '00000000-0000-0000-0000-000000000010',
                        '00000000-0000-0000-0000-000000000023',
                        1,
                        'test_data',
                        1.0,
                        NOW(),
                        NOW()
                    ),
                    (
                        '00000000-0000-0000-0000-000000000034',
                        '00000000-0000-0000-0000-000000000011',
                        '00000000-0000-0000-0000-000000000020',
                        -1.76,
                        'test_data',
                        0.9,
                        NOW(),
                        NOW()
                    ),
                    (
                        '00000000-0000-0000-0000-000000000035',
                        '00000000-0000-0000-0000-000000000011',
                        '00000000-0000-0000-0000-000000000021',
                        92.09,
                        'test_data',
                        1.0,
                        NOW(),
                        NOW()
                    ),
                    (
                        '00000000-0000-0000-0000-000000000036',
                        '00000000-0000-0000-0000-000000000012',
                        '00000000-0000-0000-0000-000000000020',
                        -3.77,
                        'test_data',
                        0.85,
                        NOW(),
                        NOW()
                    ),
                    (
                        '00000000-0000-0000-0000-000000000037',
                        '00000000-0000-0000-0000-000000000012',
                        '00000000-0000-0000-0000-000000000021',
                        342.30,
                        'test_data',
                        1.0,
                        NOW(),
                        NOW()
                    )
                ON CONFLICT (id) DO NOTHING
            """)
            
            # Create test mixtures
            cursor.execute("""
                INSERT INTO mixtures (id, name, description, creator_id, created_at, updated_at)
                VALUES
                    (
                        '00000000-0000-0000-0000-000000000040',
                        'DMSO-Glycerol Mix',
                        'A test mixture of DMSO and glycerol',
                        '00000000-0000-0000-0000-000000000001',
                        NOW(),
                        NOW()
                    ),
                    (
                        '00000000-0000-0000-0000-000000000041',
                        'Trehalose Solution',
                        'A test mixture with trehalose',
                        '00000000-0000-0000-0000-000000000001',
                        NOW(),
                        NOW()
                    )
                ON CONFLICT (id) DO NOTHING
            """)
            
            # Create test mixture components
            cursor.execute("""
                INSERT INTO mixture_components (
                    id, mixture_id, molecule_id, concentration, concentration_unit, created_at, updated_at
                )
                VALUES
                    (
                        '00000000-0000-0000-0000-000000000050',
                        '00000000-0000-0000-0000-000000000040',
                        '00000000-0000-0000-0000-000000000010',
                        10.0,
                        'percent_v_v',
                        NOW(),
                        NOW()
                    ),
                    (
                        '00000000-0000-0000-0000-000000000051',
                        '00000000-0000-0000-0000-000000000040',
                        '00000000-0000-0000-0000-000000000011',
                        5.0,
                        'percent_v_v',
                        NOW(),
                        NOW()
                    ),
                    (
                        '00000000-0000-0000-0000-000000000052',
                        '00000000-0000-0000-0000-000000000041',
                        '00000000-0000-0000-0000-000000000012',
                        0.3,
                        'molar',
                        NOW(),
                        NOW()
                    )
                ON CONFLICT (id) DO NOTHING
            """)
            
            # Commit transaction
            conn.commit()
            
        logger.info("Test data created successfully")
        return True
    except Exception as e:
        if conn:
            conn.rollback()
        logger.error(f"Error creating test data: {str(e)}")
        return False
    finally:
        if conn:
            conn.close()

def main():
    """Main function."""
    args = parse_args()
    
    # Get connection parameters
    params, db_name = get_connection_params()
    
    # Create database
    created = create_database(params, db_name, args.reset)
    
    # Apply migrations
    if not args.skip_migrations:
        if not apply_migrations(params, db_name):
            logger.error("Failed to apply migrations")
            return 1
            
    # Create test data
    if not args.skip_test_data:
        if not create_test_data(params, db_name):
            logger.error("Failed to create test data")
            return 1
            
    logger.info("Database initialization completed successfully")
    return 0

if __name__ == '__main__':
    sys.exit(main())