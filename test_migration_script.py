#!/usr/bin/env python3
"""
CryoProtect v2 Migration Test Script

This script tests the migrate_to_plural_tables.py script in a controlled environment.
It creates a test database with singular-named tables, runs the migration,
and verifies that the migration was successful.

Usage:
    python test_migration_script.py

Requirements:
    - psycopg2
    - dotenv
"""

import os
import sys
import uuid
import logging
import subprocess
from datetime import datetime
import psycopg2
from psycopg2 import sql
from psycopg2.extras import RealDictCursor
from dotenv import load_dotenv

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler(f"migration_test_log_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

# Load environment variables
load_dotenv()

# Test database connection parameters
# Use a test database if possible, otherwise use the main database with a test schema
DB_URL = os.getenv("TEST_DATABASE_URL") or os.getenv("DATABASE_URL") or os.getenv("SUPABASE_DB_URL")
if not DB_URL:
    logger.error("Database URL not found in environment variables")
    sys.exit(1)

# Test schema name
TEST_SCHEMA = f"migration_test_{datetime.now().strftime('%Y%m%d_%H%M%S')}"

def connect_to_db():
    """
    Establish a connection to the database.
    
    Returns:
        tuple: (connection, cursor) or (None, None) if connection fails
    """
    try:
        conn = psycopg2.connect(DB_URL)
        conn.autocommit = False  # We'll manage transactions explicitly
        cursor = conn.cursor(cursor_factory=RealDictCursor)
        logger.info("Successfully connected to the database")
        return conn, cursor
    except Exception as e:
        logger.error(f"Error connecting to the database: {str(e)}")
        return None, None

def setup_test_environment(conn, cursor):
    """
    Set up a test environment with singular-named tables.
    
    Args:
        conn: Database connection
        cursor: Database cursor
        
    Returns:
        bool: True if setup was successful, False otherwise
    """
    try:
        logger.info(f"Creating test schema: {TEST_SCHEMA}")
        
        # Create test schema
        cursor.execute(
            sql.SQL("CREATE SCHEMA IF NOT EXISTS {}").format(
                sql.Identifier(TEST_SCHEMA)
            )
        )
        
        # Set search path to test schema
        cursor.execute(
            sql.SQL("SET search_path TO {}, public").format(
                sql.Identifier(TEST_SCHEMA)
            )
        )
        
        # Create singular-named tables
        logger.info("Creating singular-named tables")
        
        # Create molecule table
        cursor.execute("""
            CREATE TABLE molecule (
                id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
                pubchem_cid INTEGER UNIQUE,
                name VARCHAR(255) NOT NULL,
                formula VARCHAR(100),
                smiles TEXT,
                inchikey VARCHAR(27) UNIQUE NOT NULL,
                molecular_weight NUMERIC,
                created_at TIMESTAMPTZ DEFAULT NOW(),
                updated_at TIMESTAMPTZ DEFAULT NOW(),
                created_by UUID
            )
        """)
        
        # Create mixture table
        cursor.execute("""
            CREATE TABLE mixture (
                id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
                name VARCHAR(255) NOT NULL,
                description TEXT,
                created_at TIMESTAMPTZ DEFAULT NOW(),
                updated_at TIMESTAMPTZ DEFAULT NOW(),
                created_by UUID
            )
        """)
        
        # Create experiment table
        cursor.execute("""
            CREATE TABLE experiment (
                id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
                mixture_id UUID,
                molecule_id UUID,
                property_type_id UUID NOT NULL,
                numeric_value NUMERIC,
                text_value TEXT,
                boolean_value BOOLEAN,
                experimental_conditions TEXT,
                date_performed DATE,
                created_at TIMESTAMPTZ DEFAULT NOW(),
                updated_at TIMESTAMPTZ DEFAULT NOW(),
                created_by UUID,
                CONSTRAINT experiment_molecule_or_mixture_check CHECK (
                    (molecule_id IS NOT NULL AND mixture_id IS NULL) OR
                    (mixture_id IS NOT NULL AND molecule_id IS NULL)
                )
            )
        """)
        
        # Create prediction table
        cursor.execute("""
            CREATE TABLE prediction (
                id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
                molecule_id UUID,
                mixture_id UUID,
                property_type_id UUID NOT NULL,
                calculation_method_id UUID NOT NULL,
                numeric_value NUMERIC,
                text_value TEXT,
                boolean_value BOOLEAN,
                confidence NUMERIC,
                created_at TIMESTAMPTZ DEFAULT NOW(),
                updated_at TIMESTAMPTZ DEFAULT NOW(),
                created_by UUID,
                CONSTRAINT prediction_molecule_or_mixture_check CHECK (
                    (molecule_id IS NOT NULL AND mixture_id IS NULL) OR
                    (mixture_id IS NOT NULL AND molecule_id IS NULL)
                )
            )
        """)
        
        # Create project table
        cursor.execute("""
            CREATE TABLE project (
                id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
                name VARCHAR(255) NOT NULL,
                description TEXT,
                team_id UUID,
                created_at TIMESTAMPTZ DEFAULT NOW(),
                updated_at TIMESTAMPTZ DEFAULT NOW(),
                created_by UUID
            )
        """)
        
        # Create property_types table
        cursor.execute("""
            CREATE TABLE property_types (
                id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
                name VARCHAR(100) UNIQUE NOT NULL,
                data_type VARCHAR(50) NOT NULL,
                description TEXT,
                units VARCHAR(50),
                created_at TIMESTAMPTZ DEFAULT NOW(),
                updated_at TIMESTAMPTZ DEFAULT NOW(),
                created_by UUID
            )
        """)
        
        # Create mixture_component table
        cursor.execute("""
            CREATE TABLE mixture_component (
                id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
                mixture_id UUID NOT NULL REFERENCES mixture(id) ON DELETE CASCADE,
                molecule_id UUID NOT NULL REFERENCES molecule(id) ON DELETE CASCADE,
                concentration NUMERIC NOT NULL,
                concentration_unit VARCHAR(50) NOT NULL,
                created_at TIMESTAMPTZ DEFAULT NOW(),
                updated_at TIMESTAMPTZ DEFAULT NOW(),
                created_by UUID,
                UNIQUE(mixture_id, molecule_id)
            )
        """)
        
        # Create molecular_property table
        cursor.execute("""
            CREATE TABLE molecular_property (
                id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
                molecule_id UUID NOT NULL REFERENCES molecule(id) ON DELETE CASCADE,
                property_type_id UUID NOT NULL REFERENCES property_types(id),
                numeric_value NUMERIC,
                text_value TEXT,
                boolean_value BOOLEAN,
                created_at TIMESTAMPTZ DEFAULT NOW(),
                updated_at TIMESTAMPTZ DEFAULT NOW(),
                created_by UUID,
                CONSTRAINT molecular_property_value_check CHECK (
                    (numeric_value IS NOT NULL AND text_value IS NULL AND boolean_value IS NULL) OR
                    (numeric_value IS NULL AND text_value IS NOT NULL AND boolean_value IS NULL) OR
                    (numeric_value IS NULL AND text_value IS NULL AND boolean_value IS NOT NULL)
                )
            )
        """)
        
        # Create calculation_method table
        cursor.execute("""
            CREATE TABLE calculation_method (
                id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
                name VARCHAR(100) UNIQUE NOT NULL,
                description TEXT,
                version VARCHAR(50),
                created_at TIMESTAMPTZ DEFAULT NOW(),
                updated_at TIMESTAMPTZ DEFAULT NOW(),
                created_by UUID
            )
        """)
        
        # Create project_membership table
        cursor.execute("""
            CREATE TABLE project_membership (
                id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
                project_id UUID NOT NULL REFERENCES project(id) ON DELETE CASCADE,
                user_id UUID NOT NULL,
                role VARCHAR(50) NOT NULL,
                created_at TIMESTAMPTZ DEFAULT NOW(),
                updated_at TIMESTAMPTZ DEFAULT NOW(),
                UNIQUE(project_id, user_id)
            )
        """)
        
        # Create views
        logger.info("Creating views")
        
        # Create molecule_with_properties view
        cursor.execute("""
            CREATE VIEW molecule_with_properties AS
            SELECT 
                m.*,
                jsonb_agg(
                    jsonb_build_object(
                        'property_id', mp.id,
                        'property_type_id', mp.property_type_id,
                        'property_name', pt.name,
                        'property_data_type', pt.data_type,
                        'numeric_value', mp.numeric_value,
                        'text_value', mp.text_value,
                        'boolean_value', mp.boolean_value,
                        'units', pt.units
                    )
                ) AS properties
            FROM 
                molecule m
            LEFT JOIN 
                molecular_property mp ON m.id = mp.molecule_id
            LEFT JOIN 
                property_types pt ON mp.property_type_id = pt.id
            GROUP BY 
                m.id
        """)
        
        # Create mixture_with_components view
        cursor.execute("""
            CREATE VIEW mixture_with_components AS
            SELECT 
                mix.*,
                jsonb_agg(
                    jsonb_build_object(
                        'component_id', mc.id,
                        'molecule_id', mc.molecule_id,
                        'molecule_name', m.name,
                        'concentration', mc.concentration,
                        'concentration_unit', mc.concentration_unit
                    )
                ) AS components
            FROM 
                mixture mix
            LEFT JOIN 
                mixture_component mc ON mix.id = mc.mixture_id
            LEFT JOIN 
                molecule m ON mc.molecule_id = m.id
            GROUP BY 
                mix.id
        """)
        
        # Commit the transaction
        conn.commit()
        logger.info("Test environment setup completed successfully")
        
        return True
    except Exception as e:
        logger.error(f"Error setting up test environment: {str(e)}")
        conn.rollback()
        return False

def populate_test_data(conn, cursor):
    """
    Populate test data in the singular-named tables.
    
    Args:
        conn: Database connection
        cursor: Database cursor
        
    Returns:
        bool: True if population was successful, False otherwise
    """
    try:
        logger.info("Populating test data")
        
        # Set search path to test schema
        cursor.execute(
            sql.SQL("SET search_path TO {}, public").format(
                sql.Identifier(TEST_SCHEMA)
            )
        )
        
        # Insert property types
        property_type_ids = {}
        for prop_name in ["Vitrification Temperature", "Glass Transition Temperature", "Toxicity"]:
            prop_id = str(uuid.uuid4())
            cursor.execute("""
                INSERT INTO property_types (id, name, data_type, units)
                VALUES (%s, %s, %s, %s)
            """, (prop_id, prop_name, "numeric", "K" if "Temperature" in prop_name else ""))
            property_type_ids[prop_name] = prop_id
        
        # Insert calculation methods
        calc_method_ids = {}
        for method_name in ["Molecular Dynamics", "Machine Learning", "Quantum Chemistry"]:
            method_id = str(uuid.uuid4())
            cursor.execute("""
                INSERT INTO calculation_method (id, name, description)
                VALUES (%s, %s, %s)
            """, (method_id, method_name, f"Description for {method_name}"))
            calc_method_ids[method_name] = method_id
        
        # Insert molecules
        molecule_ids = {}
        molecules = [
            {
                "name": "Glycerol",
                "formula": "C3H8O3",
                "pubchem_cid": 753,
                "inchikey": "PEDCQBHIVMGVHV-UHFFFAOYSA-N",
                "smiles": "C(C(CO)O)O",
                "molecular_weight": 92.09
            },
            {
                "name": "DMSO",
                "formula": "C2H6OS",
                "pubchem_cid": 679,
                "inchikey": "IAZDPXIOMUYVGZ-UHFFFAOYSA-N",
                "smiles": "CS(=O)C",
                "molecular_weight": 78.13
            },
            {
                "name": "Trehalose",
                "formula": "C12H22O11",
                "pubchem_cid": 7427,
                "inchikey": "KGBXLFKZBHKPEV-UHFFFAOYSA-N",
                "smiles": "C(C1C(C(C(C(O1)OC2C(C(C(C(O2)CO)O)O)O)O)O)O)O",
                "molecular_weight": 342.30
            }
        ]
        
        for molecule in molecules:
            mol_id = str(uuid.uuid4())
            cursor.execute("""
                INSERT INTO molecule (id, name, formula, pubchem_cid, inchikey, smiles, molecular_weight)
                VALUES (%s, %s, %s, %s, %s, %s, %s)
            """, (
                mol_id,
                molecule["name"],
                molecule["formula"],
                molecule["pubchem_cid"],
                molecule["inchikey"],
                molecule["smiles"],
                molecule["molecular_weight"]
            ))
            molecule_ids[molecule["name"]] = mol_id
        
        # Insert molecular properties
        for mol_name, mol_id in molecule_ids.items():
            for prop_name, prop_id in property_type_ids.items():
                cursor.execute("""
                    INSERT INTO molecular_property (molecule_id, property_type_id, numeric_value)
                    VALUES (%s, %s, %s)
                """, (
                    mol_id,
                    prop_id,
                    round(float(hash(mol_name + prop_name) % 1000) / 10, 2)  # Generate a pseudo-random value
                ))
        
        # Insert mixtures
        mixture_ids = {}
        mixtures = [
            {"name": "CryoMix A", "description": "Standard cryoprotectant mixture"},
            {"name": "CryoMix B", "description": "Low toxicity cryoprotectant mixture"}
        ]
        
        for mixture in mixtures:
            mix_id = str(uuid.uuid4())
            cursor.execute("""
                INSERT INTO mixture (id, name, description)
                VALUES (%s, %s, %s)
            """, (
                mix_id,
                mixture["name"],
                mixture["description"]
            ))
            mixture_ids[mixture["name"]] = mix_id
        
        # Insert mixture components
        components = [
            {"mixture": "CryoMix A", "molecule": "Glycerol", "concentration": 30.0, "unit": "%w/v"},
            {"mixture": "CryoMix A", "molecule": "DMSO", "concentration": 10.0, "unit": "%w/v"},
            {"mixture": "CryoMix B", "molecule": "Trehalose", "concentration": 20.0, "unit": "%w/v"},
            {"mixture": "CryoMix B", "molecule": "Glycerol", "concentration": 15.0, "unit": "%w/v"}
        ]
        
        for comp in components:
            cursor.execute("""
                INSERT INTO mixture_component (mixture_id, molecule_id, concentration, concentration_unit)
                VALUES (%s, %s, %s, %s)
            """, (
                mixture_ids[comp["mixture"]],
                molecule_ids[comp["molecule"]],
                comp["concentration"],
                comp["unit"]
            ))
        
        # Insert predictions
        for mol_name, mol_id in molecule_ids.items():
            for prop_name, prop_id in property_type_ids.items():
                for method_name, method_id in calc_method_ids.items():
                    cursor.execute("""
                        INSERT INTO prediction (molecule_id, property_type_id, calculation_method_id, numeric_value, confidence)
                        VALUES (%s, %s, %s, %s, %s)
                    """, (
                        mol_id,
                        prop_id,
                        method_id,
                        round(float(hash(mol_name + prop_name + method_name) % 1000) / 10, 2),  # Generate a pseudo-random value
                        round(float(hash(method_name) % 100) / 100, 2)  # Generate a pseudo-random confidence
                    ))
        
        for mix_name, mix_id in mixture_ids.items():
            for prop_name, prop_id in property_type_ids.items():
                for method_name, method_id in calc_method_ids.items():
                    cursor.execute("""
                        INSERT INTO prediction (mixture_id, property_type_id, calculation_method_id, numeric_value, confidence)
                        VALUES (%s, %s, %s, %s, %s)
                    """, (
                        mix_id,
                        prop_id,
                        method_id,
                        round(float(hash(mix_name + prop_name + method_name) % 1000) / 10, 2),  # Generate a pseudo-random value
                        round(float(hash(method_name) % 100) / 100, 2)  # Generate a pseudo-random confidence
                    ))
        
        # Insert experiments
        for mol_name, mol_id in molecule_ids.items():
            for prop_name, prop_id in property_type_ids.items():
                cursor.execute("""
                    INSERT INTO experiment (molecule_id, property_type_id, numeric_value, date_performed)
                    VALUES (%s, %s, %s, %s)
                """, (
                    mol_id,
                    prop_id,
                    round(float(hash(mol_name + prop_name + "exp") % 1000) / 10, 2),  # Generate a pseudo-random value
                    datetime.now().date()
                ))
        
        for mix_name, mix_id in mixture_ids.items():
            for prop_name, prop_id in property_type_ids.items():
                cursor.execute("""
                    INSERT INTO experiment (mixture_id, property_type_id, numeric_value, date_performed)
                    VALUES (%s, %s, %s, %s)
                """, (
                    mix_id,
                    prop_id,
                    round(float(hash(mix_name + prop_name + "exp") % 1000) / 10, 2),  # Generate a pseudo-random value
                    datetime.now().date()
                ))
        
        # Insert projects
        project_ids = {}
        projects = [
            {"name": "CryoProtect Research", "description": "Main research project"}
        ]
        
        for project in projects:
            proj_id = str(uuid.uuid4())
            cursor.execute("""
                INSERT INTO project (id, name, description)
                VALUES (%s, %s, %s)
            """, (
                proj_id,
                project["name"],
                project["description"]
            ))
            project_ids[project["name"]] = proj_id
        
        # Insert project memberships
        cursor.execute("""
            INSERT INTO project_membership (project_id, user_id, role)
            VALUES (%s, %s, %s)
        """, (
            list(project_ids.values())[0],
            str(uuid.uuid4()),  # Random user ID
            "admin"
        ))
        
        # Commit the transaction
        conn.commit()
        logger.info("Test data population completed successfully")
        
        return True
    except Exception as e:
        logger.error(f"Error populating test data: {str(e)}")
        conn.rollback()
        return False

def run_migration_script():
    """
    Run the migration script in the test environment.
    
    Returns:
        bool: True if migration was successful, False otherwise
    """
    try:
        logger.info("Running migration script")
        
        # Set environment variable for the test schema
        os.environ["MIGRATION_TEST_SCHEMA"] = TEST_SCHEMA
        
        # Run the migration script with dry run first
        logger.info("Running migration script in dry run mode")
        dry_run_result = subprocess.run(
            ["python", "migrate_to_plural_tables.py", "--dry-run"],
            capture_output=True,
            text=True
        )
        
        if dry_run_result.returncode != 0:
            logger.error(f"Migration dry run failed: {dry_run_result.stderr}")
            return False
        
        logger.info("Migration dry run completed successfully")
        logger.info(dry_run_result.stdout)
        
        # Run the actual migration
        logger.info("Running actual migration")
        migration_result = subprocess.run(
            ["python", "migrate_to_plural_tables.py"],
            capture_output=True,
            text=True
        )
        
        if migration_result.returncode != 0:
            logger.error(f"Migration failed: {migration_result.stderr}")
            return False
        
        logger.info("Migration completed successfully")
        logger.info(migration_result.stdout)
        
        return True
    except Exception as e:
        logger.error(f"Error running migration script: {str(e)}")
        return False

def verify_migration(conn, cursor):
    """
    Verify that the migration was successful.
    
    Args:
        conn: Database connection
        cursor: Database cursor
        
    Returns:
        bool: True if verification was successful, False otherwise
    """
    try:
        logger.info("Verifying migration")
        
        # Set search path to test schema
        cursor.execute(
            sql.SQL("SET search_path TO {}, public").format(
                sql.Identifier(TEST_SCHEMA)
            )
        )
        
        # Check that plural tables exist
        plural_tables = ["molecules", "mixtures", "experiments", "predictions", "projects"]
        for table in plural_tables:
            cursor.execute("""
                SELECT EXISTS (
                    SELECT FROM information_schema.tables 
                    WHERE table_name = %s AND table_schema = %s
                ) as exists
            """, (table, TEST_SCHEMA))
            
            table_exists = cursor.fetchone()['exists']
            if not table_exists:
                logger.error(f"Plural table {table} does not exist")
                return False
            
            logger.info(f"Plural table {table} exists")
        
        # Check that data was copied correctly
        for singular, plural in [
            ("molecule", "molecules"),
            ("mixture", "mixtures"),
            ("experiment", "experiments"),
            ("prediction", "predictions"),
            ("project", "projects")
        ]:
            # Count rows in both tables
            cursor.execute(
                sql.SQL("SELECT COUNT(*) as count FROM {}").format(
                    sql.Identifier(singular)
                )
            )
            singular_count = cursor.fetchone()['count']
            
            cursor.execute(
                sql.SQL("SELECT COUNT(*) as count FROM {}").format(
                    sql.Identifier(plural)
                )
            )
            plural_count = cursor.fetchone()['count']
            
            if singular_count != plural_count:
                logger.error(f"Row count mismatch: {singular}={singular_count}, {plural}={plural_count}")
                return False
            
            logger.info(f"Row count match: {singular}={singular_count}, {plural}={plural_count}")
        
        # Check that foreign keys were updated
        for rel in [
            {"table": "mixture_component", "column": "mixture_id", "references": "mixtures"},
            {"table": "mixture_component", "column": "molecule_id", "references": "molecules"},
            {"table": "molecular_property", "column": "molecule_id", "references": "molecules"},
            {"table": "project_membership", "column": "project_id", "references": "projects"}
        ]:
            cursor.execute("""
                SELECT EXISTS (
                    SELECT FROM information_schema.table_constraints AS tc
                    JOIN information_schema.key_column_usage AS kcu
                        ON tc.constraint_name = kcu.constraint_name
                    JOIN information_schema.constraint_column_usage AS ccu
                        ON ccu.constraint_name = tc.constraint_name
                    WHERE tc.constraint_type = 'FOREIGN KEY'
                        AND tc.table_schema = %s
                        AND tc.table_name = %s
                        AND kcu.column_name = %s
                        AND ccu.table_name = %s
                ) as exists
            """, (TEST_SCHEMA, rel["table"], rel["column"], rel["references"]))
            
            fk_exists = cursor.fetchone()['exists']
            if not fk_exists:
                logger.error(f"Foreign key {rel['table']}.{rel['column']} -> {rel['references']} does not exist")
                return False
            
            logger.info(f"Foreign key {rel['table']}.{rel['column']} -> {rel['references']} exists")
        
        # Check that views were updated
        for view in ["molecule_with_properties", "mixture_with_components"]:
            # Get view definition
            cursor.execute("""
                SELECT pg_get_viewdef(c.oid, true) as view_definition
                FROM pg_class c
                JOIN pg_namespace n ON n.oid = c.relnamespace
                WHERE c.relkind = 'v'
                AND n.nspname = %s
                AND c.relname = %s
            """, (TEST_SCHEMA, view))
            
            view_def = cursor.fetchone()
            if not view_def:
                logger.error(f"View {view} does not exist")
                return False
            
            view_def = view_def['view_definition']
            
            # Check that view references plural tables
            for singular, plural in [
                ("molecule", "molecules"),
                ("mixture", "mixtures")
            ]:
                if f" {singular} " in view_def or f" {singular}." in view_def:
                    logger.error(f"View {view} still references singular table {singular}")
                    return False
            
            logger.info(f"View {view} references plural tables")
        
        logger.info("Migration verification completed successfully")
        return True
    except Exception as e:
        logger.error(f"Error verifying migration: {str(e)}")
        return False

def cleanup_test_environment(conn, cursor):
    """
    Clean up the test environment.
    
    Args:
        conn: Database connection
        cursor: Database cursor
        
    Returns:
        bool: True if cleanup was successful, False otherwise
    """
    try:
        logger.info(f"Cleaning up test schema: {TEST_SCHEMA}")
        
        # Drop the test schema
        cursor.execute(
            sql.SQL("DROP SCHEMA IF EXISTS {} CASCADE").format(
                sql.Identifier(TEST_SCHEMA)
            )
        )
        
        # Commit the transaction
        conn.commit()
        logger.info("Test environment cleanup completed successfully")
        
        return True
    except Exception as e:
        logger.error(f"Error cleaning up test environment: {str(e)}")
        conn.rollback()
        return False

def main():
    """Main function to run the migration test."""
    logger.info("Starting migration test")
    
    # Connect to the database
    conn, cursor = connect_to_db()
    if not conn or not cursor:
        sys.exit(1)
    
    try:
        # Setup test environment
        if not setup_test_environment(conn, cursor):
            logger.error("Failed to set up test environment")
            sys.exit(1)
        
        # Populate test data
        if not populate_test_data(conn, cursor):
            logger.error("Failed to populate test data")
            sys.exit(1)
        
        # Run migration script
        if not run_migration_script():
            logger.error("Failed to run migration script")
            sys.exit(1)
        
        # Verify migration
        if not verify_migration(conn, cursor):
            logger.error("Migration verification failed")
            sys.exit(1)
        
        logger.info("Migration test completed successfully")
    finally:
        # Clean up test environment
        if not cleanup_test_environment(conn, cursor):
            logger.warning("Failed to clean up test environment")
        
        # Close database connection
        if conn:
            conn.close()

if __name__ == "__main__":
    main()