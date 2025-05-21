#!/usr/bin/env python3
"""
Create experimental linkage schema between molecules and results.

This script implements a schema for linking molecules and mixtures to
experimental results, allowing for better tracking of:
1. Cryopreservation outcomes for specific molecules/mixtures
2. Experimental conditions used with each molecule/mixture
3. Tissue or cell types tested with each molecule/mixture
"""

import os
import sys
import uuid
import json
import psycopg2
from psycopg2.extras import RealDictCursor
from dotenv import load_dotenv

# Load environment variables
load_dotenv()

# Table definitions for experimental linkage schema
EXPERIMENTAL_TABLES = [
    {
        "name": "experiment_types",
        "sql": """
            CREATE TABLE IF NOT EXISTS experiment_types (
                id UUID PRIMARY KEY,
                name VARCHAR(255) NOT NULL,
                description TEXT,
                protocol_details JSONB,
                created_by UUID,
                created_at TIMESTAMPTZ NOT NULL DEFAULT NOW(),
                updated_at TIMESTAMPTZ NOT NULL DEFAULT NOW()
            )
        """,
        "indexes": [
            "CREATE INDEX IF NOT EXISTS idx_experiment_types_name ON experiment_types (name)"
        ],
        "sample_data": [
            {
                "id": "e8c1a982-9e4d-4b8a-a71d-192a0c29d4a2",
                "name": "Slow Freezing",
                "description": "Controlled rate freezing protocol with cryoprotectants"
            },
            {
                "id": "f3b25c7d-6e8a-4c9b-b013-8a3d7c4e6df5",
                "name": "Vitrification",
                "description": "Ultra-rapid cooling to achieve glass transition"
            },
            {
                "id": "a1c34e5f-7b8d-4c9e-a012-3f4e5d6c7b8a",
                "name": "Cryoprotectant Toxicity",
                "description": "Assessment of cell survival after cryoprotectant exposure"
            }
        ]
    },
    {
        "name": "tissue_types",
        "sql": """
            CREATE TABLE IF NOT EXISTS tissue_types (
                id UUID PRIMARY KEY,
                name VARCHAR(255) NOT NULL,
                description TEXT,
                species VARCHAR(255),
                taxonomy_id INTEGER,
                created_by UUID,
                created_at TIMESTAMPTZ NOT NULL DEFAULT NOW(),
                updated_at TIMESTAMPTZ NOT NULL DEFAULT NOW()
            )
        """,
        "indexes": [
            "CREATE INDEX IF NOT EXISTS idx_tissue_types_name ON tissue_types (name)",
            "CREATE INDEX IF NOT EXISTS idx_tissue_types_species ON tissue_types (species)"
        ],
        "sample_data": [
            {
                "id": "b2c3d4e5-f6g7-4h8i-9j0k-0123456789ab",
                "name": "Human Oocyte",
                "description": "Human female reproductive cells",
                "species": "Homo sapiens",
                "taxonomy_id": 9606
            },
            {
                "id": "c3d4e5f6-g7h8-4i9j-0k1l-0123456789ab",
                "name": "Mouse Embryo",
                "description": "Preimplantation mouse embryo",
                "species": "Mus musculus",
                "taxonomy_id": 10090
            },
            {
                "id": "d4e5f6g7-h8i9-4j0k-1l2m-0123456789ab",
                "name": "Bovine Sperm",
                "description": "Bull sperm cells",
                "species": "Bos taurus",
                "taxonomy_id": 9913
            }
        ]
    },
    {
        "name": "experiments",
        "sql": """
            CREATE TABLE IF NOT EXISTS experiments (
                id UUID PRIMARY KEY,
                experiment_type_id UUID NOT NULL REFERENCES experiment_types(id),
                title VARCHAR(255) NOT NULL,
                description TEXT,
                date_performed DATE,
                temperature NUMERIC,
                cooling_rate NUMERIC,
                thawing_rate NUMERIC,
                additional_parameters JSONB,
                created_by UUID,
                created_at TIMESTAMPTZ NOT NULL DEFAULT NOW(),
                updated_at TIMESTAMPTZ NOT NULL DEFAULT NOW()
            )
        """,
        "indexes": [
            "CREATE INDEX IF NOT EXISTS idx_experiments_experiment_type ON experiments (experiment_type_id)",
            "CREATE INDEX IF NOT EXISTS idx_experiments_date ON experiments (date_performed)"
        ],
        "sample_data": [
            {
                "id": "e5f6g7h8-i9j0-4k1l-2m3n-0123456789ab",
                "experiment_type_id": "e8c1a982-9e4d-4b8a-a71d-192a0c29d4a2",
                "title": "Standard Slow Freezing Protocol",
                "description": "Controlled rate freezing with DMSO",
                "date_performed": "2025-01-15",
                "temperature": -80.0,
                "cooling_rate": -1.0,
                "thawing_rate": 5.0,
                "additional_parameters": {}
            },
            {
                "id": "f6g7h8i9-j0k1-4l2m-3n4o-0123456789ab",
                "experiment_type_id": "f3b25c7d-6e8a-4c9b-b013-8a3d7c4e6df5",
                "title": "Oocyte Vitrification Protocol",
                "description": "Ultra-rapid cooling using ethylene glycol and DMSO mixture",
                "date_performed": "2025-02-20",
                "temperature": -196.0,
                "cooling_rate": -20000.0,
                "thawing_rate": 2000.0,
                "additional_parameters": {"straws_used": True, "liquid_nitrogen_submersion": True}
            },
            {
                "id": "g7h8i9j0-k1l2-4m3n-4o5p-0123456789ab",
                "experiment_type_id": "a1c34e5f-7b8d-4c9e-a012-3f4e5d6c7b8a",
                "title": "Toxicity Assessment of Novel Cryoprotectants",
                "description": "Evaluation of cell viability after exposure to new cryoprotectant compounds",
                "date_performed": "2025-03-10",
                "temperature": 22.0,
                "cooling_rate": None,
                "thawing_rate": None,
                "additional_parameters": {"exposure_time_minutes": 30, "concentrations_tested": [5, 10, 15, 20]}
            }
        ]
    },
    {
        "name": "experiment_results",
        "sql": """
            CREATE TABLE IF NOT EXISTS experiment_results (
                id UUID PRIMARY KEY,
                experiment_id UUID NOT NULL REFERENCES experiments(id),
                molecule_id UUID REFERENCES molecules(id),
                mixture_id UUID REFERENCES mixtures(id),
                tissue_type_id UUID NOT NULL REFERENCES tissue_types(id),
                concentration NUMERIC,
                concentration_unit VARCHAR(50),
                viability_percentage NUMERIC,
                recovery_rate NUMERIC,
                functionality_score NUMERIC,
                result_details JSONB,
                notes TEXT,
                created_by UUID,
                created_at TIMESTAMPTZ NOT NULL DEFAULT NOW(),
                updated_at TIMESTAMPTZ NOT NULL DEFAULT NOW(),
                CONSTRAINT check_molecule_or_mixture CHECK (
                    (molecule_id IS NOT NULL AND mixture_id IS NULL) OR
                    (molecule_id IS NULL AND mixture_id IS NOT NULL)
                )
            )
        """,
        "indexes": [
            "CREATE INDEX IF NOT EXISTS idx_experiment_results_experiment ON experiment_results (experiment_id)",
            "CREATE INDEX IF NOT EXISTS idx_experiment_results_molecule ON experiment_results (molecule_id)",
            "CREATE INDEX IF NOT EXISTS idx_experiment_results_mixture ON experiment_results (mixture_id)",
            "CREATE INDEX IF NOT EXISTS idx_experiment_results_tissue ON experiment_results (tissue_type_id)",
            "CREATE INDEX IF NOT EXISTS idx_experiment_results_viability ON experiment_results (viability_percentage)"
        ]
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

def check_existing_tables(conn):
    """Check which tables already exist."""
    existing_tables = set()

    with conn.cursor() as cursor:
        cursor.execute("""
            SELECT table_name
            FROM information_schema.tables
            WHERE table_schema = 'public'
        """)
        for row in cursor.fetchall():
            existing_tables.add(row[0])

    return existing_tables

def get_sample_molecule_ids(conn, limit=3):
    """Get sample molecule IDs for experiment results."""
    molecule_ids = []

    try:
        with conn.cursor(cursor_factory=RealDictCursor) as cursor:
            # Try to get molecules that are known cryoprotectants first
            cursor.execute("""
                SELECT id, name
                FROM molecules
                WHERE is_cryoprotectant = true
                LIMIT %s
            """, (limit,))

            rows = cursor.fetchall()

            # If we don't have enough cryoprotectants, get any molecules
            if len(rows) < limit:
                cursor.execute("""
                    SELECT id, name
                    FROM molecules
                    LIMIT %s
                """, (limit - len(rows),))
                rows.extend(cursor.fetchall())

            molecule_ids = [{'id': row['id'], 'name': row['name']} for row in rows]
            print(f"Found {len(molecule_ids)} sample molecules for experiment results")

    except Exception as e:
        print(f"Error getting sample molecule IDs: {e}")
        # Generate mock UUIDs if we can't get real ones
        molecule_ids = [
            {'id': str(uuid.uuid4()), 'name': "DMSO"},
            {'id': str(uuid.uuid4()), 'name': "Glycerol"},
            {'id': str(uuid.uuid4()), 'name': "Ethylene Glycol"}
        ]
        print(f"Using mock molecule IDs for sample data")

    return molecule_ids

def create_table(conn, table_def):
    """Create a table and its indexes."""
    table_name = table_def["name"]
    table_sql = table_def["sql"]
    indexes = table_def.get("indexes", [])
    
    with conn.cursor() as cursor:
        try:
            # Create table
            cursor.execute(table_sql)
            print(f"Created table: {table_name}")
            
            # Create indexes
            for index_sql in indexes:
                cursor.execute(index_sql)
            
            return True
        except Exception as e:
            print(f"Error creating table {table_name}: {e}")
            return False

def insert_sample_data(conn, table_def):
    """Insert sample data into a table."""
    table_name = table_def["name"]
    sample_data = table_def.get("sample_data", [])

    if not sample_data:
        return True

    try:
        with conn.cursor() as cursor:
            # Check if table already has data
            cursor.execute(f"SELECT COUNT(*) FROM {table_name}")
            count = cursor.fetchone()[0]

            if count > 0:
                print(f"Table {table_name} already has {count} rows, skipping sample data")
                return True

            # Insert sample data
            for item in sample_data:
                item_copy = item.copy()

                # Validate UUID format if present, or generate a new one
                if "id" in item_copy:
                    try:
                        # Try to parse the UUID to validate it
                        uuid.UUID(item_copy["id"])
                    except ValueError:
                        # If invalid, generate a new UUID
                        print(f"Invalid UUID format in {table_name} sample data, generating new UUID")
                        item_copy["id"] = str(uuid.uuid4())

                # Generate field names and placeholders
                fields = ", ".join(item_copy.keys())
                placeholders = ", ".join(["%s"] * len(item_copy))

                # Build and execute INSERT statement
                sql = f"INSERT INTO {table_name} ({fields}) VALUES ({placeholders})"
                cursor.execute(sql, list(item_copy.values()))

            print(f"Inserted {len(sample_data)} sample rows into {table_name}")
            return True

    except Exception as e:
        print(f"Error inserting sample data into {table_name}: {e}")
        return False

def create_views(conn):
    """Create views for common experiment result queries."""
    views = [
        {
            "name": "molecule_experiment_results",
            "sql": """
                CREATE OR REPLACE VIEW molecule_experiment_results AS
                SELECT 
                    er.id,
                    m.name as molecule_name,
                    m.smiles as molecule_smiles,
                    et.name as experiment_type,
                    tt.name as tissue_type,
                    tt.species,
                    e.temperature,
                    e.cooling_rate,
                    e.thawing_rate,
                    er.concentration,
                    er.concentration_unit,
                    er.viability_percentage,
                    er.recovery_rate,
                    er.functionality_score,
                    er.result_details,
                    e.date_performed
                FROM 
                    experiment_results er
                JOIN 
                    experiments e ON er.experiment_id = e.id
                JOIN 
                    experiment_types et ON e.experiment_type_id = et.id
                JOIN 
                    tissue_types tt ON er.tissue_type_id = tt.id
                JOIN 
                    molecules m ON er.molecule_id = m.id
                WHERE 
                    er.molecule_id IS NOT NULL
            """
        },
        {
            "name": "mixture_experiment_results",
            "sql": """
                CREATE OR REPLACE VIEW mixture_experiment_results AS
                SELECT 
                    er.id,
                    mx.name as mixture_name,
                    et.name as experiment_type,
                    tt.name as tissue_type,
                    tt.species,
                    e.temperature,
                    e.cooling_rate,
                    e.thawing_rate,
                    er.concentration,
                    er.concentration_unit,
                    er.viability_percentage,
                    er.recovery_rate,
                    er.functionality_score,
                    er.result_details,
                    e.date_performed
                FROM 
                    experiment_results er
                JOIN 
                    experiments e ON er.experiment_id = e.id
                JOIN 
                    experiment_types et ON e.experiment_type_id = et.id
                JOIN 
                    tissue_types tt ON er.tissue_type_id = tt.id
                JOIN 
                    mixtures mx ON er.mixture_id = mx.id
                WHERE 
                    er.mixture_id IS NOT NULL
            """
        },
        {
            "name": "experiment_summary",
            "sql": """
                CREATE OR REPLACE VIEW experiment_summary AS
                SELECT 
                    et.name as experiment_type,
                    tt.name as tissue_type,
                    CASE 
                        WHEN er.molecule_id IS NOT NULL THEN m.name
                        WHEN er.mixture_id IS NOT NULL THEN mx.name
                        ELSE NULL
                    END as substance_name,
                    CASE 
                        WHEN er.molecule_id IS NOT NULL THEN 'molecule'
                        WHEN er.mixture_id IS NOT NULL THEN 'mixture'
                        ELSE NULL
                    END as substance_type,
                    AVG(er.viability_percentage) as avg_viability,
                    COUNT(*) as experiment_count,
                    MIN(e.date_performed) as first_experiment,
                    MAX(e.date_performed) as last_experiment
                FROM 
                    experiment_results er
                JOIN 
                    experiments e ON er.experiment_id = e.id
                JOIN 
                    experiment_types et ON e.experiment_type_id = et.id
                JOIN 
                    tissue_types tt ON er.tissue_type_id = tt.id
                LEFT JOIN 
                    molecules m ON er.molecule_id = m.id
                LEFT JOIN 
                    mixtures mx ON er.mixture_id = mx.id
                GROUP BY 
                    et.name, tt.name, substance_name, substance_type
            """
        }
    ]
    
    for view in views:
        try:
            with conn.cursor() as cursor:
                cursor.execute(view["sql"])
                print(f"Created view: {view['name']}")
        except Exception as e:
            print(f"Error creating view {view['name']}: {e}")

def create_experiment_results_sample_data(conn):
    """Create sample data for experiment results."""
    print("Creating sample experiment results data...")

    # Get sample molecule IDs
    molecule_samples = get_sample_molecule_ids(conn, limit=3)

    # Create sample experiment results data
    experiment_results_data = []

    # Check if experiment_results table exists
    with conn.cursor() as cursor:
        cursor.execute("SELECT EXISTS (SELECT 1 FROM information_schema.tables WHERE table_name = 'experiment_results')")
        table_exists = cursor.fetchone()[0]

        if not table_exists:
            print("Experiment results table does not exist yet, skipping sample data")
            return

        # Check if table already has data
        cursor.execute("SELECT COUNT(*) FROM experiment_results")
        count = cursor.fetchone()[0]

        if count > 0:
            print(f"Table experiment_results already has {count} rows, skipping sample data")
            return

    # Get experiment IDs
    with conn.cursor() as cursor:
        cursor.execute("SELECT id FROM experiments")
        experiment_ids = [row[0] for row in cursor.fetchall()]

        if not experiment_ids:
            print("No experiments found, skipping experiment results sample data")
            return

    # Get tissue type IDs
    with conn.cursor() as cursor:
        cursor.execute("SELECT id FROM tissue_types")
        tissue_type_ids = [row[0] for row in cursor.fetchall()]

        if not tissue_type_ids:
            print("No tissue types found, skipping experiment results sample data")
            return

    # Create sample data for each molecule
    for i, molecule in enumerate(molecule_samples):
        # Create results for each experiment
        for j, experiment_id in enumerate(experiment_ids):
            # Use different tissue types if available
            tissue_type_id = tissue_type_ids[j % len(tissue_type_ids)]

            # Create result with varying parameters
            result = {
                "id": str(uuid.uuid4()),
                "experiment_id": experiment_id,
                "molecule_id": molecule['id'],
                "tissue_type_id": tissue_type_id,
                "concentration": 5 + (i * 5),  # 5%, 10%, 15%
                "concentration_unit": "%w/v",
                "viability_percentage": 70 + (i * 10) - (j * 5),  # Varying viability
                "recovery_rate": 65 + (i * 8) - (j * 4),
                "functionality_score": 3 + (i * 0.5),
                "notes": f"Sample experiment with {molecule['name']}"
            }

            experiment_results_data.append(result)

    # Insert sample data
    with conn.cursor() as cursor:
        for result in experiment_results_data:
            fields = ", ".join(result.keys())
            placeholders = ", ".join(["%s"] * len(result))

            sql = f"INSERT INTO experiment_results ({fields}) VALUES ({placeholders})"
            cursor.execute(sql, list(result.values()))

    print(f"Created {len(experiment_results_data)} sample experiment results")

def main():
    """Create experimental linkage schema between molecules and results."""
    print("Creating experimental linkage schema...")

    # Connect to database
    conn = connect_to_db()
    conn.autocommit = False

    try:
        # Check existing tables
        existing_tables = check_existing_tables(conn)

        # Create tables
        tables_created = []
        tables_skipped = []

        for table_def in EXPERIMENTAL_TABLES:
            table_name = table_def["name"]

            if table_name in existing_tables:
                print(f"Table already exists: {table_name}")
                tables_skipped.append(table_name)
                continue

            if create_table(conn, table_def):
                tables_created.append(table_name)

                # Insert sample data for reference tables
                if "sample_data" in table_def:
                    insert_sample_data(conn, table_def)

        # Create experiment results sample data after all tables are created
        if "experiment_results" in tables_created and "experiments" in tables_created:
            create_experiment_results_sample_data(conn)

        # Create views for convenient querying
        if tables_created:
            create_views(conn)

        # Commit changes
        conn.commit()

        # Print summary
        print("\nExperimental Linkage Schema Creation Summary:")
        print(f"Created {len(tables_created)} tables: {', '.join(tables_created)}")
        print(f"Skipped {len(tables_skipped)} existing tables: {', '.join(tables_skipped)}")

        # Save results to file
        results = {
            "tables_created": tables_created,
            "tables_skipped": tables_skipped,
            "schema_definition": EXPERIMENTAL_TABLES
        }
        with open("experimental_linkage_schema_report.json", "w") as f:
            json.dump(results, f, indent=2)

        print("\nExperimental linkage schema creation completed successfully")
        print("Report saved to 'experimental_linkage_schema_report.json'")

    except Exception as e:
        conn.rollback()
        print(f"Error creating experimental linkage schema: {e}")
        sys.exit(1)

    finally:
        conn.close()

if __name__ == "__main__":
    main()