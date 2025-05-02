"""
Initial database schema migration.

Creates the initial database schema.
"""

import logging

# Configure logging
logger = logging.getLogger(__name__)

def apply(conn, environment):
    """
    Apply the migration.

    Args:
        conn: Database connection object
        environment: Target environment (development, staging, production)
    """
    logger.info(f"Applying migration 001_initial_schema in {environment} environment")

    # Example: Create tables
    conn.sql("""
        -- Create molecules table
        CREATE TABLE IF NOT EXISTS molecules (
            id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
            name TEXT NOT NULL,
            smiles TEXT NOT NULL,
            inchi TEXT,
            created_at TIMESTAMP WITH TIME ZONE DEFAULT NOW(),
            updated_at TIMESTAMP WITH TIME ZONE DEFAULT NOW()
        );

        -- Create mixtures table
        CREATE TABLE IF NOT EXISTS mixtures (
            id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
            name TEXT NOT NULL,
            description TEXT,
            created_at TIMESTAMP WITH TIME ZONE DEFAULT NOW(),
            updated_at TIMESTAMP WITH TIME ZONE DEFAULT NOW()
        );

        -- Create mixture_components table
        CREATE TABLE IF NOT EXISTS mixture_components (
            id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
            mixture_id UUID NOT NULL REFERENCES mixtures(id) ON DELETE CASCADE,
            molecule_id UUID NOT NULL REFERENCES molecules(id) ON DELETE CASCADE,
            concentration NUMERIC NOT NULL,
            units TEXT NOT NULL,
            created_at TIMESTAMP WITH TIME ZONE DEFAULT NOW(),
            updated_at TIMESTAMP WITH TIME ZONE DEFAULT NOW(),
            UNIQUE(mixture_id, molecule_id)
        );
    """).execute()

def rollback(conn, environment):
    """
    Roll back the migration.

    Args:
        conn: Database connection object
        environment: Target environment (development, staging, production)
    """
    logger.info(f"Rolling back migration 001_initial_schema in {environment} environment")

    # Example: Drop tables in reverse order
    conn.sql("""
        DROP TABLE IF EXISTS mixture_components;
        DROP TABLE IF EXISTS mixtures;
        DROP TABLE IF EXISTS molecules;
    """).execute()