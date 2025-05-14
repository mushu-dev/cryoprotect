#!/usr/bin/env python3
"""
Create basic database tables for Heroku deployment.
This script creates the minimum required tables without using Supabase-specific features.
"""

import os
import sys
import logging
from urllib.parse import urlparse
import psycopg2

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[logging.StreamHandler(sys.stdout)]
)
logger = logging.getLogger('create_heroku_tables')

def create_tables(connection_params):
    """Create the basic tables required for the CryoProtect application."""
    try:
        # Connect to the database
        conn = psycopg2.connect(**connection_params)
        cursor = conn.cursor()
        
        # Set autocommit to True to execute each statement immediately
        conn.autocommit = True
        
        # Enable the pgcrypto extension for UUID generation
        logger.info("Enabling pgcrypto extension...")
        cursor.execute("CREATE EXTENSION IF NOT EXISTS pgcrypto;")
        
        # Create molecules table
        logger.info("Creating molecules table...")
        cursor.execute("""
        CREATE TABLE IF NOT EXISTS molecules (
            id SERIAL PRIMARY KEY,
            name VARCHAR(255),
            smiles TEXT,
            pubchem_cid VARCHAR(50),
            molecular_formula VARCHAR(100),
            molecular_weight FLOAT,
            created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
            updated_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
        );
        """)
        
        # Create property_types table
        logger.info("Creating property_types table...")
        cursor.execute("""
        CREATE TABLE IF NOT EXISTS property_types (
            id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
            name VARCHAR(255) NOT NULL,
            description TEXT,
            unit VARCHAR(50),
            data_type VARCHAR(50),
            created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
        );
        """)
        
        # Create molecular_properties table
        logger.info("Creating molecular_properties table...")
        cursor.execute("""
        CREATE TABLE IF NOT EXISTS molecular_properties (
            id SERIAL PRIMARY KEY,
            molecule_id INTEGER REFERENCES molecules(id),
            property_type_id UUID REFERENCES property_types(id),
            value_text TEXT,
            value_numeric FLOAT,
            value_boolean BOOLEAN,
            source VARCHAR(100),
            calculation_method VARCHAR(100),
            created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
            updated_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
        );
        """)
        
        # Create mixtures table
        logger.info("Creating mixtures table...")
        cursor.execute("""
        CREATE TABLE IF NOT EXISTS mixtures (
            id SERIAL PRIMARY KEY,
            name VARCHAR(255),
            description TEXT,
            created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
            updated_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
        );
        """)
        
        # Create mixture_components table
        logger.info("Creating mixture_components table...")
        cursor.execute("""
        CREATE TABLE IF NOT EXISTS mixture_components (
            id SERIAL PRIMARY KEY,
            mixture_id INTEGER REFERENCES mixtures(id),
            molecule_id INTEGER REFERENCES molecules(id),
            concentration FLOAT,
            concentration_unit VARCHAR(50),
            created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
        );
        """)
        
        # Create experiments table
        logger.info("Creating experiments table...")
        cursor.execute("""
        CREATE TABLE IF NOT EXISTS experiments (
            id SERIAL PRIMARY KEY,
            name VARCHAR(255),
            description TEXT,
            protocol TEXT,
            user_id VARCHAR(255),
            created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
            updated_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
        );
        """)
        
        # Create experiment_properties table
        logger.info("Creating experiment_properties table...")
        cursor.execute("""
        CREATE TABLE IF NOT EXISTS experiment_properties (
            id SERIAL PRIMARY KEY,
            experiment_id INTEGER REFERENCES experiments(id),
            property_type_id UUID REFERENCES property_types(id),
            value_text TEXT,
            value_numeric FLOAT,
            value_boolean BOOLEAN,
            created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
            updated_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
        );
        """)
        
        # Create users table (simplified version for Heroku)
        logger.info("Creating users table...")
        cursor.execute("""
        CREATE TABLE IF NOT EXISTS users (
            id VARCHAR(255) PRIMARY KEY,
            email VARCHAR(255) UNIQUE,
            name VARCHAR(255),
            role VARCHAR(50) DEFAULT 'user',
            created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
        );
        """)
        
        # Create teams table
        logger.info("Creating teams table...")
        cursor.execute("""
        CREATE TABLE IF NOT EXISTS teams (
            id SERIAL PRIMARY KEY,
            name VARCHAR(255) NOT NULL,
            description TEXT,
            created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
        );
        """)
        
        # Create user_teams table (many-to-many relation)
        logger.info("Creating user_teams table...")
        cursor.execute("""
        CREATE TABLE IF NOT EXISTS user_teams (
            id SERIAL PRIMARY KEY,
            user_id VARCHAR(255) REFERENCES users(id),
            team_id INTEGER REFERENCES teams(id),
            created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
            UNIQUE(user_id, team_id)
        );
        """)

        # Add useful indexes on commonly queried fields
        logger.info("Creating indexes...")
        cursor.execute("CREATE INDEX IF NOT EXISTS idx_molecules_pubchem_cid ON molecules(pubchem_cid);")
        cursor.execute("CREATE INDEX IF NOT EXISTS idx_molecules_name ON molecules(name);")
        cursor.execute("CREATE INDEX IF NOT EXISTS idx_molecular_properties_molecule_id ON molecular_properties(molecule_id);")
        cursor.execute("CREATE INDEX IF NOT EXISTS idx_molecular_properties_property_type_id ON molecular_properties(property_type_id);")
        cursor.execute("CREATE INDEX IF NOT EXISTS idx_mixture_components_mixture_id ON mixture_components(mixture_id);")
        cursor.execute("CREATE INDEX IF NOT EXISTS idx_mixture_components_molecule_id ON mixture_components(molecule_id);")
        
        # Seed basic property types
        logger.info("Adding basic property types...")
        property_types = [
            ("Molecular Weight", "Molecular weight of the compound", "g/mol", "numeric"),
            ("LogP", "Partition coefficient (octanol-water)", None, "numeric"),
            ("Hydrogen Bond Donors", "Number of hydrogen bond donors", None, "numeric"),
            ("Hydrogen Bond Acceptors", "Number of hydrogen bond acceptors", None, "numeric"),
            ("Rotatable Bonds", "Number of rotatable bonds", None, "numeric"),
            ("TPSA", "Topological polar surface area", "Å²", "numeric"),
            ("pKa", "Acid dissociation constant", None, "numeric"),
            ("Solubility", "Solubility in water", "mg/mL", "numeric"),
            ("Melting Point", "Melting point", "°C", "numeric"),
            ("Boiling Point", "Boiling point", "°C", "numeric")
        ]
        
        # Check if property_types already has data
        cursor.execute("SELECT COUNT(*) FROM property_types;")
        property_count = cursor.fetchone()[0]
        
        if property_count == 0:
            for name, description, unit, data_type in property_types:
                cursor.execute(
                    """
                    INSERT INTO property_types (name, description, unit, data_type)
                    VALUES (%s, %s, %s, %s);
                    """,
                    (name, description, unit, data_type)
                )
            logger.info(f"Added {len(property_types)} initial property types.")
        else:
            logger.info(f"Property types table already has {property_count} entries. Skipping seed data.")
            
        cursor.close()
        conn.close()
        
        logger.info("All tables created successfully.")
        return True
    except Exception as e:
        logger.error(f"Error creating tables: {e}")
        return False

def main():
    """Main function to create tables."""
    # Get database URL from environment
    database_url = os.environ.get('DATABASE_URL')
    if not database_url:
        logger.error("DATABASE_URL environment variable not set.")
        return False
    
    # Parse the URL to get connection parameters
    parsed_url = urlparse(database_url)
    connection_params = {
        'host': parsed_url.hostname,
        'port': parsed_url.port or 5432,
        'database': parsed_url.path.lstrip('/'),
        'user': parsed_url.username,
        'password': parsed_url.password
    }
    
    # Create the tables
    success = create_tables(connection_params)
    return success

if __name__ == '__main__':
    success = main()
    if not success:
        logger.error("Table creation failed.")
        sys.exit(1)
    sys.exit(0)