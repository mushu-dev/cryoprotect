#!/usr/bin/env python3
"""
This script checks if the consolidated_molecules table exists and creates it if needed.
"""

import os
import sys
import logging
import psycopg2
import psycopg2.extras
from dotenv import load_dotenv
from datetime import datetime

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler("database_setup.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

# Load environment variables
load_dotenv()

# Database connection parameters
DB_HOST = os.getenv("SUPABASE_DB_HOST")
DB_PORT = os.getenv("SUPABASE_DB_PORT")
DB_NAME = os.getenv("SUPABASE_DB_NAME")
DB_USER = os.getenv("SUPABASE_DB_USER")
DB_PASSWORD = os.getenv("SUPABASE_DB_PASSWORD")

def get_db_connection():
    """Get a direct database connection using psycopg2."""
    try:
        conn = psycopg2.connect(
            host=DB_HOST,
            port=DB_PORT,
            dbname=DB_NAME,
            user=DB_USER,
            password=DB_PASSWORD
        )
        logger.info("Connected to database using direct PostgreSQL connection")
        return conn
    except Exception as e:
        logger.error(f"Database connection error: {e}")
        return None

def check_and_create_consolidated_table():
    """Check if consolidated_molecules table exists, create it if it doesn't."""
    conn = get_db_connection()
    if not conn:
        logger.error("Failed to connect to database")
        return False
        
    cursor = conn.cursor()
    
    try:
        # Check if the table exists
        cursor.execute("""
            SELECT EXISTS (
                SELECT FROM information_schema.tables 
                WHERE table_schema = 'public' 
                AND table_name = 'consolidated_molecules'
            )
        """)
        
        table_exists = cursor.fetchone()[0]
        
        if table_exists:
            logger.info("consolidated_molecules table already exists")
            return True
            
        # Table doesn't exist, create it
        logger.info("Creating consolidated_molecules table...")
        
        cursor.execute("""
            CREATE TABLE public.consolidated_molecules (
                id UUID PRIMARY KEY REFERENCES public.molecules(id),
                name TEXT NOT NULL,
                smiles TEXT,
                inchi TEXT,
                inchikey TEXT,
                formula TEXT,
                molecular_weight DOUBLE PRECISION,
                pubchem_cid TEXT,
                is_public BOOLEAN DEFAULT TRUE,
                data_source TEXT,
                created_at TIMESTAMP WITH TIME ZONE DEFAULT NOW(),
                updated_at TIMESTAMP WITH TIME ZONE DEFAULT NOW(),
                is_consolidated BOOLEAN DEFAULT FALSE,
                primary_molecule_id UUID REFERENCES public.molecules(id),
                primary_molecule_name TEXT,
                primary_pubchem_cid TEXT,
                molecule_status TEXT DEFAULT 'original',
                CONSTRAINT valid_status CHECK (molecule_status IN ('original', 'primary', 'duplicate'))
            );
            
            COMMENT ON TABLE public.consolidated_molecules IS 'Stores consolidated molecule information and mappings between duplicate molecules and their primary records';
            
            -- Create indexes
            CREATE INDEX idx_consolidated_molecules_inchikey ON public.consolidated_molecules(inchikey);
            CREATE INDEX idx_consolidated_molecules_primary_id ON public.consolidated_molecules(primary_molecule_id);
            CREATE INDEX idx_consolidated_molecules_status ON public.consolidated_molecules(molecule_status);
            
            -- Enable Row Level Security
            ALTER TABLE public.consolidated_molecules ENABLE ROW LEVEL SECURITY;
            
            -- Create RLS policies
            CREATE POLICY "Consolidated molecules are viewable by everyone" 
            ON public.consolidated_molecules
            FOR SELECT USING (true);
            
            CREATE POLICY "Consolidated molecules can be updated by authenticated users with proper role" 
            ON public.consolidated_molecules
            FOR UPDATE USING (
                auth.role() = 'authenticated' AND 
                (auth.jwt() ->> 'app_role')::text = 'admin'
            );
            
            CREATE POLICY "Consolidated molecules can be inserted by authenticated users with proper role" 
            ON public.consolidated_molecules
            FOR INSERT WITH CHECK (
                auth.role() = 'authenticated' AND 
                (auth.jwt() ->> 'app_role')::text = 'admin'
            );
        """)
        
        # Create scientific_data_audit table if it doesn't exist
        cursor.execute("""
            SELECT EXISTS (
                SELECT FROM information_schema.tables 
                WHERE table_schema = 'public' 
                AND table_name = 'scientific_data_audit'
            )
        """)
        
        audit_table_exists = cursor.fetchone()[0]
        
        if not audit_table_exists:
            logger.info("Creating scientific_data_audit table...")
            
            cursor.execute("""
                CREATE TABLE public.scientific_data_audit (
                    id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
                    table_name TEXT NOT NULL,
                    record_id TEXT NOT NULL,
                    operation TEXT NOT NULL,
                    old_value JSONB,
                    new_value JSONB,
                    user_id UUID,
                    timestamp TIMESTAMP WITH TIME ZONE DEFAULT NOW()
                );
                
                COMMENT ON TABLE public.scientific_data_audit IS 'Audit trail for scientific data operations';
                
                -- Create indexes
                CREATE INDEX idx_scientific_data_audit_table ON public.scientific_data_audit(table_name);
                CREATE INDEX idx_scientific_data_audit_record ON public.scientific_data_audit(record_id);
                CREATE INDEX idx_scientific_data_audit_timestamp ON public.scientific_data_audit(timestamp);
                
                -- Enable Row Level Security
                ALTER TABLE public.scientific_data_audit ENABLE ROW LEVEL SECURITY;
                
                -- Create RLS policies
                CREATE POLICY "Scientific data audit is viewable by authenticated users" 
                ON public.scientific_data_audit
                FOR SELECT USING (auth.role() = 'authenticated');
                
                CREATE POLICY "Scientific data audit can be inserted by service role" 
                ON public.scientific_data_audit
                FOR INSERT WITH CHECK (true);
            """)
        
        conn.commit()
        logger.info("Successfully created necessary tables")
        return True
        
    except Exception as e:
        conn.rollback()
        logger.error(f"Error creating consolidated_molecules table: {e}")
        return False
    finally:
        cursor.close()
        conn.close()

if __name__ == "__main__":
    success = check_and_create_consolidated_table()
    exit(0 if success else 1)