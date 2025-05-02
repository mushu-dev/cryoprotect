#!/usr/bin/env python3
"""
CryoProtect v2 - Fix Relationship Design Issues

This script fixes relationship design issues in the CryoProtect Supabase project:
1. Replaces fan traps with proper junction tables
2. Applies 3NF normalization principles
3. Adds missing foreign key constraints with appropriate indexes
4. Includes verification and rollback mechanisms

Usage:
    python fix_relationships.py [--dry-run] [--verify-only] [--rollback]
"""

import os
import sys
import json
import time
import argparse
import logging
from datetime import datetime
from dotenv import load_dotenv

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler("fix_relationships.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

# Load environment variables
load_dotenv()

# Supabase project ID
SUPABASE_PROJECT_ID = "tsdlmynydfuypiugmkev"

# Track changes for potential rollback
applied_changes = []

def connect_to_supabase():
    """Connect to Supabase using MCP tools."""
    try:
        import json
        from supabase import create_client, Client
        
        # Get Supabase URL and key from environment variables
        supabase_url = os.getenv("SUPABASE_URL")
        supabase_key = os.getenv("SUPABASE_KEY")
        
        if not supabase_url or not supabase_key:
            logger.error("SUPABASE_URL and SUPABASE_KEY must be set in .env file")
            sys.exit(1)
        
        # Connect to Supabase
        supabase = create_client(supabase_url, supabase_key)
        logger.info("Successfully connected to Supabase")
        return supabase
    except Exception as e:
        logger.error(f"Error connecting to Supabase: {str(e)}")
        sys.exit(1)

def execute_sql(supabase, sql, description):
    """Execute SQL using Supabase MCP tools."""
    try:
        logger.info(f"Executing SQL: {description}")
        response = supabase.rpc("exec_sql", {"sql": sql}).execute()
        logger.info(f"SQL executed successfully: {description}")
        
        # Track the change for potential rollback
        applied_changes.append({
            "description": description,
            "timestamp": datetime.now().isoformat(),
            "sql": sql
        })
        
        return True
    except Exception as e:
        logger.error(f"Error executing SQL ({description}): {str(e)}")
        return False

def check_table_exists(supabase, table_name):
    """Check if a table exists in the database."""
    try:
        # Try to select a single row from the table
        response = supabase.table(table_name).select("*").limit(1).execute()
        # If we get here without an error, the table exists
        logger.info(f"Table '{table_name}' exists")
        return True
    except Exception as e:
        if "relation" in str(e) and "does not exist" in str(e):
            logger.warning(f"Table '{table_name}' does not exist")
            return False
        # If it's some other error, log it and assume the table doesn't exist
        logger.error(f"Error checking if table {table_name} exists: {str(e)}")
        return False

def backup_database(supabase):
    """Create a backup of the current database state."""
    logger.info("Creating database backup before making changes...")
    
    tables_to_backup = [
        "molecules", 
        "mixtures", 
        "mixture_components", 
        "predictions", 
        "experiments", 
        "property_types", 
        "calculation_methods",
        "molecular_properties"
    ]
    
    backup_data = {}
    
    for table in tables_to_backup:
        try:
            if check_table_exists(supabase, table):
                response = supabase.table(table).select("*").execute()
                if hasattr(response, 'data'):
                    backup_data[table] = response.data
                    logger.info(f"Backed up {len(response.data)} records from '{table}'")
                else:
                    backup_data[table] = []
                    logger.warning(f"No data found in '{table}'")
            else:
                backup_data[table] = []
                logger.warning(f"Table '{table}' does not exist, skipping backup")
        except Exception as e:
            logger.error(f"Error backing up table '{table}': {str(e)}")
            return None
    
    # Save backup to file
    backup_filename = f"database_backup_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json"
    try:
        with open(backup_filename, 'w') as f:
            json.dump(backup_data, f, indent=2)
        logger.info(f"Database backup saved to {backup_filename}")
        return backup_filename
    except Exception as e:
        logger.error(f"Error saving backup to file: {str(e)}")
        return None

def create_molecule_proteins_junction_table(supabase, dry_run=False):
    """Create the molecule_proteins junction table."""
    if check_table_exists(supabase, "molecule_proteins") and not dry_run:
        logger.info("Table 'molecule_proteins' already exists.")
        return True
    
    logger.info("Creating 'molecule_proteins' junction table...")
    
    sql = """
    CREATE TABLE IF NOT EXISTS public.molecule_proteins (
        id UUID PRIMARY KEY DEFAULT uuid_generate_v4(),
        molecule_id UUID NOT NULL REFERENCES public.molecules(id) ON DELETE CASCADE,
        protein_id UUID NOT NULL REFERENCES public.proteins(id) ON DELETE CASCADE,
        binding_affinity NUMERIC,
        interaction_type VARCHAR(100),
        created_at TIMESTAMPTZ NOT NULL DEFAULT NOW(),
        updated_at TIMESTAMPTZ NOT NULL DEFAULT NOW(),
        created_by UUID REFERENCES auth.users(id),
        UNIQUE(molecule_id, protein_id)
    );
    
    CREATE INDEX IF NOT EXISTS idx_molecule_proteins_molecule_id ON public.molecule_proteins(molecule_id);
    CREATE INDEX IF NOT EXISTS idx_molecule_proteins_protein_id ON public.molecule_proteins(protein_id);
    
    CREATE TRIGGER set_timestamp_molecule_proteins
    BEFORE UPDATE ON public.molecule_proteins
    FOR EACH ROW EXECUTE PROCEDURE public.trigger_set_timestamp();
    """
    
    if dry_run:
        logger.info("[DRY RUN] Would create 'molecule_proteins' junction table")
        return True
    
    return execute_sql(supabase, sql, "Create molecule_proteins junction table")

def create_molecule_experiments_junction_table(supabase, dry_run=False):
    """Create the molecule_experiments junction table."""
    if check_table_exists(supabase, "molecule_experiments") and not dry_run:
        logger.info("Table 'molecule_experiments' already exists.")
        return True
    
    logger.info("Creating 'molecule_experiments' junction table...")
    
    sql = """
    CREATE TABLE IF NOT EXISTS public.molecule_experiments (
        id UUID PRIMARY KEY DEFAULT uuid_generate_v4(),
        molecule_id UUID NOT NULL REFERENCES public.molecules(id) ON DELETE CASCADE,
        experiment_id UUID NOT NULL REFERENCES public.experiments(id) ON DELETE CASCADE,
        role VARCHAR(100),
        concentration NUMERIC,
        concentration_unit VARCHAR(50),
        created_at TIMESTAMPTZ NOT NULL DEFAULT NOW(),
        updated_at TIMESTAMPTZ NOT NULL DEFAULT NOW(),
        created_by UUID REFERENCES auth.users(id),
        UNIQUE(molecule_id, experiment_id)
    );
    
    CREATE INDEX IF NOT EXISTS idx_molecule_experiments_molecule_id ON public.molecule_experiments(molecule_id);
    CREATE INDEX IF NOT EXISTS idx_molecule_experiments_experiment_id ON public.molecule_experiments(experiment_id);
    
    CREATE TRIGGER set_timestamp_molecule_experiments
    BEFORE UPDATE ON public.molecule_experiments
    FOR EACH ROW EXECUTE PROCEDURE public.trigger_set_timestamp();
    """
    
    if dry_run:
        logger.info("[DRY RUN] Would create 'molecule_experiments' junction table")
        return True
    
    return execute_sql(supabase, sql, "Create molecule_experiments junction table")

def fix_predictions_table(supabase, dry_run=False):
    """Fix the predictions table to properly handle molecule and mixture relationships."""
    logger.info("Fixing 'predictions' table to handle both molecule and mixture relationships...")
    
    sql = """
    -- First, make sure the trigger_set_timestamp function exists
    CREATE OR REPLACE FUNCTION public.trigger_set_timestamp()
    RETURNS TRIGGER AS $$
    BEGIN
      NEW.updated_at = NOW();
      RETURN NEW;
    END;
    $$ LANGUAGE plpgsql;

    -- Modify the predictions table to allow either molecule_id or mixture_id (but not both null)
    ALTER TABLE public.predictions 
    DROP CONSTRAINT IF EXISTS predictions_mixture_id_fkey,
    ALTER COLUMN mixture_id DROP NOT NULL,
    ADD COLUMN IF NOT EXISTS molecule_id UUID REFERENCES public.molecules(id) ON DELETE CASCADE,
    DROP CONSTRAINT IF EXISTS check_prediction_target,
    ADD CONSTRAINT check_prediction_target CHECK (
        (molecule_id IS NOT NULL AND mixture_id IS NULL) OR 
        (molecule_id IS NULL AND mixture_id IS NOT NULL)
    );
    
    -- Re-add the mixture_id foreign key
    ALTER TABLE public.predictions
    ADD CONSTRAINT predictions_mixture_id_fkey 
    FOREIGN KEY (mixture_id) REFERENCES public.mixtures(id) ON DELETE CASCADE;
    
    -- Add index on molecule_id
    CREATE INDEX IF NOT EXISTS idx_predictions_molecule_id ON public.predictions(molecule_id);
    
    -- Update the unique constraint to account for either molecule_id or mixture_id
    ALTER TABLE public.predictions
    DROP CONSTRAINT IF EXISTS predictions_mixture_id_property_type_id_calculation_method_id_key,
    ADD CONSTRAINT predictions_unique_constraint UNIQUE (
        COALESCE(molecule_id, uuid_nil()),
        COALESCE(mixture_id, uuid_nil()),
        property_type_id,
        calculation_method_id
    );
    """
    
    if dry_run:
        logger.info("[DRY RUN] Would fix 'predictions' table")
        return True
    
    return execute_sql(supabase, sql, "Fix predictions table")

def fix_experiments_table(supabase, dry_run=False):
    """Fix the experiments table to add missing constraints."""
    logger.info("Fixing 'experiments' table to add missing constraints...")
    
    sql = """
    -- Add unique constraint to experiments table
    ALTER TABLE public.experiments
    DROP CONSTRAINT IF EXISTS experiments_unique_constraint,
    ADD CONSTRAINT experiments_unique_constraint UNIQUE (
        mixture_id,
        property_type_id,
        date_performed
    );
    """
    
    if dry_run:
        logger.info("[DRY RUN] Would fix 'experiments' table")
        return True
    
    return execute_sql(supabase, sql, "Fix experiments table")

def fix_molecular_properties_table(supabase, dry_run=False):
    """Fix the molecular_properties table to use property_type_id instead of property_type string."""
    logger.info("Fixing 'molecular_properties' table to use property_type_id...")
    
    # First check if property_types table exists
    if not check_table_exists(supabase, "property_types") and not dry_run:
        logger.info("Creating 'property_types' table...")
        
        create_property_types_sql = """
        CREATE TABLE IF NOT EXISTS public.property_types (
            id UUID PRIMARY KEY DEFAULT uuid_generate_v4(),
            name VARCHAR(255) NOT NULL UNIQUE,
            data_type VARCHAR(50) NOT NULL,
            description TEXT,
            created_at TIMESTAMPTZ NOT NULL DEFAULT NOW(),
            updated_at TIMESTAMPTZ NOT NULL DEFAULT NOW()
        );
        
        CREATE TRIGGER set_timestamp_property_types
        BEFORE UPDATE ON public.property_types
        FOR EACH ROW EXECUTE PROCEDURE public.trigger_set_timestamp();
        """
        
        if dry_run:
            logger.info("[DRY RUN] Would create 'property_types' table")
        else:
            if not execute_sql(supabase, create_property_types_sql, "Create property_types table"):
                return False
    
    # Now modify the molecular_properties table
    sql = """
    -- First, migrate existing property types to the property_types table
    INSERT INTO public.property_types (name, data_type)
    SELECT DISTINCT property_type, 'numeric'
    FROM public.molecular_properties
    WHERE property_type NOT IN (SELECT name FROM public.property_types)
    ON CONFLICT (name) DO NOTHING;
    
    -- Add property_type_id column
    ALTER TABLE public.molecular_properties
    ADD COLUMN IF NOT EXISTS property_type_id UUID;
    
    -- Update property_type_id based on property_type
    UPDATE public.molecular_properties mp
    SET property_type_id = pt.id
    FROM public.property_types pt
    WHERE mp.property_type = pt.name
    AND mp.property_type_id IS NULL;
    
    -- Make property_type_id NOT NULL and add foreign key constraint
    ALTER TABLE public.molecular_properties
    ALTER COLUMN property_type_id SET NOT NULL,
    ADD CONSTRAINT molecular_properties_property_type_id_fkey
    FOREIGN KEY (property_type_id) REFERENCES public.property_types(id);
    
    -- Add index on property_type_id
    CREATE INDEX IF NOT EXISTS idx_molecular_properties_property_type_id 
    ON public.molecular_properties(property_type_id);
    
    -- Add unique constraint
    ALTER TABLE public.molecular_properties
    DROP CONSTRAINT IF EXISTS molecular_properties_unique_constraint,
    ADD CONSTRAINT molecular_properties_unique_constraint UNIQUE (
        molecule_id,
        property_type_id
    );
    """
    
    if dry_run:
        logger.info("[DRY RUN] Would fix 'molecular_properties' table")
        return True
    
    return execute_sql(supabase, sql, "Fix molecular_properties table")

def migrate_data_to_junction_tables(supabase, dry_run=False):
    """Migrate data to the new junction tables."""
    logger.info("Migrating data to junction tables...")
    
    # Check if proteins table exists, if not, we can't migrate to molecule_proteins
    if not check_table_exists(supabase, "proteins") and not dry_run:
        logger.warning("Table 'proteins' does not exist. Cannot migrate data to molecule_proteins.")
    else:
        # We would add code here to migrate data to molecule_proteins if needed
        pass
    
    # Migrate data to molecule_experiments
    migrate_to_molecule_experiments_sql = """
    -- Insert data into molecule_experiments from experiments and mixture_components
    INSERT INTO public.molecule_experiments (molecule_id, experiment_id, role, concentration, concentration_unit)
    SELECT mc.molecule_id, e.id, 'component', mc.amount, mc.amount_unit
    FROM public.experiments e
    JOIN public.mixture_components mc ON e.mixture_id = mc.mixture_id
    WHERE NOT EXISTS (
        SELECT 1 FROM public.molecule_experiments me
        WHERE me.molecule_id = mc.molecule_id AND me.experiment_id = e.id
    );
    """
    
    if dry_run:
        logger.info("[DRY RUN] Would migrate data to junction tables")
        return True
    
    if check_table_exists(supabase, "molecule_experiments") and not dry_run:
        if not execute_sql(supabase, migrate_to_molecule_experiments_sql, "Migrate data to molecule_experiments"):
            return False
    
    return True

def verify_relationships(supabase):
    """Verify that the relationships are correctly set up."""
    logger.info("Verifying relationships...")
    
    verification_results = {}
    all_valid = True
    
    # Check if junction tables exist
    for table in ["molecule_proteins", "molecule_experiments"]:
        exists = check_table_exists(supabase, table)
        verification_results[f"{table}_exists"] = exists
        if not exists:
            all_valid = False
    
    # Verify predictions table structure
    predictions_sql = """
    SELECT column_name, data_type, is_nullable
    FROM information_schema.columns
    WHERE table_schema = 'public' AND table_name = 'predictions'
    ORDER BY ordinal_position;
    """
    
    try:
        response = supabase.rpc("exec_sql", {"sql": predictions_sql}).execute()
        if hasattr(response, 'data'):
            columns = {col['column_name']: col for col in response.data}
            
            # Check if molecule_id exists
            if 'molecule_id' not in columns:
                verification_results["predictions_has_molecule_id"] = False
                all_valid = False
            else:
                verification_results["predictions_has_molecule_id"] = True
            
            # Check if mixture_id is nullable
            if 'mixture_id' in columns and columns['mixture_id']['is_nullable'] != 'YES':
                verification_results["predictions_mixture_id_nullable"] = False
                all_valid = False
            else:
                verification_results["predictions_mixture_id_nullable"] = True
        else:
            verification_results["predictions_structure_check"] = False
            all_valid = False
    except Exception as e:
        logger.error(f"Error verifying predictions table structure: {str(e)}")
        verification_results["predictions_structure_check"] = False
        all_valid = False
    
    # Verify molecular_properties table structure
    molecular_properties_sql = """
    SELECT column_name, data_type, is_nullable
    FROM information_schema.columns
    WHERE table_schema = 'public' AND table_name = 'molecular_properties'
    ORDER BY ordinal_position;
    """
    
    try:
        response = supabase.rpc("exec_sql", {"sql": molecular_properties_sql}).execute()
        if hasattr(response, 'data'):
            columns = {col['column_name']: col for col in response.data}
            
            # Check if property_type_id exists
            if 'property_type_id' not in columns:
                verification_results["molecular_properties_has_property_type_id"] = False
                all_valid = False
            else:
                verification_results["molecular_properties_has_property_type_id"] = True
        else:
            verification_results["molecular_properties_structure_check"] = False
            all_valid = False
    except Exception as e:
        logger.error(f"Error verifying molecular_properties table structure: {str(e)}")
        verification_results["molecular_properties_structure_check"] = False
        all_valid = False
    
    # Verify constraints
    constraints_sql = """
    SELECT conname, contype, conrelid::regclass AS table_name
    FROM pg_constraint
    WHERE connamespace = 'public'::regnamespace
    AND conrelid::regclass::text IN ('predictions', 'experiments', 'molecular_properties');
    """
    
    try:
        response = supabase.rpc("exec_sql", {"sql": constraints_sql}).execute()
        if hasattr(response, 'data'):
            constraints = response.data
            
            # Check for check constraint on predictions
            has_check_constraint = any(
                c['conname'] == 'check_prediction_target' and c['table_name'] == 'predictions' 
                for c in constraints
            )
            verification_results["predictions_has_check_constraint"] = has_check_constraint
            if not has_check_constraint:
                all_valid = False
            
            # Check for unique constraints
            tables_with_unique = set(
                c['table_name'] for c in constraints if c['contype'] == 'u'
            )
            
            for table in ['predictions', 'experiments', 'molecular_properties']:
                has_unique = table in tables_with_unique
                verification_results[f"{table}_has_unique_constraint"] = has_unique
                if not has_unique:
                    all_valid = False
        else:
            verification_results["constraints_check"] = False
            all_valid = False
    except Exception as e:
        logger.error(f"Error verifying constraints: {str(e)}")
        verification_results["constraints_check"] = False
        all_valid = False
    
    return all_valid, verification_results

def rollback_changes(supabase, backup_filename):
    """Rollback changes using the backup file."""
    logger.info(f"Rolling back changes using backup file: {backup_filename}")
    
    try:
        with open(backup_filename, 'r') as f:
            backup_data = json.load(f)
        
        # Drop the new tables first
        for table in ["molecule_proteins", "molecule_experiments"]:
            if check_table_exists(supabase, table):
                drop_sql = f"DROP TABLE IF EXISTS public.{table} CASCADE;"
                execute_sql(supabase, drop_sql, f"Drop {table} table")
        
        # Restore the original tables
        for table, data in backup_data.items():
            if check_table_exists(supabase, table) and data:
                # Clear the table
                truncate_sql = f"TRUNCATE TABLE public.{table} CASCADE;"
                execute_sql(supabase, truncate_sql, f"Truncate {table} table")
                
                # Restore data
                for record in data:
                    try:
                        supabase.table(table).upsert(record).execute()
                    except Exception as e:
                        logger.error(f"Error restoring record to {table}: {str(e)}")
        
        logger.info("Rollback completed successfully")
        return True
    except Exception as e:
        logger.error(f"Error during rollback: {str(e)}")
        return False

def main():
    """Main function to fix relationship design issues."""
    parser = argparse.ArgumentParser(description='Fix relationship design issues in CryoProtect Supabase project')
    parser.add_argument('--dry-run', action='store_true', help='Show what would be done without making changes')
    parser.add_argument('--verify-only', action='store_true', help='Only verify the relationships without making changes')
    parser.add_argument('--rollback', action='store_true', help='Rollback changes using the most recent backup')
    args = parser.parse_args()
    
    print("\n" + "=" * 80)
    print(f"CryoProtect v2 - Fix Relationship Design Issues (Project ID: {SUPABASE_PROJECT_ID})")
    print("=" * 80)
    
    try:
        # Connect to Supabase
        print("\nConnecting to Supabase...")
        supabase = connect_to_supabase()
        
        # If rollback is requested, find the most recent backup and rollback
        if args.rollback:
            backup_files = [f for f in os.listdir('.') if f.startswith('database_backup_') and f.endswith('.json')]
            if not backup_files:
                logger.error("No backup files found for rollback")
                return 1
            
            most_recent_backup = sorted(backup_files)[-1]
            print(f"\nRolling back changes using backup: {most_recent_backup}")
            
            if rollback_changes(supabase, most_recent_backup):
                print("Rollback completed successfully")
                return 0
            else:
                print("Rollback failed. Check the logs for details.")
                return 1
        
        # If verify-only is requested, just verify the relationships
        if args.verify_only:
            print("\nVerifying relationships...")
            all_valid, verification_results = verify_relationships(supabase)
            
            print("\nVerification Results:")
            for check, result in verification_results.items():
                print(f"  - {check}: {'Valid' if result else 'Invalid'}")
            
            if all_valid:
                print("\nAll relationships are valid.")
                return 0
            else:
                print("\nSome relationships are invalid. Check the logs for details.")
                return 1
        
        # Create a backup before making changes
        if not args.dry_run:
            backup_filename = backup_database(supabase)
            if not backup_filename:
                logger.error("Failed to create database backup. Aborting.")
                return 1
        
        # Fix relationship design issues
        success = True
        
        # 1. Create junction tables
        print("\nCreating junction tables...")
        if not create_molecule_proteins_junction_table(supabase, args.dry_run):
            success = False
        
        if not create_molecule_experiments_junction_table(supabase, args.dry_run):
            success = False
        
        # 2. Fix existing tables
        print("\nFixing existing tables...")
        if not fix_predictions_table(supabase, args.dry_run):
            success = False
        
        if not fix_experiments_table(supabase, args.dry_run):
            success = False
        
        if not fix_molecular_properties_table(supabase, args.dry_run):
            success = False
        
        # 3. Migrate data to junction tables
        print("\nMigrating data to junction tables...")
        if not migrate_data_to_junction_tables(supabase, args.dry_run):
            success = False
        
        # 4. Verify the changes
        if not args.dry_run and success:
            print("\nVerifying relationships...")
            all_valid, verification_results = verify_relationships(supabase)
            
            print("\nVerification Results:")
            for check, result in verification_results.items():
                print(f"  - {check}: {'Valid' if result else 'Invalid'}")
            
            if all_valid:
                print("\nAll relationships are valid.")
            else:
                print("\nSome relationships are invalid. Check the logs for details.")
                print("Consider using --rollback to revert changes.")
                success = False
        
        if args.dry_run:
            print("\n[DRY RUN] No changes were made to the database.")
        elif success:
            print("\n" + "=" * 60)
            print("Relationship Design Fixes Complete")
            print("=" * 60)
            print("\nSuccessfully fixed the following issues:")
            print("1. Created molecule_proteins junction table")
            print("2. Created molecule_experiments junction table")
            print("3. Fixed predictions table to handle both molecules and mixtures")
            print("4. Added missing constraints to experiments table")
            print("5. Fixed molecular_properties table to use property_type_id")
            print("6. Migrated data to junction tables")
            print("\nAll relationships are now properly designed according to 3NF principles.")
        else:
            print("\nFailed to fix some relationship design issues. Check the logs for details.")
            print("Consider using --rollback to revert changes.")
            return 1
        
        return 0
    
    except Exception as e:
        logger.error(f"Error fixing relationship design issues: {str(e)}")
        print(f"\nError fixing relationship design issues: {str(e)}")
        return 1

if __name__ == "__main__":
    sys.exit(main())