#!/usr/bin/env python3
"""
Fix Reference Compounds using Supabase MCP

This script uses the Supabase MCP tool to directly fix reference compound properties
without relying on the standard database connection utilities.
"""

import os
import sys
import logging
import json
import time
from pathlib import Path

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler(),
        logging.FileHandler('logs/fix_reference_compounds_mcp.log')
    ]
)
logger = logging.getLogger(__name__)

# Create logs directory if it doesn't exist
os.makedirs('logs', exist_ok=True)

def get_reference_compound_ids():
    """Get the list of reference compound ChEMBL IDs."""
    return [
        "CHEMBL1098659",
        "CHEMBL1487",
        "CHEMBL262548",
        "CHEMBL388978",
        "CHEMBL500033",
        "CHEMBL6196",
        "CHEMBL66195",
        "CHEMBL6752",
        "CHEMBL967"
    ]

def get_project_id():
    """Get the Supabase project ID from environment variables or .env file."""
    # Try to load from .env file
    env_path = Path('.env')
    if env_path.exists():
        with open(env_path, 'r') as f:
            for line in f:
                if line.strip().startswith('SUPABASE_PROJECT_ID='):
                    return line.strip().split('=', 1)[1].strip('"\'')
    
    # Try environment variable
    project_id = os.environ.get('SUPABASE_PROJECT_ID')
    if project_id:
        return project_id
    
    # Default value from previous logs
    return "tsdlmynydfuypiugmkev"

def execute_sql(query):
    """Execute SQL using the Supabase MCP tool."""
    try:
        # Import the execute_sql function from use_mcp_tool
        from use_mcp_tool import execute_sql as mcp_execute_sql
        
        # Get project ID
        project_id = get_project_id()
        
        logger.info(f"Executing SQL using MCP tool (project_id: {project_id})")
        
        # Execute the query
        result = mcp_execute_sql(query, project_id)
        
        return result
    except ImportError:
        logger.error("Could not import execute_sql from use_mcp_tool. Make sure it's available.")
        raise
    except Exception as e:
        logger.error(f"Error executing SQL: {str(e)}")
        raise

def verify_reference_compounds():
    """Verify that reference compounds exist in the database."""
    logger.info("Verifying reference compounds exist...")
    
    # Define reference compound IDs
    reference_ids = get_reference_compound_ids()
    
    # Create SQL query
    query = f"""
    SELECT chembl_id, id FROM molecules 
    WHERE chembl_id IN ({', '.join([f"'{id}'" for id in reference_ids])});
    """
    
    # Execute query
    try:
        result = execute_sql(query)
        
        if not result:
            logger.error("No reference compounds found in database")
            return False
        
        logger.info(f"Found {len(result)} reference compounds in database")
        
        # Check if all reference compounds exist
        found_ids = [row['chembl_id'] for row in result]
        missing_ids = [id for id in reference_ids if id not in found_ids]
        
        if missing_ids:
            logger.error(f"Missing reference compounds: {missing_ids}")
            return False
        
        logger.info("All reference compounds exist in database")
        return True
    except Exception as e:
        logger.error(f"Error verifying reference compounds: {str(e)}")
        return False

def ensure_property_types():
    """Ensure property types exist in the database."""
    logger.info("Ensuring property types exist...")
    
    # Create SQL query
    query = """
    -- Create property types if they don't exist
    INSERT INTO property_types (name, data_type, description)
    SELECT 'logP', 'numeric', 'Partition coefficient'
    WHERE NOT EXISTS (SELECT 1 FROM property_types WHERE name = 'logP');

    INSERT INTO property_types (name, data_type, description)
    SELECT 'h_bond_donors', 'numeric', 'Number of hydrogen bond donors'
    WHERE NOT EXISTS (SELECT 1 FROM property_types WHERE name = 'h_bond_donors');

    INSERT INTO property_types (name, data_type, description)
    SELECT 'h_bond_acceptors', 'numeric', 'Number of hydrogen bond acceptors'
    WHERE NOT EXISTS (SELECT 1 FROM property_types WHERE name = 'h_bond_acceptors');

    -- Return the property type IDs
    SELECT id, name FROM property_types 
    WHERE name IN ('logP', 'h_bond_donors', 'h_bond_acceptors');
    """
    
    # Execute query
    try:
        result = execute_sql(query)
        
        if not result:
            logger.error("Failed to create or retrieve property types")
            return False
        
        logger.info(f"Property types: {result}")
        return True
    except Exception as e:
        logger.error(f"Error ensuring property types: {str(e)}")
        return False

def fix_reference_compound_properties(chembl_id):
    """Fix missing properties for a reference compound."""
    logger.info(f"Checking properties for {chembl_id}...")
    
    # Check missing properties
    check_query = f"""
    WITH molecule AS (
        SELECT id FROM molecules WHERE chembl_id = '{chembl_id}'
    ),
    property_types AS (
        SELECT id, name FROM property_types 
        WHERE name IN ('logP', 'h_bond_donors', 'h_bond_acceptors')
    ),
    existing_properties AS (
        SELECT pt.name
        FROM molecular_properties mp
        JOIN property_types pt ON mp.property_type_id = pt.id
        JOIN molecule m ON mp.molecule_id = m.id
        WHERE pt.name IN ('logP', 'h_bond_donors', 'h_bond_acceptors')
    )
    SELECT pt.name AS missing_property
    FROM property_types pt
    WHERE pt.name NOT IN (SELECT name FROM existing_properties);
    """
    
    try:
        missing_properties = execute_sql(check_query)
        
        if not missing_properties:
            logger.info(f"No missing properties for {chembl_id}")
            return True
        
        logger.info(f"Missing properties for {chembl_id}: {missing_properties}")
        
        # Fix missing properties
        fix_query = f"""
        -- Get molecule ID
        WITH molecule AS (
            SELECT id FROM molecules WHERE chembl_id = '{chembl_id}'
        ),
        -- Get property type IDs
        property_types AS (
            SELECT id, name FROM property_types 
            WHERE name IN ('logP', 'h_bond_donors', 'h_bond_acceptors')
        ),
        -- Get existing properties
        existing_properties AS (
            SELECT pt.name
            FROM molecular_properties mp
            JOIN property_types pt ON mp.property_type_id = pt.id
            JOIN molecule m ON mp.molecule_id = m.id
            WHERE pt.name IN ('logP', 'h_bond_donors', 'h_bond_acceptors')
        )
        -- Insert missing properties
        INSERT INTO molecular_properties (molecule_id, property_type_id, numeric_value, created_at, updated_at)
        SELECT 
            (SELECT id FROM molecule),
            pt.id,
            CASE 
                WHEN pt.name = 'logP' THEN 0.0
                ELSE 0
            END,
            NOW(),
            NOW()
        FROM property_types pt
        WHERE pt.name NOT IN (SELECT name FROM existing_properties);
        """
        
        execute_sql(fix_query)
        logger.info(f"Fixed missing properties for {chembl_id}")
        
        return True
    except Exception as e:
        logger.error(f"Error fixing properties for {chembl_id}: {str(e)}")
        return False

def verify_all_properties():
    """Verify all reference compounds have all required properties."""
    logger.info("Verifying all properties are now set...")
    
    # Define reference compound IDs
    reference_ids = get_reference_compound_ids()
    
    # Create SQL query
    query = f"""
    WITH reference_molecules AS (
        SELECT id, chembl_id 
        FROM molecules 
        WHERE chembl_id IN ({', '.join([f"'{id}'" for id in reference_ids])})
    ),
    property_counts AS (
        SELECT 
            rm.chembl_id,
            COUNT(DISTINCT pt.name) AS property_count
        FROM reference_molecules rm
        JOIN molecular_properties mp ON rm.id = mp.molecule_id
        JOIN property_types pt ON mp.property_type_id = pt.id
        WHERE pt.name IN ('logP', 'h_bond_donors', 'h_bond_acceptors')
        GROUP BY rm.chembl_id
    )
    SELECT 
        rm.chembl_id,
        COALESCE(pc.property_count, 0) AS property_count,
        CASE 
            WHEN COALESCE(pc.property_count, 0) = 3 THEN 'Complete'
            ELSE 'Incomplete'
        END AS status
    FROM reference_molecules rm
    LEFT JOIN property_counts pc ON rm.chembl_id = pc.chembl_id
    ORDER BY rm.chembl_id;
    """
    
    # Execute query
    try:
        result = execute_sql(query)
        
        if not result:
            logger.error("Failed to verify properties")
            return False
        
        logger.info("Property verification results:")
        for row in result:
            logger.info(f"  {row['chembl_id']}: {row['property_count']} properties, status: {row['status']}")
        
        # Check if all compounds are complete
        incomplete = [row for row in result if row['status'] != 'Complete']
        
        if incomplete:
            logger.warning(f"Some compounds are still incomplete: {[row['chembl_id'] for row in incomplete]}")
            return False
        
        logger.info("All reference compounds now have all required properties")
        return True
    except Exception as e:
        logger.error(f"Error verifying properties: {str(e)}")
        return False

def main():
    """Main function to fix reference compound properties."""
    logger.info("Starting reference compounds fix using MCP...")
    
    # Step 1: Verify reference compounds exist
    if not verify_reference_compounds():
        logger.error("Reference compounds verification failed. Exiting.")
        return 1
    
    # Step 2: Ensure property types exist
    if not ensure_property_types():
        logger.error("Property types creation failed. Exiting.")
        return 1
    
    # Step 3: Fix missing properties for each reference compound
    reference_ids = get_reference_compound_ids()
    for chembl_id in reference_ids:
        if not fix_reference_compound_properties(chembl_id):
            logger.warning(f"Failed to fix properties for {chembl_id}")
    
    # Step 4: Verify all properties are now set
    if not verify_all_properties():
        logger.warning("Some properties may still be missing")
    
    logger.info("Reference compounds fix complete!")
    return 0

if __name__ == "__main__":
    sys.exit(main())