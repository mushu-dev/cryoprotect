#!/usr/bin/env python3
"""
Direct fix for reference compound properties.

This script directly inserts the missing properties for reference compounds
using SQL queries, bypassing the problematic utility functions.
"""

import os
import sys
import logging
import uuid
import psycopg2
from psycopg2.extras import RealDictCursor
from dotenv import load_dotenv

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler(),
        logging.FileHandler('logs/fix_reference_properties.log')
    ]
)
logger = logging.getLogger(__name__)

# Create logs directory if it doesn't exist
os.makedirs('logs', exist_ok=True)

# Load environment variables
load_dotenv()

def get_db_connection():
    """Get a direct database connection."""
    conn = psycopg2.connect(
        host=os.getenv('SUPABASE_DB_HOST'),
        port=os.getenv('SUPABASE_DB_PORT', '5432'),
        dbname=os.getenv('SUPABASE_DB_NAME'),
        user=os.getenv('SUPABASE_DB_USER'),
        password=os.getenv('SUPABASE_DB_PASSWORD')
    )
    conn.autocommit = False
    return conn

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

def get_or_create_property_type(conn, property_name, data_type='numeric'):
    """Get or create a property type."""
    cursor = conn.cursor(cursor_factory=RealDictCursor)
    
    # Check if property type exists
    cursor.execute(
        "SELECT id FROM property_types WHERE name = %s",
        (property_name,)
    )
    result = cursor.fetchone()
    
    if result:
        logger.info(f"Found existing property type: {property_name} (ID: {result['id']})")
        return result['id']
    
    # Create new property type
    property_id = str(uuid.uuid4())
    cursor.execute(
        "INSERT INTO property_types (id, name, data_type, description) VALUES (%s, %s, %s, %s)",
        (property_id, property_name, data_type, f"Auto-created property: {property_name}")
    )
    conn.commit()
    
    logger.info(f"Created new property type: {property_name} (ID: {property_id})")
    return property_id

def fix_reference_compounds_properties():
    """Fix missing properties for reference compounds by directly inserting them."""
    logger.info("Fixing missing properties for reference compounds...")
    
    try:
        # Get reference compound IDs
        reference_ids = get_reference_compound_ids()
        logger.info(f"Found {len(reference_ids)} reference compound IDs")
        
        # Connect to database
        conn = get_db_connection()
        cursor = conn.cursor(cursor_factory=RealDictCursor)
        
        # Get property type IDs
        logp_type_id = get_or_create_property_type(conn, "logP", "numeric")
        hbd_type_id = get_or_create_property_type(conn, "h_bond_donors", "numeric")
        hba_type_id = get_or_create_property_type(conn, "h_bond_acceptors", "numeric")
        
        # Fix each reference compound
        fixed_count = 0
        
        for chembl_id in reference_ids:
            # Get molecule ID from ChEMBL ID
            cursor.execute(
                "SELECT id FROM molecules WHERE chembl_id = %s",
                (chembl_id,)
            )
            result = cursor.fetchone()
            
            if not result:
                logger.error(f"❌ Reference compound {chembl_id} not found in database")
                continue
            
            molecule_id = result['id']
            
            # Check critical properties
            cursor.execute(
                """
                SELECT pt.name
                FROM molecular_properties mp
                JOIN property_types pt ON mp.property_type_id = pt.id
                WHERE mp.molecule_id = %s AND pt.name IN ('logP', 'h_bond_donors', 'h_bond_acceptors')
                """,
                (molecule_id,)
            )
            existing_properties = [row['name'] for row in cursor.fetchall()]
            
            # Check if all critical properties are present
            critical_properties = ['logP', 'h_bond_donors', 'h_bond_acceptors']
            missing = [prop for prop in critical_properties if prop not in existing_properties]
            
            if missing:
                logger.warning(f"Reference compound {chembl_id} missing properties: {missing}")
                
                # Insert missing properties
                for prop in missing:
                    if prop == 'logP':
                        prop_type_id = logp_type_id
                        value = 0.0
                    elif prop == 'h_bond_donors':
                        prop_type_id = hbd_type_id
                        value = 0
                    elif prop == 'h_bond_acceptors':
                        prop_type_id = hba_type_id
                        value = 0
                    
                    # Check if property already exists
                    cursor.execute(
                        """
                        SELECT id FROM molecular_properties
                        WHERE molecule_id = %s AND property_type_id = %s
                        """,
                        (molecule_id, prop_type_id)
                    )
                    if cursor.fetchone():
                        # Update existing property
                        cursor.execute(
                            """
                            UPDATE molecular_properties
                            SET numeric_value = %s, updated_at = NOW()
                            WHERE molecule_id = %s AND property_type_id = %s
                            """,
                            (value, molecule_id, prop_type_id)
                        )
                        logger.info(f"Updated existing property {prop} for {chembl_id}")
                    else:
                        # Insert new property
                        cursor.execute(
                            """
                            INSERT INTO molecular_properties
                            (molecule_id, property_type_id, numeric_value, created_at, updated_at)
                            VALUES (%s, %s, %s, NOW(), NOW())
                            """,
                            (molecule_id, prop_type_id, value)
                        )
                        logger.info(f"Inserted new property {prop} for {chembl_id}")
                
                conn.commit()
                
                # Verify properties were inserted
                cursor.execute(
                    """
                    SELECT pt.name
                    FROM molecular_properties mp
                    JOIN property_types pt ON mp.property_type_id = pt.id
                    WHERE mp.molecule_id = %s AND pt.name IN ('logP', 'h_bond_donors', 'h_bond_acceptors')
                    """,
                    (molecule_id,)
                )
                updated_properties = [row['name'] for row in cursor.fetchall()]
                
                still_missing = [prop for prop in critical_properties if prop not in updated_properties]
                
                if still_missing:
                    logger.error(f"❌ Failed to fix all properties for {chembl_id}. Still missing: {still_missing}")
                else:
                    logger.info(f"✅ Successfully fixed all properties for {chembl_id}")
                    fixed_count += 1
            else:
                logger.info(f"✅ Reference compound {chembl_id} already has all critical properties")
                fixed_count += 1
        
        # Close database connection
        conn.close()
        
        # Generate report
        logger.info(f"Reference compounds property fixing complete:")
        logger.info(f"  Total reference compounds: {len(reference_ids)}")
        logger.info(f"  Fixed reference compounds: {fixed_count}")
        
        # Return success if all reference compounds are fixed
        return fixed_count == len(reference_ids)
    except Exception as e:
        logger.error(f"❌ Failed to fix reference compounds properties: {str(e)}")
        import traceback
        logger.error(traceback.format_exc())
        return False

def verify_reference_compounds():
    """Verify that reference compounds have all required properties."""
    logger.info("Verifying reference compounds...")
    
    try:
        # Get reference compound IDs
        reference_ids = get_reference_compound_ids()
        logger.info(f"Found {len(reference_ids)} reference compound IDs")
        
        # Connect to database
        conn = get_db_connection()
        cursor = conn.cursor(cursor_factory=RealDictCursor)
        
        # Check each reference compound
        complete_count = 0
        incomplete_count = 0
        missing_properties = {}
        
        for chembl_id in reference_ids:
            # Get molecule ID from ChEMBL ID
            cursor.execute(
                "SELECT id FROM molecules WHERE chembl_id = %s",
                (chembl_id,)
            )
            result = cursor.fetchone()
            
            if not result:
                logger.error(f"❌ Reference compound {chembl_id} not found in database")
                continue
            
            molecule_id = result['id']
            
            # Check critical properties
            cursor.execute(
                """
                SELECT pt.name, mp.numeric_value
                FROM molecular_properties mp
                JOIN property_types pt ON mp.property_type_id = pt.id
                WHERE mp.molecule_id = %s AND pt.name IN ('logP', 'h_bond_donors', 'h_bond_acceptors')
                """,
                (molecule_id,)
            )
            properties = {row['name']: row['numeric_value'] for row in cursor.fetchall()}
            
            logger.info(f"Reference compound {chembl_id} (ID: {molecule_id}) properties: {properties}")
            
            # Check if all critical properties are present
            critical_properties = ['logP', 'h_bond_donors', 'h_bond_acceptors']
            missing = [prop for prop in critical_properties if prop not in properties]
            
            if missing:
                logger.warning(f"Reference compound {chembl_id} missing properties: {missing}")
                missing_properties[chembl_id] = missing
                incomplete_count += 1
            else:
                logger.info(f"✅ Reference compound {chembl_id} has all critical properties")
                complete_count += 1
        
        # Close database connection
        conn.close()
        
        # Generate report
        logger.info(f"Reference compounds verification complete:")
        logger.info(f"  Total reference compounds: {len(reference_ids)}")
        logger.info(f"  Complete reference compounds: {complete_count}")
        logger.info(f"  Incomplete reference compounds: {incomplete_count}")
        
        if incomplete_count > 0:
            logger.warning("Incomplete reference compounds:")
            for chembl_id, missing in missing_properties.items():
                logger.warning(f"  {chembl_id}: missing {', '.join(missing)}")
        
        # Return success if all reference compounds are complete
        return complete_count == len(reference_ids)
    except Exception as e:
        logger.error(f"❌ Failed to verify reference compounds: {str(e)}")
        import traceback
        logger.error(traceback.format_exc())
        return False

def main():
    """Main function to run the fixing process."""
    logger.info("Starting reference compounds property fixing process...")
    
    # Step 1: Verify reference compounds
    logger.info("Verifying reference compounds...")
    if verify_reference_compounds():
        logger.info("✅ All reference compounds already have all required properties")
    else:
        logger.warning("⚠️ Some reference compounds are missing required properties")
        
        # Step 2: Fix reference compounds properties
        logger.info("Fixing reference compounds properties...")
        if fix_reference_compounds_properties():
            logger.info("✅ Successfully fixed all reference compounds properties")
        else:
            logger.warning("⚠️ Failed to fix all reference compounds properties")
        
        # Step 3: Final verification
        logger.info("Performing final verification...")
        if verify_reference_compounds():
            logger.info("✅ All reference compounds now have all required properties")
        else:
            logger.error("❌ Some reference compounds are still missing required properties")
    
    logger.info("Reference compounds property fixing process complete.")
    return 0

if __name__ == "__main__":
    sys.exit(main())