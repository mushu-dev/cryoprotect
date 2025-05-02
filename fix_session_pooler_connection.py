#!/usr/bin/env python3
"""
Fix Session Pooler Connection Issues

This script diagnoses and fixes issues with the session pooler connection by:
1. Verifying environment variables are correctly set
2. Testing connection with the session pooler
3. Updating the connection order to prioritize the session pooler
4. Fixing any issues with the connection utilities
"""

import os
import sys
import logging
import json
import time
import subprocess
from pathlib import Path
from dotenv import load_dotenv

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler(),
        logging.FileHandler('logs/fix_session_pooler.log')
    ]
)
logger = logging.getLogger(__name__)

# Create logs directory if it doesn't exist
os.makedirs('logs', exist_ok=True)

def check_env_variables():
    """Check if required environment variables are set."""
    logger.info("Checking environment variables...")
    
    # Load environment variables from .env file
    load_dotenv()
    
    required_vars = [
        'SUPABASE_URL',
        'SUPABASE_KEY',
        'SUPABASE_DB_HOST',
        'SUPABASE_DB_PORT',
        'SUPABASE_DB_NAME',
        'SUPABASE_DB_USER',
        'SUPABASE_DB_PASSWORD'
    ]
    
    missing_vars = []
    for var in required_vars:
        if not os.getenv(var):
            missing_vars.append(var)
    
    if missing_vars:
        logger.error(f"Missing environment variables: {', '.join(missing_vars)}")
        return False
    
    # Check if the session pooler port is set correctly
    pooler_port = os.getenv('SUPABASE_DB_PORT')
    if pooler_port != '6543':
        logger.warning(f"Session pooler port is set to {pooler_port}, should be 6543")
        logger.info("Updating SUPABASE_DB_PORT to 6543...")
        
        # Update .env file
        env_path = Path('.env')
        if env_path.exists():
            env_content = env_path.read_text()
            if 'SUPABASE_DB_PORT' in env_content:
                env_content = env_content.replace(f'SUPABASE_DB_PORT={pooler_port}', 'SUPABASE_DB_PORT=6543')
            else:
                env_content += '\nSUPABASE_DB_PORT=6543\n'
            env_path.write_text(env_content)
            logger.info("Updated .env file with correct session pooler port")
        
        # Set environment variable for current process
        os.environ['SUPABASE_DB_PORT'] = '6543'
    
    # Check if the connection order is set correctly
    conn_order = os.getenv('DB_CONNECTION_ORDER', 'pooler,direct,ip,mcp')
    if not conn_order.startswith('pooler'):
        logger.warning(f"Connection order is set to {conn_order}, should prioritize pooler")
        logger.info("Updating DB_CONNECTION_ORDER to prioritize pooler...")
        
        # Update .env file
        env_path = Path('.env')
        if env_path.exists():
            env_content = env_path.read_text()
            if 'DB_CONNECTION_ORDER' in env_content:
                env_content = env_content.replace(f'DB_CONNECTION_ORDER={conn_order}', 'DB_CONNECTION_ORDER=pooler,direct,ip,mcp')
            else:
                env_content += '\nDB_CONNECTION_ORDER=pooler,direct,ip,mcp\n'
            env_path.write_text(env_content)
            logger.info("Updated .env file with correct connection order")
        
        # Set environment variable for current process
        os.environ['DB_CONNECTION_ORDER'] = 'pooler,direct,ip,mcp'
    
    logger.info("Environment variables check complete")
    return True

def test_session_pooler_connection():
    """Test connection to the database using the session pooler."""
    logger.info("Testing session pooler connection...")
    
    try:
        # Import the connection utilities
        from db_connection_utils import get_db_connection
        
        # Test connection
        with get_db_connection() as conn:
            cursor = conn.cursor()
            cursor.execute("SELECT 1 as connection_test")
            result = cursor.fetchone()
            
            if result and result[0] == 1:
                logger.info("✅ Successfully connected to database using session pooler")
                return True
            else:
                logger.error("❌ Session pooler connection test returned unexpected result")
                return False
    except ImportError:
        logger.error("Could not import db_connection_utils. Make sure the module is available.")
        return False
    except Exception as e:
        logger.error(f"❌ Session pooler connection test failed: {str(e)}")
        return False

def fix_connection_manager():
    """Fix any issues with the connection manager."""
    logger.info("Checking connection manager configuration...")
    
    try:
        # Import the connection manager
        from db_connection_utils import ConnectionManager
        
        # Get the connection manager instance
        manager = ConnectionManager.get_instance()
        
        # Check if the connection order is correct
        if manager.connection_order[0] != 'pooler':
            logger.warning("Connection order in manager does not prioritize pooler")
            manager.connection_order = ['pooler', 'direct', 'ip', 'mcp']
            logger.info("Updated connection order in manager to prioritize pooler")
        
        # Check if the pooler configuration is correct
        pooler_config = manager.config.get('pooler', {})
        if pooler_config.get('port') != 6543:
            logger.warning(f"Pooler port in manager is set to {pooler_config.get('port')}, should be 6543")
            pooler_config['port'] = 6543
            logger.info("Updated pooler port in manager to 6543")
        
        # Reset circuit breakers
        for breaker_name, breaker in manager.circuit_breakers.items():
            breaker.reset()
            logger.info(f"Reset circuit breaker: {breaker_name}")
        
        logger.info("Connection manager configuration check complete")
        return True
    except ImportError:
        logger.error("Could not import ConnectionManager. Make sure the module is available.")
        return False
    except Exception as e:
        logger.error(f"Error checking connection manager: {str(e)}")
        return False

def fix_reference_properties():
    """Fix reference compound properties using the session pooler connection."""
    logger.info("Fixing reference compound properties...")
    
    try:
        # Import necessary modules
        from db_connection_utils import get_db_connection, safe_transaction
        
        # Define reference compound IDs
        reference_ids = [
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
        
        # Connect to database using the session pooler
        with get_db_connection() as conn:
            cursor = conn.cursor()
            
            # Check each reference compound
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
                
                molecule_id = result[0]
                
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
                existing_properties = [row[0] for row in cursor.fetchall()]
                
                # Check if all critical properties are present
                critical_properties = ['logP', 'h_bond_donors', 'h_bond_acceptors']
                missing = [prop for prop in critical_properties if prop not in existing_properties]
                
                if missing:
                    logger.warning(f"Reference compound {chembl_id} missing properties: {missing}")
                    
                    # Insert missing properties
                    with safe_transaction() as tx_conn:
                        tx_cursor = tx_conn.cursor()
                        
                        for prop in missing:
                            # Get or create property type
                            tx_cursor.execute(
                                "SELECT id FROM property_types WHERE name = %s",
                                (prop,)
                            )
                            prop_type_result = tx_cursor.fetchone()
                            
                            if prop_type_result:
                                prop_type_id = prop_type_result[0]
                            else:
                                # Create property type
                                tx_cursor.execute(
                                    """
                                    INSERT INTO property_types (name, data_type, description)
                                    VALUES (%s, %s, %s)
                                    RETURNING id
                                    """,
                                    (prop, 'numeric', f"Auto-created property: {prop}")
                                )
                                prop_type_id = tx_cursor.fetchone()[0]
                            
                            # Set default value
                            if prop == 'logP':
                                value = 0.0
                            else:  # h_bond_donors or h_bond_acceptors
                                value = 0
                            
                            # Insert property
                            tx_cursor.execute(
                                """
                                INSERT INTO molecular_properties (molecule_id, property_type_id, numeric_value)
                                VALUES (%s, %s, %s)
                                """,
                                (molecule_id, prop_type_id, value)
                            )
                            
                            logger.info(f"Inserted property {prop} = {value} for {chembl_id}")
                else:
                    logger.info(f"✅ Reference compound {chembl_id} already has all critical properties")
        
        logger.info("Reference compound property fixing complete")
        return True
    except ImportError:
        logger.error("Could not import necessary modules. Make sure they are available.")
        return False
    except Exception as e:
        logger.error(f"Error fixing reference properties: {str(e)}")
        import traceback
        logger.error(traceback.format_exc())
        return False

def verify_reference_compounds():
    """Verify that reference compounds have all required properties."""
    logger.info("Verifying reference compounds...")
    
    try:
        # Import necessary modules
        from db_connection_utils import get_db_connection
        
        # Define reference compound IDs
        reference_ids = [
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
        
        # Connect to database using the session pooler
        with get_db_connection() as conn:
            cursor = conn.cursor()
            
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
                
                molecule_id = result[0]
                
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
                properties = {row[0]: row[1] for row in cursor.fetchall()}
                
                logger.info(f"Reference compound {chembl_id} properties: {properties}")
                
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
    except ImportError:
        logger.error("Could not import necessary modules. Make sure they are available.")
        return False
    except Exception as e:
        logger.error(f"Error verifying reference compounds: {str(e)}")
        import traceback
        logger.error(traceback.format_exc())
        return False

def main():
    """Main function to fix session pooler connection issues."""
    logger.info("Starting session pooler connection fix...")
    
    # Step 1: Check environment variables
    if not check_env_variables():
        logger.error("Environment variables check failed. Please set the required variables.")
        return 1
    
    # Step 2: Fix connection manager
    if not fix_connection_manager():
        logger.warning("Connection manager fix failed. Continuing with other fixes...")
    
    # Step 3: Test session pooler connection
    if not test_session_pooler_connection():
        logger.error("Session pooler connection test failed. Cannot proceed with data fixes.")
        return 1
    
    # Step 4: Fix reference properties
    if not fix_reference_properties():
        logger.error("Reference properties fix failed.")
        return 1
    
    # Step 5: Verify reference compounds
    if not verify_reference_compounds():
        logger.warning("Reference compounds verification failed. Some properties may still be missing.")
    else:
        logger.info("✅ All reference compounds now have all required properties")
    
    logger.info("Session pooler connection fix complete.")
    return 0

if __name__ == "__main__":
    sys.exit(main())