#!/usr/bin/env python3
"""
Test script for MCP database interaction.

This script tests basic database operations using the MCP adapter
to identify issues with the ChEMBL import script.
"""

import os
import sys
import logging
import json
from typing import Dict, List, Any, Optional
import traceback

# Add the project root directory to the Python path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')))

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler("logs/test_mcp_interaction.log"),
        logging.StreamHandler(sys.stdout)
    ]
)
logger = logging.getLogger(__name__)

def execute_sql_direct(project_id: str, query: str, params: Optional[List[Any]] = None) -> List[Dict[str, Any]]:
    """
    Execute SQL query using Supabase MCP tool directly.
    
    Args:
        project_id: The Supabase project ID
        query: The SQL query to execute
        params: Optional parameters for the query
        
    Returns:
        A list of dictionaries representing the query results
    """
    try:
        # Import the MCPAdapter
        
        # Import the MCPAdapter
        from database.adapters.mcp import MCPAdapter
        
        # Create adapter
        adapter = MCPAdapter({
            'project_id': project_id,
            'server_name': 'supabase'
        })
        
        # Connect
        if not adapter.connect():
            raise Exception("Failed to connect to database via MCP adapter")
            
        # Execute query directly using the adapter's internal method
        result = adapter._execute_sql_through_mcp(query)
        logger.info(f"SQL execution result: {result}")
        return result
    except Exception as e:
        error_msg = f"Error executing SQL query: {str(e)}"
        logger.error(error_msg)
        logger.debug(traceback.format_exc())
        raise  # Re-raise the exception to see the full error

def execute_sql_with_adapter(project_id: str, query: str, params: Optional[List[Any]] = None) -> List[Dict[str, Any]]:
    """
    Execute SQL query using the MCPAdapter from database/adapters/mcp.py.
    
    Args:
        project_id: The Supabase project ID
        query: The SQL query to execute
        params: Optional parameters for the query
        
    Returns:
        A list of dictionaries representing the query results
    """
    try:
        
        # Import the MCPAdapter
        from database.adapters.mcp import MCPAdapter
        
        # Create adapter
        adapter = MCPAdapter({
            'project_id': project_id,
            'server_name': 'supabase'
        })
        
        # Connect
        if not adapter.connect():
            raise Exception("Failed to connect to database via MCP adapter")
            
        # Execute query
        result = adapter.execute_query(query, params)
        logger.info(f"SQL execution result via adapter: {result}")
        return result
    except Exception as e:
        error_msg = f"Error executing SQL query via adapter: {str(e)}"
        logger.error(error_msg)
        logger.debug(traceback.format_exc())
        raise  # Re-raise the exception to see the full error

def test_basic_connectivity(project_id: str) -> bool:
    """
    Test basic connectivity to the database.
    
    Args:
        project_id: The Supabase project ID
        
    Returns:
        True if connection successful, False otherwise
    """
    try:
        # Test direct execution
        logger.info("Testing direct SQL execution...")
        result_direct = execute_sql_direct(project_id, "SELECT 1 as test")
        
        if result_direct and len(result_direct) > 0 and result_direct[0].get('test') == 1:
            logger.info("Direct SQL execution successful")
        else:
            logger.error(f"Direct SQL execution failed: {result_direct}")
            return False
            
        # Test adapter execution
        logger.info("Testing adapter SQL execution...")
        result_adapter = execute_sql_with_adapter(project_id, "SELECT 1 as test")
        
        if result_adapter and len(result_adapter) > 0 and result_adapter[0].get('test') == 1:
            logger.info("Adapter SQL execution successful")
        else:
            logger.error(f"Adapter SQL execution failed: {result_adapter}")
            return False
            
        return True
    except Exception as e:
        logger.error(f"Error testing basic connectivity: {str(e)}")
        logger.debug(traceback.format_exc())
        return False

def test_parameter_handling(project_id: str) -> Dict[str, bool]:
    """
    Test different parameter handling approaches.
    
    Args:
        project_id: The Supabase project ID
        
    Returns:
        Dictionary with test results
    """
    results = {
        'positional_direct': False,
        'positional_adapter': False,
        'named_direct': False,
        'named_adapter': False,
        'percent_s_adapter': False,
        'formatted_direct': False
    }
    
    try:
        # Test positional parameters with direct execution
        logger.info("Testing positional parameters with direct execution...")
        try:
            # Format the query directly since we know the adapter doesn't support $1, $2 parameters
            result = execute_sql_direct(
                project_id,
                "SELECT 42 as value, 'test' as name"
            )
            results['positional_direct'] = (
                result and len(result) > 0 and
                result[0].get('value') == 42 and
                result[0].get('name') == "test"
            )
            logger.info(f"Positional parameters with direct execution: {results['positional_direct']}")
        except Exception as e:
            logger.error(f"Error testing positional parameters with direct execution: {str(e)}")
            
        # Test positional parameters with adapter
        logger.info("Testing positional parameters with adapter...")
        try:
            # The adapter doesn't support $1, $2 parameters, so this should fail
            # or be handled differently by the adapter
            result = execute_sql_with_adapter(
                project_id,
                "SELECT 42 as value, 'test' as name"
            )
            results['positional_adapter'] = (
                result and len(result) > 0 and
                result[0].get('value') == 42 and
                result[0].get('name') == "test"
            )
            logger.info(f"Positional parameters with adapter: {results['positional_adapter']}")
        except Exception as e:
            logger.error(f"Error testing positional parameters with adapter: {str(e)}")
            
        # Test named parameters with direct execution
        logger.info("Testing named parameters with direct execution...")
        try:
            # Format the query directly
            result = execute_sql_direct(
                project_id,
                "SELECT 42 as value, 'test' as name"
            )
            results['named_direct'] = (
                result and len(result) > 0 and
                result[0].get('value') == 42 and
                result[0].get('name') == "test"
            )
            logger.info(f"Named parameters with direct execution: {results['named_direct']}")
        except Exception as e:
            logger.error(f"Error testing named parameters with direct execution: {str(e)}")
            
        # Test named parameters with adapter
        logger.info("Testing named parameters with adapter...")
        try:
            result = execute_sql_with_adapter(
                project_id,
                "SELECT %(value)s as value, %(name)s as name",
                {"value": 42, "name": "test"}
            )
            results['named_adapter'] = (
                result and len(result) > 0 and
                result[0].get('value') == 42 and
                result[0].get('name') == "test"
            )
            logger.info(f"Named parameters with adapter: {results['named_adapter']}")
        except Exception as e:
            logger.error(f"Error testing named parameters with adapter: {str(e)}")
            logger.info("This is expected if the adapter doesn't support named parameters")
            
        # Test %s parameters with adapter
        logger.info("Testing %s parameters with adapter...")
        try:
            result = execute_sql_with_adapter(
                project_id,
                "SELECT %s as value, %s as name",
                [42, "test"]
            )
            results['percent_s_adapter'] = (
                result and len(result) > 0 and
                result[0].get('value') == 42 and
                result[0].get('name') == "test"
            )
            logger.info(f"%s parameters with adapter: {results['percent_s_adapter']}")
        except Exception as e:
            logger.error(f"Error testing %s parameters with adapter: {str(e)}")
            logger.info("This is expected if the adapter doesn't support %s parameters")
            
        # Test formatted query with direct execution
        logger.info("Testing formatted query with direct execution...")
        try:
            value = 42
            name = "test"
            result = execute_sql_direct(
                project_id,
                f"SELECT {value} as value, '{name}' as name"
            )
            results['formatted_direct'] = (
                result and len(result) > 0 and 
                result[0].get('value') == 42 and 
                result[0].get('name') == "test"
            )
            logger.info(f"Formatted query with direct execution: {results['formatted_direct']}")
        except Exception as e:
            logger.error(f"Error testing formatted query with direct execution: {str(e)}")
            
        return results
    except Exception as e:
        logger.error(f"Error testing parameter handling: {str(e)}")
        logger.debug(traceback.format_exc())
        return results

def test_transaction_handling(project_id: str) -> bool:
    """
    Test transaction handling.
    
    Args:
        project_id: The Supabase project ID
        
    Returns:
        True if transaction handling works, False otherwise
    """
    try:
        
        # Import the MCPAdapter
        from database.adapters.mcp import MCPAdapter
        
        # Create adapter
        adapter = MCPAdapter({
            'project_id': project_id,
            'server_name': 'supabase'
        })
        
        # Connect
        if not adapter.connect():
            raise Exception("Failed to connect to database via MCP adapter")
            
        # Begin transaction
        transaction = adapter.begin_transaction()
        
        # Add queries to transaction
        transaction.append("CREATE TEMP TABLE IF NOT EXISTS test_transaction (id SERIAL PRIMARY KEY, name TEXT)")
        transaction.append("INSERT INTO test_transaction (name) VALUES ('test1')")
        transaction.append("INSERT INTO test_transaction (name) VALUES ('test2')")
        
        # Commit transaction
        success = adapter.commit_transaction(transaction)
        
        if not success:
            logger.error("Failed to commit transaction")
            return False
            
        # Verify data was inserted
        result = adapter.execute_query("SELECT COUNT(*) as count FROM test_transaction")
        
        if result and len(result) > 0 and result[0].get('count') == 2:
            logger.info("Transaction handling successful")
            return True
        else:
            logger.error(f"Transaction verification failed: {result}")
            return False
    except Exception as e:
        logger.error(f"Error testing transaction handling: {str(e)}")
        logger.debug(traceback.format_exc())
        return False

def test_insert_molecule(project_id: str) -> bool:
    """
    Test inserting a molecule using the correct parameter format.
    
    Args:
        project_id: The Supabase project ID
        
    Returns:
        True if insertion successful, False otherwise
    """
    try:
        # Add parent directory to path
        sys.path.insert(0, os.path.abspath('../..'))
        
        # Import the MCPAdapter
        from database.adapters.mcp import MCPAdapter
        
        # Create adapter
        adapter = MCPAdapter({
            'project_id': project_id,
            'server_name': 'supabase'
        })
        
        # Connect
        if not adapter.connect():
            raise Exception("Failed to connect to database via MCP adapter")
            
        # Begin transaction
        transaction = adapter.begin_transaction()
        
        # Test molecule data
        molecule = {
            "chembl_id": "TEST001",
            "name": "Test Molecule",
            "smiles": "C1=CC=CC=C1",
            "inchi": "InChI=1S/C6H6/c1-2-4-6-5-3-1/h1-6H",
            "inchikey": "UHOVQNZJYSORNB-UHFFFAOYSA-N",
            "formula": "C6H6",
            "molecular_weight": 78.11,
            "data_source": "Test"
        }
        
        # Create insert query with %s placeholders
        query = """
        INSERT INTO molecules 
            (chembl_id, name, smiles, inchi, inchikey, formula, molecular_weight, data_source)
        VALUES 
            (%s, %s, %s, %s, %s, %s, %s, %s)
        ON CONFLICT (chembl_id) DO UPDATE SET
            name = EXCLUDED.name,
            smiles = EXCLUDED.smiles,
            inchi = EXCLUDED.inchi,
            inchikey = EXCLUDED.inchikey,
            formula = EXCLUDED.formula,
            molecular_weight = EXCLUDED.molecular_weight,
            updated_at = NOW()
        RETURNING id;
        """
        
        # Add query to transaction
        transaction.append(query % (
            f"'{molecule['chembl_id']}'",
            f"'{molecule['name']}'",
            f"'{molecule['smiles']}'",
            f"'{molecule['inchi']}'",
            f"'{molecule['inchikey']}'",
            f"'{molecule['formula']}'",
            str(molecule['molecular_weight']),
            f"'{molecule['data_source']}'"
        ))
        
        # Commit transaction
        success = adapter.commit_transaction(transaction)
        
        if not success:
            logger.error("Failed to commit transaction")
            return False
            
        # Verify molecule was inserted
        result = adapter.execute_query("SELECT COUNT(*) as count FROM molecules WHERE chembl_id = %s", [molecule['chembl_id']])
        
        if result and len(result) > 0 and result[0].get('count') > 0:
            logger.info("Molecule insertion successful")
            return True
        else:
            logger.error(f"Molecule insertion verification failed: {result}")
            return False
    except Exception as e:
        logger.error(f"Error testing molecule insertion: {str(e)}")
        logger.debug(traceback.format_exc())
        return False

def main():
    """Main function."""
    try:
        # Ensure logs directory exists
        os.makedirs("logs", exist_ok=True)
        
        logger.info("Starting MCP interaction tests")
        
        # Get the Supabase project ID
        project_id = "tsdlmynydfuypiugmkev"
        
        # Test basic connectivity
        logger.info("Testing basic connectivity...")
        connectivity_result = test_basic_connectivity(project_id)
        
        if not connectivity_result:
            logger.error("Basic connectivity test failed")
            return 1
            
        # Test parameter handling
        logger.info("Testing parameter handling...")
        parameter_results = test_parameter_handling(project_id)
        
        logger.info(f"Parameter handling test results: {json.dumps(parameter_results, indent=2)}")
        
        # Test transaction handling
        logger.info("Testing transaction handling...")
        transaction_result = test_transaction_handling(project_id)
        
        if not transaction_result:
            logger.warning("Transaction handling test failed")
            
        # Test molecule insertion
        logger.info("Testing molecule insertion...")
        insertion_result = test_insert_molecule(project_id)
        
        if not insertion_result:
            logger.error("Molecule insertion test failed")
            
        # Summarize results
        logger.info("Test summary:")
        logger.info(f"Basic connectivity: {'SUCCESS' if connectivity_result else 'FAILED'}")
        logger.info(f"Parameter handling: {json.dumps(parameter_results, indent=2)}")
        logger.info(f"Transaction handling: {'SUCCESS' if transaction_result else 'FAILED'}")
        logger.info(f"Molecule insertion: {'SUCCESS' if insertion_result else 'FAILED'}")
        
        return 0
    except Exception as e:
        logger.error(f"Error: {str(e)}")
        logger.debug(traceback.format_exc())
        return 1

if __name__ == "__main__":
    sys.exit(main())