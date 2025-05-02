"""
Example usage of PostgreSQL direct connection and SQL executor utilities.

This script demonstrates how to use the PostgreSQL direct connection and SQL executor
utilities for common database operations.
"""

import os
import sys
import logging
from dotenv import load_dotenv

# Add the parent directory to the path
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# Import the modules
from postgres_direct import PostgresDirectConnection
import sql_executor

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger('postgres_direct_example')

def setup_example_table():
    """Create an example table for testing."""
    # Check if the table exists
    check_query = """
    SELECT EXISTS (
        SELECT FROM information_schema.tables 
        WHERE table_schema = 'public' 
        AND table_name = 'example_compounds'
    )
    """
    
    table_exists = sql_executor.execute_query(check_query, fetch_one=True)
    
    if table_exists and table_exists.get('exists', False):
        logger.info("Example table already exists")
        return
    
    # Create the table
    create_table_query = """
    CREATE TABLE example_compounds (
        id SERIAL PRIMARY KEY,
        name VARCHAR(255) NOT NULL,
        smiles TEXT,
        inchi_key VARCHAR(27),
        molecular_weight NUMERIC(10, 4),
        logp NUMERIC(6, 3),
        created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
    )
    """
    
    sql_executor.execute_query(create_table_query)
    logger.info("Created example_compounds table")

def insert_example_data():
    """Insert example data into the table."""
    # Example compound data
    compounds = [
        {
            'name': 'Glycerol',
            'smiles': 'C(C(CO)O)O',
            'inchi_key': 'PEDCQBHIVMGVHV-UHFFFAOYSA-N',
            'molecular_weight': 92.0938,
            'logp': -1.76
        },
        {
            'name': 'Dimethyl sulfoxide',
            'smiles': 'CS(=O)C',
            'inchi_key': 'IAZDPXIOMUYVGZ-UHFFFAOYSA-N',
            'molecular_weight': 78.1334,
            'logp': -1.35
        },
        {
            'name': 'Ethylene glycol',
            'smiles': 'C(CO)O',
            'inchi_key': 'LYCAIKOWRPUZTN-UHFFFAOYSA-N',
            'molecular_weight': 62.0678,
            'logp': -1.69
        }
    ]
    
    # Insert the compounds
    inserted_ids = sql_executor.bulk_insert(
        'example_compounds',
        compounds,
        return_ids=True
    )
    
    logger.info(f"Inserted {len(inserted_ids)} compounds with IDs: {inserted_ids}")
    return inserted_ids

def query_example_data():
    """Query and display the example data."""
    # Query all compounds
    query = "SELECT * FROM example_compounds ORDER BY molecular_weight"
    compounds = sql_executor.execute_query(query)
    
    logger.info(f"Found {len(compounds)} compounds:")
    for compound in compounds:
        logger.info(f"  {compound['name']} (MW: {compound['molecular_weight']}, LogP: {compound['logp']})")
    
    return compounds

def update_example_data(compound_id, new_logp):
    """Update a compound's logP value."""
    # Update the compound
    update_query = "UPDATE example_compounds SET logp = %s WHERE id = %s"
    sql_executor.execute_query(update_query, [new_logp, compound_id])
    
    # Verify the update
    verify_query = "SELECT name, logp FROM example_compounds WHERE id = %s"
    result = sql_executor.execute_query(verify_query, [compound_id], fetch_one=True)
    
    logger.info(f"Updated {result['name']} logP to {result['logp']}")

def batch_update_example_data():
    """Demonstrate batch update operations."""
    # Get all compounds
    compounds = sql_executor.execute_query("SELECT id, name FROM example_compounds")
    
    # Prepare batch update parameters
    update_query = "UPDATE example_compounds SET name = %s || ' (Updated)' WHERE id = %s"
    params_list = [[compound['name'], compound['id']] for compound in compounds]
    
    # Execute batch update
    sql_executor.execute_batch(update_query, params_list)
    
    logger.info(f"Batch updated {len(compounds)} compounds")

def demonstrate_transaction():
    """Demonstrate transaction handling."""
    @sql_executor.with_transaction
    def transaction_example(conn):
        # Create a cursor
        with conn.cursor() as cursor:
            # Execute multiple operations in a transaction
            cursor.execute("UPDATE example_compounds SET molecular_weight = molecular_weight * 1.01 WHERE id = 1")
            cursor.execute("UPDATE example_compounds SET molecular_weight = molecular_weight * 1.01 WHERE id = 2")
            cursor.execute("UPDATE example_compounds SET molecular_weight = molecular_weight * 1.01 WHERE id = 3")
    
    # Execute the transaction
    transaction_example()
    logger.info("Transaction completed successfully")

def demonstrate_retry():
    """Demonstrate retry logic."""
    # Define a function that fails on first attempt
    attempt_count = 0
    
    @sql_executor.with_retry(max_retries=3, retry_delay=1)
    def retry_example():
        nonlocal attempt_count
        attempt_count += 1
        
        if attempt_count < 2:
            logger.info("Simulating a failure...")
            raise Exception("Simulated failure")
        
        logger.info("Operation succeeded!")
        return "Success"
    
    # Execute the function with retry logic
    result = retry_example()
    logger.info(f"Retry result: {result} (took {attempt_count} attempts)")

def demonstrate_batch_processing():
    """Demonstrate batch processing."""
    # Generate a large list of items
    items = list(range(1, 101))
    
    # Process in batches
    def process_batch(batch):
        logger.info(f"Processing batch with {len(batch)} items")
        return sum(batch)
    
    results = sql_executor.process_in_batches(
        items,
        batch_size=25,
        process_func=process_batch
    )
    
    logger.info(f"Batch processing results: {results}")
    logger.info(f"Total sum: {sum(results)}")

def cleanup_example():
    """Clean up the example by dropping the table."""
    response = input("Do you want to drop the example_compounds table? (y/n): ")
    
    if response.lower() == 'y':
        sql_executor.execute_query("DROP TABLE example_compounds")
        logger.info("Dropped example_compounds table")
    else:
        logger.info("Keeping example_compounds table")

def main():
    """Main function to run the example."""
    # Load environment variables
    load_dotenv()
    
    # Check if required environment variables are set
    required_vars = ['SUPABASE_DB_HOST', 'SUPABASE_DB_PASSWORD']
    missing_vars = [var for var in required_vars if not os.getenv(var)]
    
    if missing_vars:
        logger.error(f"Missing required environment variables: {', '.join(missing_vars)}")
        logger.error("Please set these variables in a .env file or environment")
        return
    
    try:
        # Setup
        logger.info("Setting up example...")
        setup_example_table()
        
        # Insert data
        logger.info("\nInserting example data...")
        compound_ids = insert_example_data()
        
        # Query data
        logger.info("\nQuerying example data...")
        query_example_data()
        
        # Update data
        if compound_ids:
            logger.info("\nUpdating example data...")
            update_example_data(compound_ids[0], -1.80)
        
        # Batch update
        logger.info("\nDemonstrating batch update...")
        batch_update_example_data()
        
        # Transaction
        logger.info("\nDemonstrating transaction...")
        demonstrate_transaction()
        
        # Retry logic
        logger.info("\nDemonstrating retry logic...")
        demonstrate_retry()
        
        # Batch processing
        logger.info("\nDemonstrating batch processing...")
        demonstrate_batch_processing()
        
        # Display connection metrics
        logger.info("\nConnection metrics:")
        metrics = sql_executor.get_connection_metrics()
        for key, value in metrics.items():
            if isinstance(value, float):
                logger.info(f"  {key}: {value:.4f}")
            else:
                logger.info(f"  {key}: {value}")
        
        # Cleanup
        logger.info("\nExample completed successfully!")
        cleanup_example()
        
    except Exception as e:
        logger.error(f"Error in example: {str(e)}", exc_info=True)
    finally:
        # Close connections
        sql_executor.close_connections()
        logger.info("Connections closed")

if __name__ == "__main__":
    main()