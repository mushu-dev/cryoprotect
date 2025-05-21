#!/usr/bin/env python3
"""
Script to improve relationships between molecules and their properties.
This script applies the SQL migration to enhance property relationships
and then tests the new functionality.
"""

import os
import sys
import argparse
import json
import logging
from datetime import datetime
import db_utils
from psycopg2.extras import RealDictCursor

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

def execute_migration():
    """Execute the migration to improve property relationships."""
    sql_file = os.path.join(
        os.path.dirname(os.path.realpath(__file__)), 
        'migrations', 
        '027_improve_property_relationships.sql'
    )
    
    try:
        if db_utils.execute_sql_file(sql_file):
            logger.info("Successfully executed migration to improve property relationships")
            return True
        else:
            logger.error("Failed to execute migration")
            return False
    except Exception as e:
        logger.error(f"Error executing migration: {e}")
        return False

def test_molecule_property_functions():
    """Test the new molecule property functions."""
    tests = []
    
    try:
        # Get a molecule ID to test with
        query = "SELECT id, name FROM molecules LIMIT 1"
        result = db_utils.execute_query(query, cursor_factory=RealDictCursor)
        
        if not result:
            logger.error("No molecules found to test with")
            return []
        
        molecule_id = result[0]['id']
        molecule_name = result[0]['name']
        logger.info(f"Testing with molecule: {molecule_name} (ID: {molecule_id})")
        
        # Test setting a property
        test_property = f"test_property_{datetime.now().strftime('%Y%m%d%H%M%S')}"
        test_value = "42.5"
        
        query = f"""
        SELECT set_molecule_property('{molecule_id}', '{test_property}', '{test_value}', 'numeric');
        """
        db_utils.execute_query(query, fetch=False)
        logger.info(f"Set property {test_property} = {test_value} for molecule {molecule_name}")
        
        # Test getting the property
        query = f"""
        SELECT get_molecule_property('{molecule_id}', '{test_property}');
        """
        result = db_utils.execute_query(query, cursor_factory=RealDictCursor)[0]
        retrieved_value = result['get_molecule_property']
        logger.info(f"Retrieved property {test_property} = {retrieved_value}")
        
        tests.append({
            'test': 'set_and_get_property',
            'molecule_id': str(molecule_id),
            'property_name': test_property,
            'value_set': test_value,
            'value_retrieved': retrieved_value,
            'success': retrieved_value == test_value
        })
        
        # Test searching for molecules by property
        query = f"""
        SELECT * FROM search_molecules_by_property('{test_property}', '{test_value}');
        """
        result = db_utils.execute_query(query, cursor_factory=RealDictCursor)
        logger.info(f"Found {len(result)} molecules with {test_property} = {test_value}")
        
        tests.append({
            'test': 'search_by_property',
            'property_name': test_property,
            'property_value': test_value,
            'molecules_found': len(result),
            'success': len(result) > 0
        })
        
        # Test the molecule_with_properties view
        query = f"""
        SELECT * FROM molecule_with_properties WHERE id = '{molecule_id}';
        """
        result = db_utils.execute_query(query, cursor_factory=RealDictCursor)
        
        if result and len(result) > 0:
            properties = result[0]['properties']
            logger.info(f"Retrieved {len(properties) if properties else 0} properties for molecule {molecule_name} from the view")
            
            tests.append({
                'test': 'molecule_with_properties_view',
                'molecule_id': str(molecule_id),
                'properties_retrieved': properties is not None,
                'property_count': len(properties) if properties else 0,
                'success': properties is not None and test_property in properties
            })
        else:
            logger.error(f"Failed to retrieve molecule {molecule_name} from the view")
            tests.append({
                'test': 'molecule_with_properties_view',
                'molecule_id': str(molecule_id),
                'success': False,
                'error': 'No results returned from view'
            })
        
        return tests
    except Exception as e:
        logger.error(f"Error testing molecule property functions: {e}")
        return []

def save_report(tests, filename_prefix="property_relationship_improvement_report"):
    """Save the test report to a file."""
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    filename = f"{filename_prefix}_{timestamp}.json"
    
    all_success = all(test['success'] for test in tests)
    
    report = {
        'timestamp': datetime.now().isoformat(),
        'tests': tests,
        'all_tests_passed': all_success,
        'total_tests': len(tests),
        'successful_tests': sum(1 for test in tests if test['success'])
    }
    
    try:
        with open(filename, 'w') as f:
            json.dump(report, f, indent=2)
        logger.info(f"Report saved to {filename}")
        return filename
    except Exception as e:
        logger.error(f"Error saving report: {e}")
        return None

def main():
    """Main function to run the script."""
    parser = argparse.ArgumentParser(description="Improve relationships between molecules and their properties.")
    parser.add_argument("--test-only", action="store_true", help="Only test the existing functions without applying the migration")
    args = parser.parse_args()
    
    # Check database connection
    if not db_utils.test_connection():
        logger.error("Database connection failed. Exiting.")
        sys.exit(1)
    
    if not args.test_only:
        # Execute the migration
        logger.info("Executing migration to improve property relationships...")
        success = execute_migration()
        
        if not success:
            logger.error("Migration failed. Exiting.")
            sys.exit(1)
    
    # Test the molecule property functions
    logger.info("Testing molecule property functions...")
    tests = test_molecule_property_functions()
    
    if not tests:
        logger.error("Failed to run tests. Exiting.")
        sys.exit(1)
    
    # Print test results
    all_passed = all(test['success'] for test in tests)
    logger.info(f"Test results: {sum(1 for test in tests if test['success'])}/{len(tests)} tests passed")
    
    if all_passed:
        logger.info("All tests passed!")
    else:
        logger.warning("Some tests failed. Check the report for details.")
    
    # Save report
    report_file = save_report(tests)
    if report_file:
        logger.info(f"Test report saved to {report_file}")
    
    if not all_passed and not args.test_only:
        logger.warning("Migration was applied but some tests failed.")
    
    sys.exit(0 if all_passed else 1)

if __name__ == "__main__":
    main()