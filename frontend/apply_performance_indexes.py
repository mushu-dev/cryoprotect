#!/usr/bin/env python3
"""
Script to apply performance-optimizing indexes to the database.
This script executes the 033_add_performance_indexes.sql migration.
"""

import os
import sys
import logging
import json
from datetime import datetime
from dotenv import load_dotenv
import db_utils

# Load environment variables
load_dotenv()

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

def execute_migration():
    """Execute the migration to add performance indexes."""
    sql_file = os.path.join(
        os.path.dirname(os.path.realpath(__file__)), 
        'migrations', 
        '033_add_performance_indexes.sql'
    )
    
    try:
        if db_utils.execute_sql_file(sql_file):
            logger.info("Successfully executed migration to add performance indexes")
            return True
        else:
            logger.error("Failed to execute migration")
            return False
    except Exception as e:
        logger.error(f"Error executing migration: {e}")
        return False

def verify_indexes():
    """Verify that indexes were created successfully."""
    try:
        query = """
        SELECT 
            tablename, 
            indexname, 
            indexdef
        FROM 
            pg_indexes 
        WHERE 
            schemaname = 'public' 
            AND indexname LIKE 'idx_%'
        ORDER BY 
            tablename, 
            indexname;
        """
        
        indexes = db_utils.execute_query(query, cursor_factory=db_utils.RealDictCursor)
        
        # Group by table
        by_table = {}
        for idx in indexes:
            table = idx['tablename']
            if table not in by_table:
                by_table[table] = []
            
            by_table[table].append({
                'name': idx['indexname'],
                'definition': idx['indexdef']
            })
        
        return by_table
    except Exception as e:
        logger.error(f"Error verifying indexes: {e}")
        return None

def save_index_report(indexes):
    """Save the index creation report to a file."""
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    filename = f"performance_indexes_report_{timestamp}.json"
    
    # Calculate some statistics
    total_indexes = sum(len(table_indexes) for table_indexes in indexes.values())
    tables_with_indexes = len(indexes)
    
    report = {
        'timestamp': datetime.now().isoformat(),
        'total_indexes': total_indexes,
        'tables_with_indexes': tables_with_indexes,
        'indexes_by_table': indexes
    }
    
    try:
        with open(filename, 'w') as f:
            json.dump(report, f, indent=2)
        logger.info(f"Report saved to {filename}")
        return filename
    except Exception as e:
        logger.error(f"Error saving report: {e}")
        return None

def simulate_query_improvement():
    """Simulate the improvement in query performance from the added indexes."""
    # Queries to test (similar to what we use in analyze_query_performance.py)
    test_queries = [
        {
            'name': 'Molecule lookup by name',
            'query': 'EXPLAIN ANALYZE SELECT * FROM molecules WHERE name ILIKE \'%glycerol%\';',
            'expected_improvement': 'Uses index idx_molecules_name_trgm for pattern matching'
        },
        {
            'name': 'Property lookup by molecule',
            'query': 'EXPLAIN ANALYZE SELECT * FROM molecular_properties WHERE molecule_id = uuid_generate_v4();',
            'expected_improvement': 'Uses index idx_molecular_properties_molecule_id'
        },
        {
            'name': 'Property lookup by name',
            'query': 'EXPLAIN ANALYZE SELECT * FROM molecular_properties WHERE property_name = \'molecular_weight\';',
            'expected_improvement': 'Uses index idx_molecular_properties_property_name'
        },
        {
            'name': 'Molecule search by formula',
            'query': 'EXPLAIN ANALYZE SELECT * FROM molecules WHERE formula = \'C3H8O3\';',
            'expected_improvement': 'Uses index idx_molecules_formula'
        },
        {
            'name': 'Join molecules and properties with filtering',
            'query': 'EXPLAIN ANALYZE SELECT m.*, mp.numeric_value FROM molecules m JOIN molecular_properties mp ON m.id = mp.molecule_id WHERE mp.property_name = \'melting_point\' AND m.is_public = TRUE ORDER BY mp.numeric_value DESC LIMIT 10;',
            'expected_improvement': 'Uses indexes on molecule_id, property_name, is_public, and numeric_value'
        }
    ]
    
    # Results for each query
    results = []
    
    # Run each test query
    for test in test_queries:
        try:
            query_plan = db_utils.execute_query(test['query'], cursor_factory=db_utils.RealDictCursor)
            
            # Extract execution time and index usage
            execution_time = None
            used_indexes = []
            
            for row in query_plan:
                for key, value in row.items():
                    if isinstance(value, str):
                        # Look for execution time
                        if "execution time" in value.lower():
                            import re
                            match = re.search(r"execution time: (\d+\.\d+)", value.lower())
                            if match:
                                execution_time = float(match.group(1))
                        
                        # Look for index usage
                        if "index scan" in value.lower() or "bitmap index scan" in value.lower():
                            for idx_pattern in ["idx_", "pkey"]:
                                if idx_pattern in value.lower():
                                    import re
                                    match = re.search(rf"on (\w+) .*?({idx_pattern}\w+)", value.lower())
                                    if match:
                                        used_indexes.append(match.group(2))
            
            results.append({
                'name': test['name'],
                'execution_time_ms': execution_time,
                'used_indexes': used_indexes,
                'expected_improvement': test['expected_improvement'],
                'query': test['query']
            })
        except Exception as e:
            logger.warning(f"Error executing test query '{test['name']}': {e}")
            results.append({
                'name': test['name'],
                'error': str(e),
                'expected_improvement': test['expected_improvement'],
                'query': test['query']
            })
    
    return results

def main():
    """Main function."""
    import argparse
    
    parser = argparse.ArgumentParser(description="Apply performance-optimizing indexes to the database.")
    parser.add_argument("--dry-run", action="store_true", help="Only verify existing indexes without creating new ones")
    parser.add_argument("--simulate", action="store_true", help="Simulate query performance improvements")
    args = parser.parse_args()
    
    # Check database connection
    if not db_utils.test_connection():
        logger.error("Database connection failed. Exiting.")
        sys.exit(1)
    
    if args.dry_run:
        logger.info("Running in dry-run mode")
        indexes = verify_indexes()
        
        if not indexes:
            logger.error("Failed to verify indexes")
            sys.exit(1)
        
        report_file = save_index_report(indexes)
        
        if report_file:
            print(f"\nCurrent indexes report saved to {report_file}")
            
            # Print summary
            total_indexes = sum(len(table_indexes) for table_indexes in indexes.values())
            tables_with_indexes = len(indexes)
            
            print(f"\nDatabase has {total_indexes} custom indexes across {tables_with_indexes} tables:")
            for table, table_indexes in indexes.items():
                print(f"- {table}: {len(table_indexes)} indexes")
        
        sys.exit(0)
    
    # Execute the migration
    logger.info("Executing migration to add performance indexes...")
    success = execute_migration()
    
    if not success:
        logger.error("Failed to execute migration. Exiting.")
        sys.exit(1)
    
    # Verify that indexes were created
    logger.info("Verifying indexes...")
    indexes = verify_indexes()
    
    if not indexes:
        logger.error("Failed to verify indexes")
        sys.exit(1)
    
    # Save index report
    report_file = save_index_report(indexes)
    
    # Simulate query improvements if requested
    if args.simulate:
        logger.info("Simulating query improvements...")
        query_results = simulate_query_improvement()
        
        # Save query results to report
        if report_file and query_results:
            # Update the existing report
            try:
                with open(report_file, 'r') as f:
                    report = json.load(f)
                
                report['query_performance'] = query_results
                
                with open(report_file, 'w') as f:
                    json.dump(report, f, indent=2)
                
                logger.info(f"Updated report with query performance results")
            except Exception as e:
                logger.error(f"Error updating report with query results: {e}")
    
    # Print summary
    if report_file:
        print(f"\nIndex creation report saved to {report_file}")
        
        # Print summary
        total_indexes = sum(len(table_indexes) for table_indexes in indexes.values())
        tables_with_indexes = len(indexes)
        
        print(f"\nCreated or verified {total_indexes} custom indexes across {tables_with_indexes} tables:")
        for table, table_indexes in sorted(indexes.items()):
            print(f"- {table}: {len(table_indexes)} indexes")
    
    if args.simulate and query_results:
        print("\nQuery Performance Simulation:")
        for result in query_results:
            print(f"\n- {result['name']}:")
            if 'error' in result:
                print(f"  Error: {result['error']}")
            else:
                print(f"  Execution time: {result.get('execution_time_ms', 'N/A')} ms")
                print(f"  Used indexes: {', '.join(result.get('used_indexes', ['None']))}")
                print(f"  Expected improvement: {result['expected_improvement']}")
    
    print("\nIndex creation completed successfully!")

if __name__ == "__main__":
    main()