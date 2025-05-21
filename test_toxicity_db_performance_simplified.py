#!/usr/bin/env python3
"""
CryoProtect - Toxicity Database Performance Test (Simplified)

A simplified version of the database performance test that works with the
existing database utilities.

Usage:
    python test_toxicity_db_performance_simplified.py [--detailed] [--export FILENAME]
"""

import os
import sys
import json
import time
import argparse
import logging
import statistics
from datetime import datetime
from typing import Dict, List, Any, Tuple

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler(),
        logging.FileHandler(f'toxicity_db_performance_{datetime.now().strftime("%Y%m%d_%H%M%S")}.log')
    ]
)
logger = logging.getLogger(__name__)

# Import database utilities
try:
    # Add your project's root directory to the Python path
    sys.path.append(os.path.dirname(os.path.abspath(__file__)))
    from database.utils import execute_sql
    logger.info("Successfully imported database utilities")
except ImportError as e:
    logger.error(f"Failed to import database utilities: {str(e)}")
    sys.exit(1)

# Function to execute SQL with timing
def execute_sql_timed(query, params=None):
    """Execute SQL and return the result and execution time."""
    start_time = time.time()
    result = execute_sql(query, params)
    execution_time = time.time() - start_time
    return result, execution_time

class DBPerformanceTest:
    """Database performance test for toxicity optimization."""
    
    def __init__(self, detailed=False):
        """Initialize the test suite."""
        self.detailed = detailed
        self.results = {
            'query_times': [],
            'materialized_view_refresh': None
        }
        self.test_molecule_ids = self._get_test_molecule_ids()
    
    def _get_test_molecule_ids(self, limit=5):
        """Get sample molecule IDs for testing."""
        try:
            query = f"SELECT id FROM molecules LIMIT {limit}"
            result = execute_sql(query)
            molecule_ids = [row['id'] for row in result]
            logger.info(f"Retrieved {len(molecule_ids)} test molecule IDs")
            return molecule_ids
        except Exception as e:
            logger.error(f"Failed to get test molecule IDs: {str(e)}")
            return []
    
    def run_all_tests(self):
        """Run all performance tests."""
        logger.info("Starting database performance tests")
        
        # Run query performance tests
        self._test_query_performance()
        
        # Test materialized view refresh
        self._test_materialized_view_refresh()
        
        # Print summary
        self._print_summary()
        
        logger.info("Database performance tests completed")
        return self.results
    
    def _test_query_performance(self):
        """Test query performance comparison between original and optimized schemas."""
        if not self.test_molecule_ids:
            logger.warning("No test molecule IDs available, skipping query performance test")
            return
        
        logger.info("Testing query performance")
        
        # Define test queries
        test_queries = [
            {
                'name': 'Basic Toxicity Data',
                'original': "SELECT * FROM toxicity_data WHERE molecule_id = %s",
                'optimized': "SELECT * FROM toxicity_summary WHERE molecule_id = %s"
            },
            {
                'name': 'Tox21 Activity',
                'original': """
                    SELECT td.*, ta.assay_name, ta.assay_target, ta.toxicological_endpoint 
                    FROM toxicity_data td
                    JOIN toxicity_assay ta ON td.assay_id = ta.id
                    WHERE td.molecule_id = %s AND td.source = 'Tox21'
                """,
                'optimized': "SELECT * FROM tox21_activity_summary WHERE molecule_id = %s"
            },
            {
                'name': 'LD50 Values',
                'original': """
                    SELECT td.* 
                    FROM toxicity_data td
                    WHERE td.molecule_id = %s AND td.toxicity_type = 'LD50'
                    ORDER BY td.value ASC
                """,
                'optimized': "SELECT * FROM ld50_summary WHERE molecule_id = %s ORDER BY ld50_value ASC"
            },
            {
                'name': 'Hazard Classifications',
                'original': """
                    SELECT tc.* 
                    FROM toxicity_classification tc
                    WHERE tc.molecule_id = %s
                """,
                'optimized': "SELECT * FROM hazard_classification_summary WHERE molecule_id = %s"
            },
            {
                'name': 'Toxicity Score Calculation',
                'original': """
                    WITH tox21_count AS (
                        SELECT COUNT(*) AS count
                        FROM toxicity_data
                        WHERE molecule_id = %s AND source = 'Tox21' AND hit_call = true
                    ),
                    ld50_min AS (
                        SELECT MIN(value) AS min_value
                        FROM toxicity_data
                        WHERE molecule_id = %s AND toxicity_type = 'LD50'
                    ),
                    classification_count AS (
                        SELECT COUNT(*) AS count
                        FROM toxicity_classification
                        WHERE molecule_id = %s AND classification_system = 'GHS'
                    )
                    SELECT 
                        COALESCE(tox21_count.count, 0) * 0.2 +
                        CASE 
                            WHEN ld50_min.min_value < 50 THEN 5
                            WHEN ld50_min.min_value < 300 THEN 3
                            WHEN ld50_min.min_value < 2000 THEN 1
                            ELSE 0
                        END +
                        COALESCE(classification_count.count, 0) * 0.5 AS toxicity_score
                    FROM tox21_count, ld50_min, classification_count
                """,
                'optimized': "SELECT calculate_toxicity_score(%s) AS toxicity_score"
            }
        ]
        
        # Run each query for each test molecule
        for query_set in test_queries:
            query_name = query_set['name']
            original_times = []
            optimized_times = []
            
            logger.info(f"Testing {query_name} queries")
            
            for molecule_id in self.test_molecule_ids:
                # Test original query
                try:
                    if query_set['name'] == 'Toxicity Score Calculation':
                        # Special case for toxicity score which needs multiple param references
                        _, original_time = execute_sql_timed(
                            query_set['original'], 
                            [molecule_id, molecule_id, molecule_id]
                        )
                    else:
                        _, original_time = execute_sql_timed(
                            query_set['original'], 
                            [molecule_id]
                        )
                    original_times.append(original_time)
                    
                    if self.detailed:
                        logger.info(f"Original {query_name} query for {molecule_id}: {original_time:.6f}s")
                except Exception as e:
                    logger.warning(f"Error in original {query_name} query: {str(e)}")
                
                # Test optimized query
                try:
                    _, optimized_time = execute_sql_timed(
                        query_set['optimized'], 
                        [molecule_id]
                    )
                    optimized_times.append(optimized_time)
                    
                    if self.detailed:
                        logger.info(f"Optimized {query_name} query for {molecule_id}: {optimized_time:.6f}s")
                except Exception as e:
                    logger.warning(f"Error in optimized {query_name} query: {str(e)}")
            
            # Calculate averages and improvement
            if original_times and optimized_times:
                avg_original = statistics.mean(original_times)
                avg_optimized = statistics.mean(optimized_times)
                improvement = (avg_original - avg_optimized) / avg_original * 100
                
                logger.info(f"{query_name} - Avg original: {avg_original:.6f}s, " +
                           f"Avg optimized: {avg_optimized:.6f}s, " +
                           f"Improvement: {improvement:.2f}%")
                
                # Record results
                self.results['query_times'].append({
                    'name': query_name,
                    'original_time': avg_original,
                    'optimized_time': avg_optimized,
                    'improvement': improvement
                })
    
    def _test_materialized_view_refresh(self):
        """Test materialized view refresh performance."""
        logger.info("Testing materialized view refresh performance")
        
        try:
            # Time the refresh
            start_time = time.time()
            execute_sql("SELECT refresh_toxicity_materialized_views()")
            refresh_time = time.time() - start_time
            
            logger.info(f"Materialized view refresh completed in {refresh_time:.2f}s")
            
            # Record result
            self.results['materialized_view_refresh'] = refresh_time
        except Exception as e:
            logger.error(f"Error refreshing materialized views: {str(e)}")
    
    def _print_summary(self):
        """Print performance test summary."""
        print("\n" + "=" * 80)
        print("TOXICITY DATABASE PERFORMANCE TEST SUMMARY")
        print("=" * 80)
        
        # Print query performance results
        print("\nQUERY PERFORMANCE:")
        print("-" * 80)
        print(f"{'Query':<30} {'Original (s)':<15} {'Optimized (s)':<15} {'Improvement':<15}")
        print("-" * 80)
        
        for result in self.results['query_times']:
            print(f"{result['name']:<30} {result['original_time']:<15.6f} "
                  f"{result['optimized_time']:<15.6f} {result['improvement']:.2f}%")
        
        # Print materialized view refresh time
        print("\nMATERIALIZED VIEW REFRESH:")
        print("-" * 80)
        
        if self.results['materialized_view_refresh'] is not None:
            print(f"Refresh time: {self.results['materialized_view_refresh']:.2f}s")
        else:
            print("Refresh test failed or was skipped")
        
        # Print average improvement
        if self.results['query_times']:
            avg_improvement = statistics.mean([r['improvement'] for r in self.results['query_times']])
            print("\nAVERAGE QUERY PERFORMANCE IMPROVEMENT: {:.2f}%".format(avg_improvement))
        
        print("\n" + "=" * 80)
    
    def export_results(self, filename):
        """Export test results to file."""
        if not filename:
            return
        
        logger.info(f"Exporting results to {filename}")
        
        try:
            with open(filename, 'w') as f:
                json.dump(self.results, f, indent=2)
            logger.info(f"Results exported to {filename}")
        except Exception as e:
            logger.error(f"Failed to export results: {str(e)}")

def main():
    parser = argparse.ArgumentParser(description="Test toxicity database performance")
    parser.add_argument("--detailed", action="store_true", help="Show detailed test output")
    parser.add_argument("--export", help="Export results to file (JSON)")
    args = parser.parse_args()
    
    try:
        # Run tests
        test = DBPerformanceTest(detailed=args.detailed)
        results = test.run_all_tests()
        
        # Export results if requested
        if args.export:
            test.export_results(args.export)
        
        # Success if we have results
        return 0 if results['query_times'] else 1
        
    except Exception as e:
        logger.error(f"Test failed: {str(e)}")
        return 1

if __name__ == "__main__":
    sys.exit(main())