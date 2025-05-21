#!/usr/bin/env python3
"""
CryoProtect - Toxicity Database Performance Test

This script tests database-level performance of the toxicity optimization by
directly executing SQL queries and comparing performance between original
and optimized schema approaches.

The tests include:
1. Query execution time comparison
2. Query plan analysis
3. Index usage verification
4. Materialized view performance

Usage:
    python test_toxicity_db_performance.py [--detailed] [--export FILENAME]
"""

import os
import sys
import json
import time
import argparse
import logging
import statistics
import pandas as pd
import matplotlib.pyplot as plt
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
    from database.utils import execute_sql, execute_sql_with_plan
except ImportError:
    logger.error("Failed to import database utilities. Make sure you're in the correct directory.")
    sys.exit(1)

class DBPerformanceTest:
    """Database performance test for toxicity optimization."""
    
    def __init__(self, detailed=False):
        """Initialize the test suite."""
        self.detailed = detailed
        self.results = {
            'query_times': [],
            'index_usage': [],
            'materialized_view_refresh': None
        }
        self.test_molecule_ids = self._get_test_molecule_ids()
    
    def _get_test_molecule_ids(self, limit=5):
        """Get sample molecule IDs for testing."""
        try:
            query = f"SELECT id FROM molecules LIMIT {limit}"
            result = execute_sql(query)
            return [row['id'] for row in result]
        except Exception as e:
            logger.error(f"Failed to get test molecule IDs: {str(e)}")
            return []
    
    def run_all_tests(self):
        """Run all performance tests."""
        logger.info("Starting database performance tests")
        
        # Run query performance tests
        self._test_query_performance()
        
        # Test index usage
        self._test_index_usage()
        
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
                    start_time = time.time()
                    if query_set['name'] == 'Toxicity Score Calculation':
                        # Special case for toxicity score which needs multiple param references
                        execute_sql(query_set['original'], [molecule_id, molecule_id, molecule_id])
                    else:
                        execute_sql(query_set['original'], [molecule_id])
                    original_time = time.time() - start_time
                    original_times.append(original_time)
                    
                    if self.detailed:
                        logger.info(f"Original {query_name} query for {molecule_id}: {original_time:.6f}s")
                except Exception as e:
                    logger.warning(f"Error in original {query_name} query: {str(e)}")
                
                # Test optimized query
                try:
                    start_time = time.time()
                    execute_sql(query_set['optimized'], [molecule_id])
                    optimized_time = time.time() - start_time
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
    
    def _test_index_usage(self):
        """Test index usage in optimized schema."""
        if not self.test_molecule_ids:
            logger.warning("No test molecule IDs available, skipping index usage test")
            return
        
        logger.info("Testing index usage")
        
        # Select first test molecule
        molecule_id = self.test_molecule_ids[0]
        
        # Define queries to test
        test_queries = [
            {
                'name': 'Toxicity Summary Index',
                'query': "SELECT * FROM toxicity_summary WHERE molecule_id = %s",
                'expected_index': 'idx_toxicity_summary_molecule_id'
            },
            {
                'name': 'LD50 Summary Species Index',
                'query': "SELECT * FROM ld50_summary WHERE species = 'Rat' AND molecule_id = %s",
                'expected_index': 'idx_ld50_summary_species'
            },
            {
                'name': 'Tox21 Activity Outcome Index',
                'query': "SELECT * FROM tox21_activity_summary WHERE activity_outcome = 'Active' AND molecule_id = %s",
                'expected_index': 'idx_tox21_summary_activity'
            },
            {
                'name': 'Hazard Classification System Index',
                'query': "SELECT * FROM hazard_classification_summary WHERE classification_system = 'GHS' AND molecule_id = %s",
                'expected_index': 'idx_hazard_summary_system'
            }
        ]
        
        # Check each query's execution plan
        for query_info in test_queries:
            try:
                # Get execution plan
                plan = execute_sql_with_plan(query_info['query'], [molecule_id])
                
                # Check for index scan in plan
                plan_str = json.dumps(plan)
                index_used = query_info['expected_index'] in plan_str
                scan_type = 'Index Scan' if 'Index Scan' in plan_str else 'Sequential Scan'
                
                # Log result
                if index_used:
                    logger.info(f"{query_info['name']}: Successfully using {query_info['expected_index']}")
                else:
                    logger.warning(f"{query_info['name']}: Not using expected index {query_info['expected_index']}")
                    if self.detailed:
                        logger.info(f"Plan: {json.dumps(plan, indent=2)}")
                
                # Record result
                self.results['index_usage'].append({
                    'name': query_info['name'],
                    'expected_index': query_info['expected_index'],
                    'index_used': index_used,
                    'scan_type': scan_type
                })
            except Exception as e:
                logger.error(f"Error testing index usage for {query_info['name']}: {str(e)}")
    
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
        
        # Print index usage results
        print("\nINDEX USAGE:")
        print("-" * 80)
        print(f"{'Query':<30} {'Expected Index':<35} {'Used':<10}")
        print("-" * 80)
        
        for result in self.results['index_usage']:
            print(f"{result['name']:<30} {result['expected_index']:<35} "
                  f"{'Yes' if result['index_used'] else 'No':<10}")
        
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
            
            # Try to generate charts if pandas and matplotlib are available
            self._generate_charts(filename.replace('.json', '_charts.png'))
            
            logger.info(f"Results exported to {filename}")
        except Exception as e:
            logger.error(f"Failed to export results: {str(e)}")
    
    def _generate_charts(self, chart_filename):
        """Generate performance charts."""
        try:
            # Create performance comparison chart
            if self.results['query_times']:
                plt.figure(figsize=(12, 8))
                
                # Query performance chart
                df = pd.DataFrame(self.results['query_times'])
                
                # Bar chart
                ax = plt.subplot(2, 1, 1)
                index = range(len(df))
                bar_width = 0.35
                
                original_bars = ax.bar(index, df['original_time'], bar_width, 
                                     label='Original', color='#ff9999')
                optimized_bars = ax.bar([i + bar_width for i in index], df['optimized_time'], 
                                      bar_width, label='Optimized', color='#99ff99')
                
                ax.set_xlabel('Query')
                ax.set_ylabel('Time (seconds)')
                ax.set_title('Query Execution Time Comparison')
                ax.set_xticks([i + bar_width/2 for i in index])
                ax.set_xticklabels(df['name'], rotation=45, ha='right')
                ax.legend()
                
                # Improvement chart
                ax2 = plt.subplot(2, 1, 2)
                ax2.bar(df['name'], df['improvement'], color='#9999ff')
                ax2.set_xlabel('Query')
                ax2.set_ylabel('Improvement (%)')
                ax2.set_title('Performance Improvement')
                plt.xticks(rotation=45, ha='right')
                
                plt.tight_layout()
                plt.savefig(chart_filename)
                
                logger.info(f"Performance charts saved to {chart_filename}")
        except Exception as e:
            logger.warning(f"Failed to generate charts: {str(e)}")

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