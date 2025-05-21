#!/usr/bin/env python3
"""
CryoProtect - Simple Toxicity Database Performance Test

This script tests database-level performance of the toxicity optimization by
directly executing SQL queries and comparing performance between original
and optimized schema approaches.

Usage:
    python test_toxicity_db_simple.py
"""

import os
import sys
import time
import json
import statistics
import logging
from datetime import datetime

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Add project root to Python path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

# Import database module
from database.db import execute_query, init_connection_pool

def time_query(query, params=None):
    """Execute a query and measure execution time."""
    start_time = time.time()
    result = execute_query(query, params)
    end_time = time.time()
    
    return result, end_time - start_time

def main():
    # Initialize database connection
    if not init_connection_pool():
        logger.error("Failed to initialize database connection")
        return 1
    
    print("Testing toxicity database performance...")
    
    # Get test molecule IDs
    try:
        molecule_result = execute_query("SELECT id FROM molecules LIMIT 5")
        if not molecule_result:
            logger.error("No molecules found in database")
            return 1
        
        molecule_ids = [row['id'] for row in molecule_result]
        logger.info(f"Found {len(molecule_ids)} test molecules")
    except Exception as e:
        logger.error(f"Error getting test molecules: {str(e)}")
        return 1
    
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
        }
    ]
    
    # Test query performance
    results = []
    for query_set in test_queries:
        print(f"Testing {query_set['name']}...")
        query_name = query_set['name']
        original_times = []
        optimized_times = []
        
        for molecule_id in molecule_ids:
            # Test original query
            try:
                _, orig_time = time_query(query_set['original'], [molecule_id])
                logger.info(f"Original query for {molecule_id}: {orig_time:.6f}s")
                original_times.append(orig_time)
            except Exception as e:
                logger.warning(f"Error in original query: {str(e)}")
            
            # Test optimized query
            try:
                _, opt_time = time_query(query_set['optimized'], [molecule_id])
                logger.info(f"Optimized query for {molecule_id}: {opt_time:.6f}s")
                optimized_times.append(opt_time)
            except Exception as e:
                logger.warning(f"Error in optimized query: {str(e)}")
        
        # Calculate averages and improvement
        if original_times and optimized_times:
            avg_original = statistics.mean(original_times)
            avg_optimized = statistics.mean(optimized_times)
            improvement = (avg_original - avg_optimized) / avg_original * 100
            
            result = {
                'name': query_name,
                'avg_original': avg_original,
                'avg_optimized': avg_optimized,
                'improvement': improvement,
                'sample_count': len(original_times)
            }
            results.append(result)
    
    # Print results
    print("\n" + "=" * 60)
    print("TOXICITY DATABASE PERFORMANCE RESULTS")
    print("=" * 60)
    print(f"{'Query':<30} {'Original (s)':<12} {'Optimized (s)':<12} {'Improvement':<12}")
    print("-" * 60)
    
    for result in results:
        print(f"{result['name']:<30} {result['avg_original']:<12.6f} "
              f"{result['avg_optimized']:<12.6f} {result['improvement']:.2f}%")
    
    # Calculate overall improvement
    if results:
        avg_improvement = statistics.mean([r['improvement'] for r in results])
        print("\nAVERAGE PERFORMANCE IMPROVEMENT: {:.2f}%".format(avg_improvement))
    
    # Test materialized view refresh
    print("\nTesting materialized view refresh...")
    try:
        start_time = time.time()
        execute_query("SELECT refresh_toxicity_materialized_views()")
        refresh_time = time.time() - start_time
        print(f"Materialized view refresh completed in {refresh_time:.2f}s")
    except Exception as e:
        logger.error(f"Error refreshing materialized views: {str(e)}")
    
    # Export results
    results_file = f"toxicity_db_results_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json"
    try:
        with open(results_file, 'w') as f:
            json.dump({
                'query_results': results,
                'timestamp': datetime.now().isoformat()
            }, f, indent=2)
        print(f"\nResults saved to {results_file}")
    except Exception as e:
        logger.error(f"Error saving results: {str(e)}")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())