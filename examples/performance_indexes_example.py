#!/usr/bin/env python3
"""
CryoProtect v2 - Database Performance Optimization Example

This script demonstrates how to use the add_performance_indexes.py script
to optimize database performance by creating appropriate indexes.

Usage:
    python examples/performance_indexes_example.py [--test-only] [--schema SCHEMA]

Options:
    --test-only     Only test performance without creating indexes
    --schema        Database schema to modify (default: public)
"""

import os
import sys
import argparse
import logging
from datetime import datetime

# Add the parent directory to the path so we can import the script
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# Import the script
import add_performance_indexes

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.StreamHandler()
    ]
)
logger = logging.getLogger("performance_example")

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="Example script for database performance optimization")
    parser.add_argument("--schema", default="public",
                        help="Database schema to modify (default: public)")
    parser.add_argument("--test-only", action="store_true",
                        help="Only test the performance without applying indexes")
    return parser.parse_args()

def main():
    """Main function."""
    try:
        # Parse arguments
        args = parse_arguments()
        
        logger.info("Starting database performance optimization example")
        logger.info(f"Schema: {args.schema}")
        logger.info(f"Test only: {args.test_only}")
        
        # Step 1: Run performance tests before optimization
        logger.info("Step 1: Running initial performance tests")
        before_metrics = add_performance_indexes.run_performance_test()
        before_summary = before_metrics.get_summary()
        
        # Print initial performance metrics
        logger.info("Initial performance metrics:")
        for operation, stats in before_summary["operations"].items():
            logger.info(f"  {operation}: avg={stats['avg_ms']:.2f}ms, p95={stats['p95_ms']:.2f}ms")
        
        # Step 2: Analyze metrics and generate recommendations
        logger.info("Step 2: Analyzing performance metrics")
        analysis = add_performance_indexes.analyze_metrics(before_summary)
        recommendations = add_performance_indexes.generate_recommendations(before_summary)
        
        # Print bottlenecks and recommendations
        if analysis["slow_operations"]:
            logger.info("Identified slow operations:")
            for op in analysis["slow_operations"]:
                logger.info(f"  {op['operation']}: p95={op['p95_ms']:.2f}ms")
        
        if recommendations:
            logger.info("Recommendations:")
            for rec in recommendations:
                if "operation" in rec:
                    logger.info(f"  {rec['operation']}: {rec['recommendation']}")
                else:
                    logger.info(f"  {rec['issue']}: {rec['recommendation']}")
        
        # Step 3: Apply performance indexes if not in test-only mode
        if not args.test_only:
            logger.info("Step 3: Applying performance indexes")
            
            # Get SQL statements for creating indexes
            index_batches = add_performance_indexes.get_performance_indexes_sql(args.schema)
            
            # Print the SQL statements that would be executed
            logger.info("SQL statements for creating indexes:")
            for batch_name, batch_sql in index_batches.items():
                logger.info(f"  {batch_name}:")
                for line in batch_sql.split("\n"):
                    if "CREATE INDEX" in line:
                        logger.info(f"    {line.strip()}")
            
            # Apply the indexes
            logger.info("Applying indexes...")
            results = add_performance_indexes.apply_performance_indexes(
                schema=args.schema,
                test_only=False
            )
            
            # Print results
            logger.info(f"Created {len(results['created_indexes'])} indexes")
            if results["failed_indexes"]:
                logger.warning(f"Failed to create {len(results['failed_indexes'])} indexes")
                for index in results["failed_indexes"]:
                    logger.warning(f"  Failed index: {index}")
            
            # Step 4: Run performance tests after optimization
            logger.info("Step 4: Running performance tests after optimization")
            after_metrics = add_performance_indexes.run_performance_test()
            after_summary = after_metrics.get_summary()
            
            # Calculate improvements
            improvements = add_performance_indexes.calculate_improvements({
                "before_metrics": before_summary,
                "after_metrics": after_summary
            })
            
            # Print performance improvements
            logger.info("Performance improvements:")
            for op_name, improvement in improvements.items():
                logger.info(f"  {op_name}:")
                logger.info(f"    Avg: {improvement['avg_ms_before']:.2f}ms -> {improvement['avg_ms_after']:.2f}ms ({improvement['avg_improvement_pct']:.2f}% improvement)")
                logger.info(f"    P95: {improvement['p95_ms_before']:.2f}ms -> {improvement['p95_ms_after']:.2f}ms ({improvement['p95_improvement_pct']:.2f}% improvement)")
        else:
            logger.info("Test-only mode, skipping index creation")
        
        logger.info("Performance optimization example completed")
        
        return 0
    
    except Exception as e:
        logger.error(f"Error in example script: {str(e)}", exc_info=True)
        return 1

if __name__ == "__main__":
    sys.exit(main())