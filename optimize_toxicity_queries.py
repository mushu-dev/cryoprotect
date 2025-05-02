#!/usr/bin/env python3
"""
Optimize database queries for toxicity data.

This script:
1. Adds additional indexes for commonly queried toxicity data
2. Creates materialized views for frequently accessed data
3. Verifies query performance with explain analyze
"""

import os
import sys
import logging
import time
import argparse
from supabase import create_client
from config import Config

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def add_indexes(supabase):
    """Add performance indexes for toxicity data."""
    logger.info("Adding performance indexes for toxicity data...")
    
    indexes = [
        # Compound indexes for commonly used query patterns
        "CREATE INDEX IF NOT EXISTS idx_toxicity_data_molecule_hit_call ON toxicity_data (molecule_id, hit_call)",
        "CREATE INDEX IF NOT EXISTS idx_toxicity_data_assay_activity ON toxicity_data (assay_id, activity_value)",
        "CREATE INDEX IF NOT EXISTS idx_toxicity_score_molecule_type ON toxicity_score (molecule_id, score_type)",
        "CREATE INDEX IF NOT EXISTS idx_assay_endpoint_relevance ON assay_endpoint_mapping (assay_id, endpoint_id, relevance_score)",
        "CREATE INDEX IF NOT EXISTS idx_toxicity_assay_endpoint ON toxicity_assay (toxicological_endpoint)"
    ]
    
    success_count = 0
    for index_sql in indexes:
        try:
            logger.info(f"Creating index: {index_sql}")
            # Execute the index creation
            supabase.postgrest.rpc('exec_sql', {'query': index_sql}).execute()
            success_count += 1
        except Exception as e:
            logger.error(f"Error creating index: {str(e)}")
    
    logger.info(f"Successfully created {success_count}/{len(indexes)} indexes")
    return success_count == len(indexes)

def create_materialized_views(supabase):
    """Create materialized views for frequently accessed toxicity data."""
    logger.info("Creating materialized views for toxicity data...")
    
    views = [
        """
        CREATE MATERIALIZED VIEW IF NOT EXISTS toxicity_molecule_summary AS
        SELECT 
            m.id AS molecule_id,
            m.name AS molecule_name,
            m.toxicity_score,
            COUNT(td.id) AS total_assays,
            SUM(CASE WHEN td.hit_call = true THEN 1 ELSE 0 END) AS active_assays,
            AVG(td.reliability_score) AS avg_reliability
        FROM 
            molecule m
        LEFT JOIN 
            toxicity_data td ON m.id = td.molecule_id
        GROUP BY 
            m.id, m.name, m.toxicity_score
        """,
        
        """
        CREATE MATERIALIZED VIEW IF NOT EXISTS toxicity_endpoint_summary AS
        SELECT 
            te.id AS endpoint_id,
            te.name AS endpoint_name,
            te.category,
            COUNT(DISTINCT ta.id) AS assay_count,
            COUNT(DISTINCT td.molecule_id) AS molecule_count
        FROM 
            toxicity_endpoint te
        LEFT JOIN 
            assay_endpoint_mapping aem ON te.id = aem.endpoint_id
        LEFT JOIN 
            toxicity_assay ta ON aem.assay_id = ta.id
        LEFT JOIN 
            toxicity_data td ON ta.id = td.assay_id
        GROUP BY 
            te.id, te.name, te.category
        """
    ]
    
    success_count = 0
    for view_sql in views:
        try:
            logger.info(f"Creating materialized view: {view_sql.split('AS')[0]}")
            # Execute the view creation
            supabase.postgrest.rpc('exec_sql', {'query': view_sql}).execute()
            success_count += 1
        except Exception as e:
            logger.error(f"Error creating materialized view: {str(e)}")
    
    logger.info(f"Successfully created {success_count}/{len(views)} materialized views")
    return success_count == len(views)

def test_query_performance(supabase):
    """Test the performance of common toxicity data queries."""
    logger.info("Testing query performance...")
    
    queries = [
        # Common query patterns
        "SELECT * FROM toxicity_data WHERE molecule_id = '00000000-0000-0000-0000-000000000000' LIMIT 1",
        "SELECT * FROM toxicity_score WHERE molecule_id = '00000000-0000-0000-0000-000000000000' LIMIT 1",
        "SELECT * FROM toxicity_assay WHERE toxicological_endpoint = 'nuclear_receptor' LIMIT 10"
    ]
    
    # Replace the dummy UUID with a real molecule ID
    try:
        molecule_response = supabase.table("molecule").select("id").limit(1).execute()
        if molecule_response.data and len(molecule_response.data) > 0:
            real_id = molecule_response.data[0]["id"]
            queries = [q.replace('00000000-0000-0000-0000-000000000000', real_id) for q in queries]
    except Exception as e:
        logger.warning(f"Couldn't get real molecule ID: {str(e)}")
    
    for query in queries:
        try:
            logger.info(f"Testing query: {query}")
            # Execute EXPLAIN ANALYZE
            explain_query = f"EXPLAIN ANALYZE {query}"
            start_time = time.time()
            response = supabase.postgrest.rpc('exec_sql', {'query': explain_query}).execute()
            duration = time.time() - start_time
            
            if hasattr(response, 'data') and response.data:
                logger.info("Query execution plan:")
                for line in response.data:
                    logger.info(f"  {line}")
            
            logger.info(f"Query executed in {duration:.4f} seconds")
        except Exception as e:
            logger.error(f"Error testing query: {str(e)}")
    
    # Also test the materialized views
    materialized_view_queries = [
        "SELECT * FROM toxicity_molecule_summary LIMIT 10",
        "SELECT * FROM toxicity_endpoint_summary LIMIT 10"
    ]
    
    for query in materialized_view_queries:
        try:
            logger.info(f"Testing materialized view query: {query}")
            start_time = time.time()
            response = supabase.postgrest.rpc('exec_sql', {'query': query}).execute()
            duration = time.time() - start_time
            
            if hasattr(response, 'data') and response.data:
                row_count = len(response.data)
                logger.info(f"Query returned {row_count} rows in {duration:.4f} seconds")
            else:
                logger.info(f"Query executed in {duration:.4f} seconds (no data returned)")
        except Exception as e:
            logger.error(f"Error testing materialized view query: {str(e)}")
    
    return True

def refresh_materialized_views(supabase):
    """Refresh the materialized views to ensure they have current data."""
    logger.info("Refreshing materialized views...")
    
    views = [
        "toxicity_molecule_summary",
        "toxicity_endpoint_summary"
    ]
    
    success_count = 0
    for view in views:
        try:
            refresh_sql = f"REFRESH MATERIALIZED VIEW {view}"
            logger.info(f"Refreshing view: {view}")
            supabase.postgrest.rpc('exec_sql', {'query': refresh_sql}).execute()
            success_count += 1
        except Exception as e:
            logger.error(f"Error refreshing materialized view {view}: {str(e)}")
    
    logger.info(f"Successfully refreshed {success_count}/{len(views)} materialized views")
    return success_count == len(views)

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='Optimize toxicity queries')
    parser.add_argument('--skip-indexes', action='store_true', help='Skip creating indexes')
    parser.add_argument('--skip-views', action='store_true', help='Skip creating materialized views')
    parser.add_argument('--skip-tests', action='store_true', help='Skip query performance tests')
    parser.add_argument('--refresh-only', action='store_true', help='Only refresh existing materialized views')
    return parser.parse_args()

def main():
    """Main function to optimize toxicity queries."""
    args = parse_arguments()
    
    # Initialize Supabase client
    config = Config()
    supabase_url = os.environ.get("SUPABASE_URL") or config.SUPABASE_URL
    supabase_key = os.environ.get("SUPABASE_KEY") or os.environ.get("SUPABASE_SERVICE_KEY") or config.SUPABASE_KEY or config.SUPABASE_SERVICE_KEY
    
    if not supabase_url or not supabase_key:
        logger.error("Supabase URL or key not found")
        logger.error("Please ensure SUPABASE_URL and either SUPABASE_KEY or SUPABASE_SERVICE_KEY are set in your environment or config.py")
        return False
    
    supabase = create_client(supabase_url, supabase_key)
    
    if args.refresh_only:
        # Only refresh the materialized views
        return refresh_materialized_views(supabase)
    
    # Add indexes
    if not args.skip_indexes:
        if not add_indexes(supabase):
            logger.warning("Failed to add all indexes")
    else:
        logger.info("Skipping index creation")
    
    # Create materialized views
    if not args.skip_views:
        if not create_materialized_views(supabase):
            logger.warning("Failed to create all materialized views")
    else:
        logger.info("Skipping materialized view creation")
    
    # Test query performance
    if not args.skip_tests:
        test_query_performance(supabase)
    else:
        logger.info("Skipping query performance tests")
    
    logger.info("Toxicity query optimization completed")
    return True

if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)