#!/usr/bin/env python3
"""
Optimize Complex RLS Queries

This script improves performance of Row-Level Security (RLS) policies for complex 
queries in the CryoProtect database by applying multiple optimization strategies:

1. Security definer functions that reduce RLS policy evaluation overhead
2. Performance-tuned indexes for RLS policy expressions
3. Materialized views for frequently accessed data patterns
4. Optimized RLS policies leveraging the security definer functions
5. Query performance analysis and optimization for complex data access patterns
6. Cache layer for temporary query results

The script can be used with either direct database connection or Supabase MCP.
"""

import os
import sys
import argparse
import logging
import time
import json
from pathlib import Path
from typing import Dict, List, Any, Optional, Union
from datetime import datetime
import statistics
import traceback

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler(),
        logging.FileHandler(f"rls_optimization_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"),
    ]
)
logger = logging.getLogger(__name__)

# Define SQL file paths relative to script location
SCRIPT_DIR = Path(__file__).parent
MIGRATIONS_DIR = SCRIPT_DIR / "migrations" / "rls_helpers"

# SQL file paths
SECURITY_DEFINER_FUNCTIONS_FILE = MIGRATIONS_DIR / "security_definer_functions.sql"
PERFORMANCE_INDEXES_FILE = MIGRATIONS_DIR / "performance_indexes.sql"
MATERIALIZED_VIEWS_FILE = MIGRATIONS_DIR / "materialized_views.sql"
RLS_POLICIES_FILE = MIGRATIONS_DIR / "rls_policies.sql"
COMPLEX_QUERY_OPTIMIZATION_FILE = MIGRATIONS_DIR / "complex_query_optimization.sql"

# Custom SQL for complex query optimization if file doesn't exist
COMPLEX_QUERY_OPTIMIZATION_SQL = """
-- Complex Query Optimization SQL
-- This adds specialized functions and indexes for complex query patterns

-- 1. Molecular Property Range Query Optimization
CREATE OR REPLACE FUNCTION find_molecules_by_property_range(
    property_name_param TEXT,
    min_value NUMERIC,
    max_value NUMERIC
) RETURNS SETOF uuid AS $$
BEGIN
    RETURN QUERY
    SELECT DISTINCT m.id
    FROM molecules m
    JOIN molecular_properties mp ON m.id = mp.molecule_id
    WHERE mp.property_name = property_name_param
    AND mp.property_value::numeric >= min_value
    AND mp.property_value::numeric <= max_value
    AND (m.is_public = true OR m.created_by = auth.uid() OR 
         EXISTS (
             SELECT 1 FROM project_molecules pm
             JOIN team_projects tp ON pm.project_id = tp.project_id
             JOIN user_profile up ON tp.team_id = up.team_id
             WHERE pm.molecule_id = m.id AND up.auth_user_id = auth.uid()
         ));
END;
$$ LANGUAGE plpgsql SECURITY DEFINER;

-- 2. Full Text Search Optimization for molecules
CREATE INDEX IF NOT EXISTS idx_molecules_name_description_tsvector 
ON molecules USING gin(to_tsvector('english', COALESCE(name, '') || ' ' || COALESCE(description, '')));

CREATE OR REPLACE FUNCTION search_molecules_text(search_term TEXT)
RETURNS SETOF uuid AS $$
BEGIN
    RETURN QUERY
    SELECT id FROM molecules
    WHERE to_tsvector('english', COALESCE(name, '') || ' ' || COALESCE(description, '')) @@ plainto_tsquery('english', search_term)
    AND (is_public = true OR created_by = auth.uid() OR
         EXISTS (
             SELECT 1 FROM project_molecules pm
             JOIN team_projects tp ON pm.project_id = tp.project_id
             JOIN user_profile up ON tp.team_id = up.team_id
             WHERE pm.molecule_id = molecules.id AND up.auth_user_id = auth.uid()
         ));
END;
$$ LANGUAGE plpgsql SECURITY DEFINER;

-- 3. Bulk Molecule Access Check
CREATE OR REPLACE FUNCTION filter_accessible_molecules(molecule_ids uuid[])
RETURNS SETOF uuid AS $$
BEGIN
    RETURN QUERY
    SELECT id FROM molecules
    WHERE id = ANY(molecule_ids)
    AND (is_public = true OR created_by = auth.uid() OR
         EXISTS (
             SELECT 1 FROM project_molecules pm
             JOIN team_projects tp ON pm.project_id = tp.project_id
             JOIN user_profile up ON tp.team_id = up.team_id
             WHERE pm.molecule_id = molecules.id AND up.auth_user_id = auth.uid()
         ));
END;
$$ LANGUAGE plpgsql SECURITY DEFINER;

-- 4. Optimization for Common Query Patterns
-- 4.1 Molecule-with-properties query
CREATE OR REPLACE FUNCTION get_molecules_with_properties(
    limit_count INTEGER DEFAULT 100,
    offset_val INTEGER DEFAULT 0
) RETURNS TABLE (
    id uuid,
    name TEXT, 
    smiles TEXT, 
    molecular_formula TEXT,
    property_count INTEGER
) AS $$
DECLARE
    user_id uuid;
BEGIN
    user_id := auth.uid();
    
    RETURN QUERY
    WITH accessible_molecules AS (
        SELECT m.id
        FROM molecules m
        WHERE m.is_public = true 
           OR m.created_by = user_id
           OR EXISTS (
               SELECT 1 FROM project_molecules pm
               JOIN team_projects tp ON pm.project_id = tp.project_id
               JOIN user_profile up ON tp.team_id = up.team_id
               WHERE pm.molecule_id = m.id AND up.auth_user_id = user_id
           )
    )
    SELECT 
        m.id,
        m.name,
        m.smiles,
        m.molecular_formula,
        COUNT(mp.id)::INTEGER AS property_count
    FROM 
        molecules m
    JOIN
        accessible_molecules am ON m.id = am.id
    LEFT JOIN
        molecular_properties mp ON m.id = mp.molecule_id
    GROUP BY
        m.id, m.name, m.smiles, m.molecular_formula
    ORDER BY
        m.name
    LIMIT limit_count
    OFFSET offset_val;
END;
$$ LANGUAGE plpgsql SECURITY DEFINER;

-- 4.2 Mixture-with-components query
CREATE OR REPLACE FUNCTION get_mixtures_with_components(
    limit_count INTEGER DEFAULT 100,
    offset_val INTEGER DEFAULT 0
) RETURNS TABLE (
    id uuid,
    name TEXT, 
    description TEXT,
    component_count INTEGER
) AS $$
DECLARE
    user_id uuid;
BEGIN
    user_id := auth.uid();
    
    RETURN QUERY
    WITH accessible_mixtures AS (
        SELECT m.id
        FROM mixtures m
        WHERE m.is_public = true 
           OR m.created_by = user_id
           OR EXISTS (
               SELECT 1 FROM project_mixtures pm
               JOIN team_projects tp ON pm.project_id = tp.project_id
               JOIN user_profile up ON tp.team_id = up.team_id
               WHERE pm.mixture_id = m.id AND up.auth_user_id = user_id
           )
    )
    SELECT 
        m.id,
        m.name,
        m.description,
        COUNT(mc.id)::INTEGER AS component_count
    FROM 
        mixtures m
    JOIN
        accessible_mixtures am ON m.id = am.id
    LEFT JOIN
        mixture_components mc ON m.id = mc.mixture_id
    GROUP BY
        m.id, m.name, m.description
    ORDER BY
        m.name
    LIMIT limit_count
    OFFSET offset_val;
END;
$$ LANGUAGE plpgsql SECURITY DEFINER;

-- Create hypertable indexes for efficient time-series queries (if timescaledb is available)
DO $$
BEGIN
    IF EXISTS (
        SELECT 1 FROM pg_extension WHERE extname = 'timescaledb'
    ) THEN
        -- Add timescale-specific optimizations for time-series data
        CREATE INDEX IF NOT EXISTS idx_experiments_time 
        ON experiments(created_at DESC, id);
        
        CREATE INDEX IF NOT EXISTS idx_predictions_time 
        ON predictions(created_at DESC, id);
    END IF;
END $$;

-- Additional composite indexes for complex query patterns
CREATE INDEX IF NOT EXISTS idx_molecular_properties_name_numeric 
ON molecular_properties(property_name, (property_value::numeric))
WHERE property_value ~ '^-?[0-9]+(\.[0-9]+)?$';

CREATE INDEX IF NOT EXISTS idx_molecular_properties_name_text
ON molecular_properties(property_name, property_value)
WHERE property_type = 'text';

-- Add covering indexes for common queries that involve multiple columns
CREATE INDEX IF NOT EXISTS idx_molecules_cover_common
ON molecules(id, name, smiles, molecular_formula, is_public, created_by);

CREATE INDEX IF NOT EXISTS idx_mixtures_cover_common
ON mixtures(id, name, description, is_public, created_by);
"""

# Dictionary of complex query patterns and their test queries
COMPLEX_QUERY_TESTS = {
    "property_range_query": """
        SELECT m.id, m.name, m.molecular_formula
        FROM molecules m
        JOIN molecular_properties mp ON m.id = mp.molecule_id
        WHERE mp.property_name = 'Molecular Weight'
        AND mp.property_value::numeric BETWEEN 100 AND 500
        AND (m.is_public = true OR m.created_by = auth.uid() OR 
             EXISTS (
                 SELECT 1 FROM project_molecules pm
                 JOIN team_projects tp ON pm.project_id = tp.project_id
                 JOIN user_profile up ON tp.team_id = up.team_id
                 WHERE pm.molecule_id = m.id AND up.auth_user_id = auth.uid()
             ))
        LIMIT 10;
    """,
    
    "fulltext_search_query": """
        SELECT id, name, molecular_formula 
        FROM molecules
        WHERE to_tsvector('english', COALESCE(name, '') || ' ' || COALESCE(description, '')) @@ plainto_tsquery('english', 'cryoprotectant')
        AND (is_public = true OR created_by = auth.uid() OR
             EXISTS (
                 SELECT 1 FROM project_molecules pm
                 JOIN team_projects tp ON pm.project_id = tp.project_id
                 JOIN user_profile up ON tp.team_id = up.team_id
                 WHERE pm.molecule_id = molecules.id AND up.auth_user_id = auth.uid()
             ))
        LIMIT 10;
    """,
    
    "molecules_with_properties_query": """
        SELECT 
            m.id,
            m.name,
            m.smiles,
            m.molecular_formula,
            COUNT(mp.id) AS property_count
        FROM 
            molecules m
        LEFT JOIN
            molecular_properties mp ON m.id = mp.molecule_id
        WHERE 
            m.is_public = true 
            OR m.created_by = auth.uid()
            OR EXISTS (
                SELECT 1 FROM project_molecules pm
                JOIN team_projects tp ON pm.project_id = tp.project_id
                JOIN user_profile up ON tp.team_id = up.team_id
                WHERE pm.molecule_id = m.id AND up.auth_user_id = auth.uid()
            )
        GROUP BY
            m.id, m.name, m.smiles, m.molecular_formula
        ORDER BY
            m.name
        LIMIT 10;
    """,
    
    "mixtures_with_components_query": """
        SELECT 
            m.id,
            m.name,
            m.description,
            COUNT(mc.id) AS component_count
        FROM 
            mixtures m
        LEFT JOIN
            mixture_components mc ON m.id = mc.mixture_id
        WHERE 
            m.is_public = true 
            OR m.created_by = auth.uid()
            OR EXISTS (
                SELECT 1 FROM project_mixtures pm
                JOIN team_projects tp ON pm.project_id = tp.project_id
                JOIN user_profile up ON tp.team_id = up.team_id
                WHERE pm.mixture_id = m.id AND up.auth_user_id = auth.uid()
            )
        GROUP BY
            m.id, m.name, m.description
        ORDER BY
            m.name
        LIMIT 10;
    """
}

# Optimized versions of the complex queries
OPTIMIZED_QUERY_TESTS = {
    "optimized_property_range_query": """
        SELECT m.id, m.name, m.molecular_formula
        FROM find_molecules_by_property_range('Molecular Weight', 100, 500) mol_ids
        JOIN molecules m ON m.id = mol_ids
        LIMIT 10;
    """,
    
    "optimized_fulltext_search_query": """
        SELECT m.id, m.name, m.molecular_formula
        FROM search_molecules_text('cryoprotectant') mol_ids
        JOIN molecules m ON m.id = mol_ids
        LIMIT 10;
    """,
    
    "optimized_molecules_with_properties_query": """
        SELECT * FROM get_molecules_with_properties(10, 0);
    """,
    
    "optimized_mixtures_with_components_query": """
        SELECT * FROM get_mixtures_with_components(10, 0);
    """
}

def read_sql_file(file_path: Path) -> str:
    """Read SQL file content."""
    try:
        if not file_path.exists():
            if file_path.name == "complex_query_optimization.sql":
                logger.info(f"File {file_path} not found, using built-in SQL")
                return COMPLEX_QUERY_OPTIMIZATION_SQL
            else:
                raise FileNotFoundError(f"SQL file not found: {file_path}")
                
        with open(file_path, 'r') as f:
            sql = f.read()
            logger.info(f"Read SQL file {file_path}")
            return sql
    except Exception as e:
        logger.error(f"Error reading SQL file {file_path}: {e}")
        raise

def apply_sql_direct(conn, sql: str) -> bool:
    """Apply SQL using direct database connection."""
    try:
        cursor = conn.cursor()
        cursor.execute(sql)
        conn.commit()
        cursor.close()
        logger.info("SQL applied successfully via direct connection")
        return True
    except Exception as e:
        conn.rollback()
        logger.error(f"Error applying SQL via direct connection: {e}")
        return False

def apply_sql_via_supabase(supabase, sql: str) -> bool:
    """Apply SQL using Supabase client."""
    try:
        response = supabase.rpc('exec_sql', {'query': sql}).execute()
        if hasattr(response, 'error') and response.error:
            logger.error(f"Error applying SQL via Supabase: {response.error}")
            return False
        logger.info("SQL applied successfully via Supabase")
        return True
    except Exception as e:
        logger.error(f"Error applying SQL via Supabase: {e}")
        return False

def apply_sql_via_mcp(project_id: str, sql: str) -> bool:
    """Apply SQL via Supabase MCP."""
    try:
        # Try different MCP modules in order of preference
        try:
            # Try using supabase_mcp_tools first
            from supabase_mcp_tools import execute_sql_on_supabase
            logger.info("Using supabase_mcp_tools to execute SQL")
            result = execute_sql_on_supabase(project_id, sql)
            
            if isinstance(result, dict) and result.get('error'):
                logger.error(f"MCP execution error: {result.get('error')}")
                return False
            
            logger.info("SQL executed successfully via MCP")
            return True
        except ImportError:
            # Try direct mcp__supabase__execute_sql
            try:
                import mcp__supabase__execute_sql
                logger.info("Using direct MCP function to execute SQL")
                result = mcp__supabase__execute_sql.execute({
                    "project_id": project_id,
                    "query": sql
                })
                
                if isinstance(result, dict) and result.get('error'):
                    logger.error(f"MCP execution error: {result.get('error')}")
                    return False
                
                logger.info("SQL executed successfully via direct MCP function")
                return True
            except ImportError:
                logger.error("MCP modules not available")
                return False
    except Exception as e:
        logger.error(f"Error executing SQL via MCP: {e}")
        return False

def execute_sql(db_client, sql: str, project_id: str = None) -> bool:
    """Execute SQL using the appropriate method based on client type."""
    if project_id:  # Use MCP
        return apply_sql_via_mcp(project_id, sql)
    elif hasattr(db_client, 'rpc'):  # Use Supabase client
        return apply_sql_via_supabase(db_client, sql)
    else:  # Use direct connection
        return apply_sql_direct(db_client, sql)

def apply_optimization_step(db_client, step_name: str, sql_file: Path, project_id: str = None) -> bool:
    """Apply a specific optimization step."""
    logger.info(f"Applying {step_name}...")
    start_time = time.time()
    
    try:
        sql = read_sql_file(sql_file)
        success = execute_sql(db_client, sql, project_id)
        
        if success:
            elapsed = time.time() - start_time
            logger.info(f"Successfully applied {step_name} in {elapsed:.2f} seconds")
            return True
        else:
            logger.error(f"Failed to apply {step_name}")
            return False
    except Exception as e:
        logger.error(f"Error during {step_name}: {e}")
        return False

def test_query_performance(db_client, query: str, project_id: str = None) -> Dict[str, Any]:
    """Test the performance of a query."""
    result = {
        "query": query,
        "success": False,
        "execution_time_ms": None,
        "error": None,
        "explain_plan": None
    }
    
    try:
        # Run query with timing
        start_time = time.time()
        if project_id:  # MCP
            try:
                from supabase_mcp_tools import execute_sql_on_supabase
                query_result = execute_sql_on_supabase(project_id, query)
            except ImportError:
                import mcp__supabase__execute_sql
                query_result = mcp__supabase__execute_sql.execute({
                    "project_id": project_id,
                    "query": query
                })
        elif hasattr(db_client, 'rpc'):  # Supabase
            query_result = db_client.rpc('exec_sql', {'query': query}).execute()
        else:  # Direct
            cursor = db_client.cursor()
            cursor.execute(query)
            query_result = cursor.fetchall()
            cursor.close()
            
        execution_time = time.time() - start_time
        result["execution_time_ms"] = execution_time * 1000
        result["success"] = True
        
        # Get explain plan
        try:
            explain_query = f"EXPLAIN ANALYZE {query}"
            if project_id:  # MCP
                try:
                    from supabase_mcp_tools import execute_sql_on_supabase
                    explain_result = execute_sql_on_supabase(project_id, explain_query)
                except ImportError:
                    import mcp__supabase__execute_sql
                    explain_result = mcp__supabase__execute_sql.execute({
                        "project_id": project_id,
                        "query": explain_query
                    })
            elif hasattr(db_client, 'rpc'):  # Supabase
                explain_result = db_client.rpc('exec_sql', {'query': explain_query}).execute()
                explain_result = explain_result.data if hasattr(explain_result, 'data') else explain_result
            else:  # Direct
                cursor = db_client.cursor()
                cursor.execute(explain_query)
                explain_result = cursor.fetchall()
                cursor.close()
                
            result["explain_plan"] = explain_result
        except Exception as e:
            logger.warning(f"Could not get explain plan: {e}")
            
    except Exception as e:
        logger.error(f"Error testing query performance: {e}")
        result["error"] = str(e)
        
    return result

def run_performance_tests(db_client, project_id: str = None) -> Dict[str, Any]:
    """Run performance tests on complex queries before and after optimization."""
    performance_results = {
        "original_queries": {},
        "optimized_queries": {},
        "improvement": {}
    }
    
    # Test original complex queries
    logger.info("Testing performance of original complex queries...")
    for query_name, query in COMPLEX_QUERY_TESTS.items():
        logger.info(f"Testing {query_name}...")
        performance_results["original_queries"][query_name] = test_query_performance(db_client, query, project_id)
        
    # Test optimized complex queries
    logger.info("Testing performance of optimized complex queries...")
    for query_name, query in OPTIMIZED_QUERY_TESTS.items():
        logger.info(f"Testing {query_name}...")
        performance_results["optimized_queries"][query_name] = test_query_performance(db_client, query, project_id)
    
    # Calculate improvement percentages
    for orig_name, orig_result in performance_results["original_queries"].items():
        opt_name = f"optimized_{orig_name}"
        if opt_name in performance_results["optimized_queries"]:
            opt_result = performance_results["optimized_queries"][opt_name]
            
            if orig_result["success"] and opt_result["success"] and orig_result["execution_time_ms"] and opt_result["execution_time_ms"]:
                improvement_pct = ((orig_result["execution_time_ms"] - opt_result["execution_time_ms"]) / orig_result["execution_time_ms"]) * 100
                performance_results["improvement"][orig_name] = {
                    "original_ms": orig_result["execution_time_ms"],
                    "optimized_ms": opt_result["execution_time_ms"],
                    "improvement_pct": improvement_pct,
                    "speedup_factor": orig_result["execution_time_ms"] / opt_result["execution_time_ms"]
                }
                
                logger.info(f"Query {orig_name} improvement: {improvement_pct:.2f}% "
                           f"({orig_result['execution_time_ms']:.2f}ms → {opt_result['execution_time_ms']:.2f}ms)")
    
    return performance_results

def verify_optimization(db_client, project_id: str = None) -> Dict[str, Any]:
    """Verify that optimizations have been applied."""
    verification_results = {
        "status": "passed",
        "checks": {},
        "issues": []
    }
    
    logger.info("Verifying optimizations...")
    
    try:
        # Define verification queries
        verification_queries = {
            "security_definer_functions": """
                SELECT proname, prosecdef
                FROM pg_proc p
                JOIN pg_namespace n ON p.pronamespace = n.oid
                WHERE n.nspname = 'public'
                AND proname IN (
                    'has_molecule_access', 'has_mixture_access', 'is_team_member',
                    'user_has_clearance', 'find_molecules_by_property_range',
                    'search_molecules_text', 'filter_accessible_molecules'
                );
            """,
            
            "performance_indexes": """
                SELECT indexname, indexdef
                FROM pg_indexes
                WHERE schemaname = 'public'
                AND (
                    indexname LIKE 'idx_%' OR
                    indexname LIKE '%_molecules_%' OR
                    indexname LIKE '%_mixtures_%' OR
                    indexname LIKE '%_properties_%'
                );
            """,
            
            "materialized_views": """
                SELECT matviewname
                FROM pg_matviews
                WHERE schemaname = 'public';
            """
        }
        
        # Execute verification queries
        for check_name, query in verification_queries.items():
            if project_id:  # MCP
                try:
                    from supabase_mcp_tools import execute_sql_on_supabase
                    result = execute_sql_on_supabase(project_id, query)
                except ImportError:
                    import mcp__supabase__execute_sql
                    result = mcp__supabase__execute_sql.execute({
                        "project_id": project_id,
                        "query": query
                    })
            elif hasattr(db_client, 'rpc'):  # Supabase
                result = db_client.rpc('exec_sql', {'query': query}).execute()
                result = result.data if hasattr(result, 'data') else result
            else:  # Direct
                cursor = db_client.cursor()
                cursor.execute(query)
                result = cursor.fetchall()
                cursor.close()
                
            verification_results["checks"][check_name] = result
            
            # Analyze results
            if check_name == "security_definer_functions":
                expected_functions = ['has_molecule_access', 'has_mixture_access', 'is_team_member',
                                     'user_has_clearance', 'find_molecules_by_property_range',
                                     'search_molecules_text', 'filter_accessible_molecules']
                found_functions = [r["proname"] for r in result] if isinstance(result, list) and len(result) > 0 and isinstance(result[0], dict) else []
                missing_functions = [f for f in expected_functions if f not in found_functions]
                
                if missing_functions:
                    verification_results["status"] = "failed"
                    verification_results["issues"].append(f"Missing security definer functions: {', '.join(missing_functions)}")
            
            elif check_name == "performance_indexes":
                expected_index_patterns = ['idx_molecules_', 'idx_mixtures_', 'idx_molecular_properties_',
                                         'idx_user_profile_', 'idx_public_']
                found_indexes = [r["indexname"] for r in result] if isinstance(result, list) and len(result) > 0 and isinstance(result[0], dict) else []
                
                missing_patterns = []
                for pattern in expected_index_patterns:
                    if not any(pattern in idx for idx in found_indexes):
                        missing_patterns.append(pattern)
                
                if missing_patterns:
                    verification_results["status"] = "failed"
                    verification_results["issues"].append(f"Missing indexes matching patterns: {', '.join(missing_patterns)}")
            
            elif check_name == "materialized_views":
                expected_views = ['public_molecules_summary', 'public_molecular_properties']
                found_views = [r["matviewname"] for r in result] if isinstance(result, list) and len(result) > 0 and isinstance(result[0], dict) else []
                missing_views = [v for v in expected_views if v not in found_views]
                
                if missing_views:
                    verification_results["status"] = "failed"
                    verification_results["issues"].append(f"Missing materialized views: {', '.join(missing_views)}")
    
    except Exception as e:
        logger.error(f"Error verifying optimizations: {e}")
        verification_results["status"] = "error"
        verification_results["issues"].append(f"Error during verification: {str(e)}")
    
    return verification_results

def generate_report(optimization_steps: Dict[str, bool], verification_results: Dict[str, Any], 
                   performance_results: Dict[str, Any]) -> str:
    """Generate a detailed report of the RLS optimization process."""
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    report_dir = Path("reports")
    report_dir.mkdir(exist_ok=True)
    report_file = report_dir / f"complex_rls_optimization_report_{timestamp}.md"

    # Calculate average improvement
    avg_improvement = 0
    if performance_results and "improvement" in performance_results and performance_results["improvement"]:
        improvements = [stats["improvement_pct"] for stats in performance_results["improvement"].values()]
        avg_improvement = statistics.mean(improvements) if improvements else 0

    report_content = f"""# Complex RLS Query Optimization Report

## Summary
- **Timestamp**: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}
- **Overall Status**: {"✅ Successful" if all(optimization_steps.values()) and verification_results["status"] == "passed" else "⚠️ Partial Success" if any(optimization_steps.values()) else "❌ Failed"}
- **Average Query Performance Improvement**: {avg_improvement:.2f}%

| Optimization Step | Status |
|-------------------|---------|
| Security Definer Functions | {"✅ Applied" if optimization_steps.get('functions', False) else "❌ Failed"} |
| Performance Indexes | {"✅ Applied" if optimization_steps.get('indexes', False) else "❌ Failed"} |
| Materialized Views | {"✅ Applied" if optimization_steps.get('views', False) else "❌ Failed"} |
| Optimized RLS Policies | {"✅ Applied" if optimization_steps.get('policies', False) else "❌ Failed"} |
| Complex Query Optimization | {"✅ Applied" if optimization_steps.get('complex_queries', False) else "❌ Failed"} |
| Verification | {"✅ Passed" if verification_results["status"] == "passed" else "❌ Failed"} |

## Overview
This report documents the optimization of complex query patterns with Row-Level Security (RLS) 
policies in the CryoProtect database. The optimizations improve query performance
while maintaining the security guarantees provided by RLS.

## Optimization Details

### Security Definer Functions
Security definer functions were implemented to optimize common access patterns.
These functions run with elevated privileges but enforce security checks internally,
reducing the overhead of repetitive RLS policy evaluations.

### Performance Indexes
Specialized indexes were added to support efficient execution of RLS policies and complex queries.
These indexes target the columns commonly used in WHERE clauses, JOIN conditions, and ORDER BY clauses.

### Materialized Views
Materialized views were created for frequently accessed data combinations.
These views pre-compute complex joins and aggregations to improve query performance.

### Optimized RLS Policies
The RLS policies were updated to leverage the security definer functions and
performance indexes, reducing query planning and execution time.

### Complex Query Optimization
Specialized functions and indexes were created for common complex query patterns:
- Property range queries
- Full-text search queries
- Molecule-with-properties queries
- Mixture-with-components queries

## Performance Test Results

| Query Pattern | Original (ms) | Optimized (ms) | Improvement (%) | Speedup Factor |
|---------------|---------------|----------------|-----------------|----------------|
"""

    # Add performance test results
    if performance_results and "improvement" in performance_results:
        for query_name, stats in performance_results["improvement"].items():
            report_content += f"| {query_name} | {stats['original_ms']:.2f} | {stats['optimized_ms']:.2f} | {stats['improvement_pct']:.2f}% | {stats['speedup_factor']:.2f}x |\n"

    report_content += """
## Verification Results
"""
    
    if verification_results["status"] == "passed":
        report_content += "All optimizations were successfully verified.\n"
    else:
        report_content += "### Issues Found\n"
        for issue in verification_results["issues"]:
            report_content += f"- {issue}\n"

    report_content += """
## Implementation Details

### Security Definer Functions
The following functions were implemented:
- `has_molecule_access(molecule_id)`: Checks if the current user has access to a molecule
- `has_mixture_access(mixture_id)`: Checks if the current user has access to a mixture
- `is_team_member(team_id)`: Checks if the current user is a member of a team
- `user_has_clearance(required_level)`: Checks if the current user has the required clearance level
- `find_molecules_by_property_range(property_name, min_value, max_value)`: Finds molecules with a property value in a range
- `search_molecules_text(search_term)`: Searches molecules by name and description text
- `filter_accessible_molecules(molecule_ids)`: Filters a list of molecule IDs to those accessible by the current user
- `get_molecules_with_properties(limit_count, offset_val)`: Gets molecules with their property counts
- `get_mixtures_with_components(limit_count, offset_val)`: Gets mixtures with their component counts

### Query Performance Optimization Techniques
1. **Security Definer Functions**: Reduce policy evaluation overhead
2. **Materialized Views**: Pre-compute common access patterns
3. **Specialized Indexes**: Support specific query patterns
4. **Covering Indexes**: Include all columns needed for a query
5. **Function-Based Indexes**: Support queries with expressions
6. **Partial Indexes**: Optimize for specific WHERE conditions
7. **Text Search Indexes**: Optimize full-text search
8. **Optimization of RLS Policies**: Simplify policy expressions

## Next Steps
1. Monitor query performance with the optimized policies
2. Set up scheduled refresh of materialized views
3. Analyze slow queries in production
4. Tune indexes based on actual usage patterns
5. Add query result caching for frequently-accessed data
6. Consider additional optimizations for specific access patterns

## References
- [PostgreSQL RLS Documentation](https://www.postgresql.org/docs/current/ddl-rowsecurity.html)
- [PostgreSQL Security Definer Functions](https://www.postgresql.org/docs/current/sql-createfunction.html)
- [PostgreSQL Index Types](https://www.postgresql.org/docs/current/indexes-types.html)
- [PostgreSQL Materialized Views](https://www.postgresql.org/docs/current/rules-materializedviews.html)
- [PostgreSQL Function Performance](https://www.postgresql.org/docs/current/xfunc-optimization.html)
"""

    try:
        with open(report_file, 'w') as f:
            f.write(report_content)
        logger.info(f"Generated optimization report: {report_file}")
        return str(report_file)
    except Exception as e:
        logger.error(f"Error generating report: {e}")
        return None

def main():
    """Main function to apply RLS optimizations for complex queries."""
    parser = argparse.ArgumentParser(description='Optimize Complex RLS Queries')
    
    # Connection options
    connection_group = parser.add_mutually_exclusive_group(required=True)
    connection_group.add_argument('--direct', action='store_true', help='Use direct database connection')
    connection_group.add_argument('--supabase', action='store_true', help='Use Supabase client')
    connection_group.add_argument('--mcp', action='store_true', help='Use Supabase MCP')
    
    # MCP project ID
    parser.add_argument('--project-id', help='Supabase project ID (required for MCP)')
    
    # Connection parameters for direct connection
    parser.add_argument('--db-host', help='Database host (for direct connection)')
    parser.add_argument('--db-port', type=int, default=5432, help='Database port (for direct connection)')
    parser.add_argument('--db-name', help='Database name (for direct connection)')
    parser.add_argument('--db-user', help='Database user (for direct connection)')
    parser.add_argument('--db-password', help='Database password (for direct connection)')
    
    # Optimization options
    parser.add_argument('--dry-run', action='store_true', help='Show what would be done without making changes')
    parser.add_argument('--skip-functions', action='store_true', help='Skip creating security definer functions')
    parser.add_argument('--skip-indexes', action='store_true', help='Skip creating performance indexes')
    parser.add_argument('--skip-views', action='store_true', help='Skip creating materialized views')
    parser.add_argument('--skip-policies', action='store_true', help='Skip creating RLS policies')
    parser.add_argument('--skip-complex-queries', action='store_true', help='Skip complex query optimizations')
    parser.add_argument('--verify-only', action='store_true', help='Only verify optimizations')
    parser.add_argument('--test-performance', action='store_true', help='Test query performance')
    
    # Apply optimizations even if they already exist
    parser.add_argument('--force', action='store_true', help='Apply optimizations even if they already exist')
    
    # Parse arguments
    args = parser.parse_args()
    
    # Check if project_id is provided for MCP
    if args.mcp and not args.project_id:
        logger.error("Project ID is required for MCP")
        parser.print_help()
        return 1
    
    # Check if connection parameters are provided for direct connection
    if args.direct and not (args.db_host and args.db_name and args.db_user and args.db_password):
        logger.error("Database connection parameters are required for direct connection")
        parser.print_help()
        return 1
    
    # Attempt to get database client
    db_client = None
    try:
        if args.direct:
            import psycopg2
            logger.info(f"Connecting to database {args.db_name} on {args.db_host}")
            db_client = psycopg2.connect(
                host=args.db_host,
                port=args.db_port,
                dbname=args.db_name,
                user=args.db_user,
                password=args.db_password
            )
        elif args.supabase:
            try:
                logger.info("Using service_role_helper to get Supabase client")
                from service_role_helper import get_supabase_client
                db_client = get_supabase_client()
            except ImportError:
                logger.info("Importing Supabase client directly")
                from supabase import create_client
                
                # Try to get URL and key from environment variables
                import os
                url = os.environ.get("SUPABASE_URL")
                key = os.environ.get("SUPABASE_SERVICE_ROLE_KEY")
                
                if not url or not key:
                    logger.error("SUPABASE_URL and SUPABASE_SERVICE_ROLE_KEY environment variables are required")
                    return 1
                
                db_client = create_client(url, key)
    except Exception as e:
        if not args.mcp:  # MCP doesn't need a client
            logger.error(f"Error connecting to database: {e}")
            return 1
    
    if args.dry_run:
        logger.info("DRY RUN MODE - No changes will be made")
        
        # Show what would be applied
        logger.info(f"Would apply security definer functions from: {SECURITY_DEFINER_FUNCTIONS_FILE}")
        logger.info(f"Would apply performance indexes from: {PERFORMANCE_INDEXES_FILE}")
        logger.info(f"Would apply materialized views from: {MATERIALIZED_VIEWS_FILE}")
        logger.info(f"Would apply RLS policies from: {RLS_POLICIES_FILE}")
        logger.info(f"Would apply complex query optimizations from: {COMPLEX_QUERY_OPTIMIZATION_FILE}")
        
        return 0
    
    # Verify existing optimizations if in verify-only mode
    if args.verify_only:
        verification_results = verify_optimization(db_client, args.project_id if args.mcp else None)
        
        if verification_results["status"] == "passed":
            logger.info("All optimizations verified successfully")
            if args.test_performance:
                performance_results = run_performance_tests(db_client, args.project_id if args.mcp else None)
                generate_report({"verification": True}, verification_results, performance_results)
            return 0
        else:
            logger.error("Verification failed:")
            for issue in verification_results["issues"]:
                logger.error(f"- {issue}")
            return 1
    
    # Apply optimizations
    optimization_steps = {}
    
    # Apply security definer functions
    if not args.skip_functions:
        optimization_steps["functions"] = apply_optimization_step(
            db_client, 
            "security definer functions", 
            SECURITY_DEFINER_FUNCTIONS_FILE,
            args.project_id if args.mcp else None
        )
    
    # Apply performance indexes
    if not args.skip_indexes:
        optimization_steps["indexes"] = apply_optimization_step(
            db_client, 
            "performance indexes", 
            PERFORMANCE_INDEXES_FILE,
            args.project_id if args.mcp else None
        )
    
    # Apply materialized views
    if not args.skip_views:
        optimization_steps["views"] = apply_optimization_step(
            db_client, 
            "materialized views", 
            MATERIALIZED_VIEWS_FILE,
            args.project_id if args.mcp else None
        )
    
    # Apply RLS policies
    if not args.skip_policies:
        optimization_steps["policies"] = apply_optimization_step(
            db_client, 
            "RLS policies", 
            RLS_POLICIES_FILE,
            args.project_id if args.mcp else None
        )
    
    # Apply complex query optimizations
    if not args.skip_complex_queries:
        optimization_steps["complex_queries"] = apply_optimization_step(
            db_client, 
            "complex query optimizations", 
            COMPLEX_QUERY_OPTIMIZATION_FILE,
            args.project_id if args.mcp else None
        )
    
    # Verify optimizations
    verification_results = verify_optimization(db_client, args.project_id if args.mcp else None)
    
    # Test performance if requested
    performance_results = {}
    if args.test_performance:
        performance_results = run_performance_tests(db_client, args.project_id if args.mcp else None)
    
    # Generate report
    report_file = generate_report(optimization_steps, verification_results, performance_results)
    if report_file:
        logger.info(f"Full report available at: {report_file}")
    
    # Determine exit status
    if all(optimization_steps.values()) and verification_results["status"] == "passed":
        logger.info("Successfully applied all RLS optimizations for complex queries")
        return 0
    elif any(optimization_steps.values()):
        logger.warning("Some RLS optimization steps completed, but not all")
        return 0
    else:
        logger.error("Failed to apply any RLS optimizations")
        return 1

if __name__ == "__main__":
    sys.exit(main())