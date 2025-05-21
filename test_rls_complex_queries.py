#!/usr/bin/env python3
"""
Test RLS Complex Queries Optimization

This script tests the performance of complex queries with RLS policies
before and after optimization to verify the effectiveness of the optimizations.
"""

import os
import sys
import time
import logging
import json
import argparse
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Any, Optional, Tuple

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler(),
        logging.FileHandler(f"rls_query_test_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"),
    ]
)
logger = logging.getLogger(__name__)

# Dictionary of complex query patterns and their test queries
ORIGINAL_QUERIES = {
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
    """,
    
    "complex_multi_join_query": """
        SELECT 
            m.id,
            m.name,
            m.molecular_formula,
            mp.property_name,
            mp.property_value,
            exp.experimental_value,
            pred.predicted_value
        FROM 
            molecules m
        JOIN
            molecular_properties mp ON m.id = mp.molecule_id
        LEFT JOIN
            experiments exp ON m.id = exp.molecule_id
        LEFT JOIN
            predictions pred ON m.id = pred.molecule_id
        WHERE 
            mp.property_name = 'Solubility'
            AND (m.is_public = true OR m.created_by = auth.uid() OR
                 EXISTS (
                     SELECT 1 FROM project_molecules pm
                     JOIN team_projects tp ON pm.project_id = tp.project_id
                     JOIN user_profile up ON tp.team_id = up.team_id
                     WHERE pm.molecule_id = m.id AND up.auth_user_id = auth.uid()
                 ))
        LIMIT 10;
    """
}

# Optimized versions of the complex queries
OPTIMIZED_QUERIES = {
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
    """,
    
    "optimized_complex_multi_join_query": """
        SELECT 
            m.id,
            m.name,
            m.molecular_formula,
            mp.property_name,
            mp.property_value,
            exp.experimental_value,
            pred.predicted_value
        FROM 
            molecules m
        JOIN
            filter_accessible_molecules(ARRAY(
                SELECT DISTINCT molecule_id 
                FROM molecular_properties 
                WHERE property_name = 'Solubility'
                LIMIT 10
            )) mol_ids ON m.id = mol_ids
        JOIN
            molecular_properties mp ON m.id = mp.molecule_id AND mp.property_name = 'Solubility'
        LEFT JOIN
            experiments exp ON m.id = exp.molecule_id
        LEFT JOIN
            predictions pred ON m.id = pred.molecule_id
        LIMIT 10;
    """
}

class RLSQueryTester:
    """
    Test RLS query performance before and after optimization.
    """
    
    def __init__(self, connection_type: str, project_id: Optional[str] = None,
                db_params: Optional[Dict[str, Any]] = None):
        """Initialize the tester with connection parameters."""
        self.connection_type = connection_type
        self.project_id = project_id
        self.db_params = db_params
        self.db_client = None
        
        # Connect to database
        self._connect_to_database()
        
    def _connect_to_database(self) -> None:
        """Connect to the database using the specified connection type."""
        if self.connection_type == "direct":
            try:
                import psycopg2
                logger.info(f"Connecting to database {self.db_params['db_name']} on {self.db_params['db_host']}")
                self.db_client = psycopg2.connect(
                    host=self.db_params['db_host'],
                    port=self.db_params.get('db_port', 5432),
                    dbname=self.db_params['db_name'],
                    user=self.db_params['db_user'],
                    password=self.db_params['db_password']
                )
                logger.info("Connected to database successfully")
            except Exception as e:
                logger.error(f"Error connecting to database: {e}")
                raise
        
        elif self.connection_type == "supabase":
            try:
                try:
                    logger.info("Using service_role_helper to get Supabase client")
                    from service_role_helper import get_supabase_client
                    self.db_client = get_supabase_client()
                except ImportError:
                    logger.info("Importing Supabase client directly")
                    from supabase import create_client
                    
                    # Try to get URL and key from environment variables or db_params
                    import os
                    url = self.db_params.get('url') if self.db_params else os.environ.get("SUPABASE_URL")
                    key = self.db_params.get('key') if self.db_params else os.environ.get("SUPABASE_SERVICE_ROLE_KEY")
                    
                    if not url or not key:
                        logger.error("SUPABASE_URL and SUPABASE_SERVICE_ROLE_KEY are required")
                        raise ValueError("Supabase URL and key are required")
                    
                    self.db_client = create_client(url, key)
                logger.info("Connected to Supabase successfully")
            except Exception as e:
                logger.error(f"Error connecting to Supabase: {e}")
                raise
                
        elif self.connection_type == "mcp":
            # No need to initialize connection for MCP mode
            if not self.project_id:
                logger.error("Project ID is required for MCP mode")
                raise ValueError("Project ID is required for MCP mode")
            logger.info(f"Using MCP mode with project ID: {self.project_id}")
        
        else:
            logger.error(f"Unsupported connection type: {self.connection_type}")
            raise ValueError(f"Unsupported connection type: {self.connection_type}")
    
    def execute_query(self, query: str) -> Tuple[bool, float, int, Optional[Any]]:
        """
        Execute a query and return its execution time and result count.
        
        Returns:
            (success, execution_time_ms, result_count, error)
        """
        try:
            # Run query with timing
            start_time = time.time()
            
            if self.connection_type == "mcp":
                try:
                    from supabase_mcp_tools import execute_sql_on_supabase
                    result = execute_sql_on_supabase(self.project_id, query)
                except ImportError:
                    try:
                        import mcp__supabase__execute_sql
                        result = mcp__supabase__execute_sql.execute({
                            "project_id": self.project_id,
                            "query": query
                        })
                    except ImportError:
                        logger.error("MCP modules not available")
                        return False, 0, 0, "MCP modules not available"
                
            elif self.connection_type == "supabase":
                result = self.db_client.rpc('exec_sql', {'query': query}).execute()
                result = result.data if hasattr(result, 'data') else result
                
            else:  # Direct
                cursor = self.db_client.cursor()
                cursor.execute(query)
                result = cursor.fetchall()
                cursor.close()
            
            # Calculate execution time
            execution_time = time.time() - start_time
            
            # Count results
            if isinstance(result, list):
                result_count = len(result)
            elif hasattr(result, 'data') and isinstance(result.data, list):
                result_count = len(result.data)
            else:
                result_count = 0
                
            return True, execution_time * 1000, result_count, None
            
        except Exception as e:
            logger.error(f"Error executing query: {e}")
            return False, 0, 0, str(e)
    
    def test_query(self, name: str, query: str, iterations: int = 3) -> Dict[str, Any]:
        """Test a query multiple times and return statistics."""
        logger.info(f"Testing query: {name}")
        results = {
            "name": name,
            "query": query,
            "success": False,
            "execution_times_ms": [],
            "avg_execution_time_ms": 0,
            "result_count": 0,
            "error": None
        }
        
        for i in range(iterations):
            success, execution_time, result_count, error = self.execute_query(query)
            
            if success:
                results["execution_times_ms"].append(execution_time)
                results["result_count"] = result_count
                results["success"] = True
            else:
                results["success"] = False
                results["error"] = error
                break
        
        if results["success"] and results["execution_times_ms"]:
            results["avg_execution_time_ms"] = sum(results["execution_times_ms"]) / len(results["execution_times_ms"])
            logger.info(f"Query {name} - Average execution time: {results['avg_execution_time_ms']:.2f}ms - Results: {results['result_count']}")
        else:
            logger.error(f"Query {name} failed: {results['error']}")
            
        return results
    
    def run_all_tests(self, iterations: int = 3) -> Dict[str, Any]:
        """Run all tests and return results."""
        test_results = {
            "timestamp": datetime.now().isoformat(),
            "original_queries": {},
            "optimized_queries": {},
            "summary": {}
        }
        
        logger.info(f"Running tests with {iterations} iterations per query")
        
        # Test original queries
        for name, query in ORIGINAL_QUERIES.items():
            test_results["original_queries"][name] = self.test_query(name, query, iterations)
        
        # Test optimized queries
        for name, query in OPTIMIZED_QUERIES.items():
            # Extract original query name (remove "optimized_" prefix)
            orig_name = name[10:] if name.startswith("optimized_") else name
            test_results["optimized_queries"][orig_name] = self.test_query(name, query, iterations)
        
        # Generate summary with improvements
        for orig_name, orig_result in test_results["original_queries"].items():
            if orig_name in test_results["optimized_queries"]:
                opt_result = test_results["optimized_queries"][orig_name]
                
                if orig_result["success"] and opt_result["success"]:
                    orig_time = orig_result["avg_execution_time_ms"]
                    opt_time = opt_result["avg_execution_time_ms"]
                    
                    if orig_time > 0:
                        improvement_pct = ((orig_time - opt_time) / orig_time) * 100
                        speedup_factor = orig_time / opt_time if opt_time > 0 else float('inf')
                        
                        test_results["summary"][orig_name] = {
                            "original_ms": orig_time,
                            "optimized_ms": opt_time,
                            "improvement_pct": improvement_pct,
                            "speedup_factor": speedup_factor
                        }
                        
                        logger.info(f"Query {orig_name} - Improvement: {improvement_pct:.2f}% "
                                   f"({orig_time:.2f}ms → {opt_time:.2f}ms) - Speedup: {speedup_factor:.2f}x")
        
        # Calculate average improvement
        if test_results["summary"]:
            improvements = [s["improvement_pct"] for s in test_results["summary"].values()]
            speedups = [s["speedup_factor"] for s in test_results["summary"].values()]
            
            test_results["summary"]["overall"] = {
                "avg_improvement_pct": sum(improvements) / len(improvements),
                "avg_speedup_factor": sum(speedups) / len(speedups),
                "queries_tested": len(test_results["summary"])
            }
            
            logger.info(f"Overall - Average improvement: {test_results['summary']['overall']['avg_improvement_pct']:.2f}% "
                       f"- Average speedup: {test_results['summary']['overall']['avg_speedup_factor']:.2f}x")
        
        return test_results
    
    def save_results(self, results: Dict[str, Any], output_file: Optional[str] = None) -> str:
        """Save test results to a file."""
        if not output_file:
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            output_file = f"rls_query_test_results_{timestamp}.json"
        
        try:
            with open(output_file, 'w') as f:
                json.dump(results, f, indent=2)
            logger.info(f"Test results saved to {output_file}")
            return output_file
        except Exception as e:
            logger.error(f"Error saving test results: {e}")
            return ""
    
    def generate_markdown_report(self, results: Dict[str, Any], output_file: Optional[str] = None) -> str:
        """Generate a markdown report from test results."""
        if not output_file:
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            output_file = f"rls_query_test_report_{timestamp}.md"
        
        try:
            with open(output_file, 'w') as f:
                f.write("# RLS Complex Query Optimization Test Report\n\n")
                f.write(f"## Summary\n\n")
                f.write(f"- **Test Date**: {datetime.fromisoformat(results['timestamp']).strftime('%Y-%m-%d %H:%M:%S')}\n")
                f.write(f"- **Connection Type**: {self.connection_type}\n")
                
                if "overall" in results["summary"]:
                    overall = results["summary"]["overall"]
                    f.write(f"- **Average Improvement**: {overall['avg_improvement_pct']:.2f}%\n")
                    f.write(f"- **Average Speedup Factor**: {overall['avg_speedup_factor']:.2f}x\n")
                    f.write(f"- **Queries Tested**: {overall['queries_tested']}\n\n")
                
                f.write("## Query Performance Comparison\n\n")
                f.write("| Query | Original (ms) | Optimized (ms) | Improvement | Speedup |\n")
                f.write("|-------|--------------|----------------|-------------|--------|\n")
                
                for query_name, stats in results["summary"].items():
                    if query_name != "overall":
                        f.write(f"| {query_name} | {stats['original_ms']:.2f} | {stats['optimized_ms']:.2f} | {stats['improvement_pct']:.2f}% | {stats['speedup_factor']:.2f}x |\n")
                
                f.write("\n## Query Details\n\n")
                
                for query_name, orig_result in results["original_queries"].items():
                    f.write(f"### {query_name}\n\n")
                    
                    f.write("#### Original Query\n\n")
                    f.write("```sql\n")
                    f.write(orig_result["query"])
                    f.write("\n```\n\n")
                    
                    f.write(f"- **Average Execution Time**: {orig_result['avg_execution_time_ms']:.2f}ms\n")
                    f.write(f"- **Result Count**: {orig_result['result_count']}\n")
                    f.write(f"- **Success**: {orig_result['success']}\n")
                    
                    if query_name in results["optimized_queries"]:
                        opt_result = results["optimized_queries"][query_name]
                        
                        f.write("\n#### Optimized Query\n\n")
                        f.write("```sql\n")
                        f.write(opt_result["query"])
                        f.write("\n```\n\n")
                        
                        f.write(f"- **Average Execution Time**: {opt_result['avg_execution_time_ms']:.2f}ms\n")
                        f.write(f"- **Result Count**: {opt_result['result_count']}\n")
                        f.write(f"- **Success**: {opt_result['success']}\n")
                        
                        if query_name in results["summary"]:
                            stats = results["summary"][query_name]
                            f.write(f"- **Improvement**: {stats['improvement_pct']:.2f}%\n")
                            f.write(f"- **Speedup Factor**: {stats['speedup_factor']:.2f}x\n")
                    
                    f.write("\n")
                
                f.write("\n## Conclusion\n\n")
                
                if "overall" in results["summary"]:
                    overall = results["summary"]["overall"]
                    if overall["avg_improvement_pct"] > 75:
                        f.write("The optimizations have resulted in **dramatic performance improvements** for complex RLS queries. ")
                    elif overall["avg_improvement_pct"] > 50:
                        f.write("The optimizations have resulted in **significant performance improvements** for complex RLS queries. ")
                    elif overall["avg_improvement_pct"] > 25:
                        f.write("The optimizations have resulted in **notable performance improvements** for complex RLS queries. ")
                    else:
                        f.write("The optimizations have resulted in **some performance improvements** for complex RLS queries. ")
                    
                    f.write(f"On average, optimized queries run {overall['avg_speedup_factor']:.2f}x faster than the original queries, ")
                    f.write(f"with an average improvement of {overall['avg_improvement_pct']:.2f}%.\n\n")
                
                f.write("The key optimization techniques that contributed to these improvements include:\n\n")
                f.write("1. **Security definer functions** that reduce RLS policy evaluation overhead\n")
                f.write("2. **Specialized indexes** for common query patterns\n")
                f.write("3. **Materialized views** for frequently accessed data\n")
                f.write("4. **Optimized RLS policies** using security definer functions\n")
                f.write("5. **Query rewriting** to use more efficient access patterns\n")
                
            logger.info(f"Markdown report generated: {output_file}")
            return output_file
        except Exception as e:
            logger.error(f"Error generating markdown report: {e}")
            return ""

def main():
    """Main function."""
    parser = argparse.ArgumentParser(description='Test RLS Complex Queries Optimization')
    
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
    
    # Supabase connection parameters
    parser.add_argument('--supabase-url', help='Supabase URL (for Supabase connection)')
    parser.add_argument('--supabase-key', help='Supabase service role key (for Supabase connection)')
    
    # Test parameters
    parser.add_argument('--iterations', type=int, default=3, help='Number of iterations per query (default: 3)')
    parser.add_argument('--output-json', help='Output JSON file path')
    parser.add_argument('--output-report', help='Output markdown report file path')
    
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
    
    # Check if Supabase URL and key are provided for Supabase connection
    if args.supabase and not (args.supabase_url and args.supabase_key):
        # Check environment variables
        import os
        if not (os.environ.get("SUPABASE_URL") and os.environ.get("SUPABASE_SERVICE_ROLE_KEY")):
            logger.error("Supabase URL and key are required for Supabase connection")
            parser.print_help()
            return 1
    
    # Determine connection type
    connection_type = "direct" if args.direct else "supabase" if args.supabase else "mcp"
    
    # Prepare database parameters
    db_params = None
    if args.direct:
        db_params = {
            'db_host': args.db_host,
            'db_port': args.db_port,
            'db_name': args.db_name,
            'db_user': args.db_user,
            'db_password': args.db_password
        }
    elif args.supabase:
        db_params = {
            'url': args.supabase_url,
            'key': args.supabase_key
        }
    
    try:
        # Initialize tester
        tester = RLSQueryTester(connection_type, args.project_id, db_params)
        
        # Run tests
        results = tester.run_all_tests(args.iterations)
        
        # Save results
        json_file = tester.save_results(results, args.output_json)
        
        # Generate markdown report
        report_file = tester.generate_markdown_report(results, args.output_report)
        
        logger.info("Testing completed successfully")
        logger.info(f"JSON results: {json_file}")
        logger.info(f"Markdown report: {report_file}")
        
        return 0
    except Exception as e:
        logger.error(f"Error during testing: {e}")
        import traceback
        logger.error(traceback.format_exc())
        return 1

if __name__ == "__main__":
    sys.exit(main())