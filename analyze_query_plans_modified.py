#!/usr/bin/env python3
"""
CryoProtect v2 - Database Query Plan Analysis (Modified)

This script analyzes the execution plans of common database queries to identify
potential optimization opportunities.
"""

import os
import json
import time
from dotenv import load_dotenv
from supabase import create_client, Client

# Set up logging
import logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler("query_plan_analysis.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

# Load environment variables
load_dotenv()

# Supabase connection
SUPABASE_URL = os.getenv("SUPABASE_URL")
SUPABASE_KEY = os.getenv("SUPABASE_KEY")

if not SUPABASE_URL or not SUPABASE_KEY:
    raise ValueError("SUPABASE_URL and SUPABASE_KEY must be set in .env file")

# Queries to analyze
QUERIES = [
    {
        "name": "Get all molecules with properties",
        "query": "SELECT * FROM molecules LIMIT 100"
    },
    {
        "name": "Get molecule by ID",
        "query": "SELECT * FROM molecules WHERE id = '{{molecule_id}}'"
    },
    {
        "name": "Get all mixtures with components",
        "query": "SELECT * FROM mixtures LIMIT 100"
    },
    {
        "name": "Get mixture by ID",
        "query": "SELECT * FROM mixtures WHERE id = '{{mixture_id}}'"
    },
    {
        "name": "Get predictions for mixture",
        "query": """
            SELECT p.*, pt.name as property_name, cm.name as calculation_method_name
            FROM predictions p
            JOIN property_types pt ON p.property_type_id = pt.id
            JOIN calculation_methods cm ON p.calculation_method_id = cm.id
            WHERE p.mixture_id = '{{mixture_id}}'
            LIMIT 100
        """
    },
    {
        "name": "Get experiments for mixture",
        "query": """
            SELECT e.*, pt.name as property_name
            FROM experiments e
            JOIN property_types pt ON e.property_type_id = pt.id
            WHERE e.mixture_id = '{{mixture_id}}'
        """
    },
    {
        "name": "Compare predictions with experiments",
        "query": """
            SELECT 
                pt.name as property_name,
                p.numeric_value as predicted_value,
                e.numeric_value as experimental_value,
                ABS(p.numeric_value - e.numeric_value) as absolute_difference,
                (ABS(p.numeric_value - e.numeric_value) / NULLIF(ABS(e.numeric_value), 0)) * 100 as percentage_difference,
                p.confidence as prediction_confidence,
                cm.name as calculation_method
            FROM predictions p
            JOIN experiments e ON p.mixture_id = e.mixture_id AND p.property_type_id = e.property_type_id
            JOIN property_types pt ON p.property_type_id = pt.id
            JOIN calculation_methods cm ON p.calculation_method_id = cm.id
            WHERE p.mixture_id = '{{mixture_id}}'
            AND (p.property_type_id = '{{property_type_id}}' OR '{{property_type_id}}' IS NULL)
        """
    },
    {
        "name": "Search molecules by name",
        "query": "SELECT * FROM molecules WHERE name ILIKE '%{{search_term}}%' LIMIT 100"
    },
    {
        "name": "Search mixtures by name",
        "query": "SELECT * FROM mixtures WHERE name ILIKE '%{{search_term}}%' LIMIT 100"
    },
    {
        "name": "Get molecules with specific property value range",
        "query": """
            SELECT m.*, mp.numeric_value
            FROM molecules m
            JOIN molecular_properties mp ON m.id = mp.molecule_id
            JOIN property_types pt ON mp.property_type_id = pt.id
            WHERE pt.name = '{{property_name}}'
            AND mp.numeric_value BETWEEN {{min_value}} AND {{max_value}}
            LIMIT 100
        """
    }
]

def connect_to_supabase():
    """Connect to Supabase using service role key."""
    supabase = create_client(SUPABASE_URL, SUPABASE_KEY)
    logger.info("Connected to Supabase using service role key")
    return supabase

def get_test_data(supabase):
    """Get test data for query parameters."""
    test_data = {}
    
    # Get a molecule ID
    response = supabase.table("molecules").select("id").limit(1).execute()
    if hasattr(response, 'data') and response.data:
        test_data["molecule_id"] = response.data[0]["id"]
    else:
        test_data["molecule_id"] = "00000000-0000-0000-0000-000000000000"  # Fallback
    
    # Get a mixture ID
    response = supabase.table("mixtures").select("id").limit(1).execute()
    if hasattr(response, 'data') and response.data:
        test_data["mixture_id"] = response.data[0]["id"]
    else:
        test_data["mixture_id"] = "00000000-0000-0000-0000-000000000000"  # Fallback
    
    # Get a property type ID
    response = supabase.table("property_types").select("id, name").limit(1).execute()
    if hasattr(response, 'data') and response.data:
        test_data["property_type_id"] = response.data[0]["id"]
        test_data["property_name"] = response.data[0]["name"]
    else:
        test_data["property_type_id"] = "00000000-0000-0000-0000-000000000000"  # Fallback
        test_data["property_name"] = "LogP"  # Fallback
    
    # Set other test data
    test_data["search_term"] = "glyc"
    test_data["min_value"] = -10
    test_data["max_value"] = 10
    
    return test_data

def replace_placeholders(query, test_data):
    """Replace placeholders in the query with test data."""
    for key, value in test_data.items():
        placeholder = "{{" + key + "}}"
        query = query.replace(placeholder, str(value))
    return query

def analyze_query_plan(supabase, query):
    """Get the execution plan for a query."""
    try:
        # Use EXPLAIN ANALYZE to get the execution plan
        explain_query = f"EXPLAIN ANALYZE {query}"
        response = supabase.rpc("exec_sql", {"query": explain_query}).execute()
        
        if hasattr(response, 'error') and response.error:
            logger.error(f"Error analyzing query plan: {response.error}")
            return None
        
        return response.data
    except Exception as e:
        logger.error(f"Error analyzing query plan: {str(e)}")
        return None

def execute_query(supabase, query):
    """Execute a query and measure its performance."""
    try:
        start_time = time.time()
        response = supabase.rpc("exec_sql", {"query": query}).execute()
        end_time = time.time()
        
        if hasattr(response, 'error') and response.error:
            logger.error(f"Error executing query: {response.error}")
            return None, None
        
        execution_time = (end_time - start_time) * 1000  # Convert to milliseconds
        return response.data, execution_time
    except Exception as e:
        logger.error(f"Error executing query: {str(e)}")
        return None, None

def analyze_plan_for_issues(plan):
    """Analyze the execution plan for potential issues."""
    issues = []
    
    if not plan:
        return issues
    
    plan_text = "\n".join(plan)
    
    # Check for sequential scans
    if "Seq Scan" in plan_text:
        issues.append({
            "type": "sequential_scan",
            "description": "Sequential scan detected. Consider adding an index to improve performance.",
            "severity": "high"
        })
    
    # Check for hash joins (might indicate missing indexes)
    if "Hash Join" in plan_text:
        issues.append({
            "type": "hash_join",
            "description": "Hash join detected. Consider adding indexes on join columns.",
            "severity": "medium"
        })
    
    # Check for high cost operations
    if "cost=" in plan_text:
        # Extract cost values
        import re
        cost_matches = re.findall(r"cost=(\d+\.\d+)\.\.(\d+\.\d+)", plan_text)
        if cost_matches:
            max_cost = max([float(match[1]) for match in cost_matches])
            if max_cost > 1000:
                issues.append({
                    "type": "high_cost",
                    "description": f"High cost operation detected (cost={max_cost}). Consider optimizing the query.",
                    "severity": "high"
                })
    
    # Check for sorting operations
    if "Sort" in plan_text:
        issues.append({
            "type": "sort",
            "description": "Sort operation detected. Consider adding an index to avoid sorting.",
            "severity": "medium"
        })
    
    # Check for temporary files
    if "temporary file" in plan_text.lower():
        issues.append({
            "type": "temporary_file",
            "description": "Query is using temporary files. Consider increasing work_mem or optimizing the query.",
            "severity": "high"
        })
    
    return issues

def suggest_optimizations(query, plan, issues):
    """Suggest optimizations based on the query and execution plan."""
    suggestions = []
    
    # Check for table scans and suggest indexes
    for issue in issues:
        if issue["type"] == "sequential_scan":
            # Extract table name from the plan
            import re
            table_matches = re.findall(r"Seq Scan on (\w+)", "\n".join(plan))
            for table in table_matches:
                # Extract WHERE conditions
                where_matches = re.findall(r"WHERE\s+([^()]+)(?:\s+AND|\s+OR|\s*$)", query, re.IGNORECASE)
                if where_matches:
                    for where_clause in where_matches:
                        # Extract column names from WHERE clause
                        column_matches = re.findall(r"(\w+)\s*[=><]", where_clause)
                        for column in column_matches:
                            suggestions.append({
                                "type": "create_index",
                                "description": f"Create an index on {table}.{column} to improve query performance.",
                                "sql": f"CREATE INDEX idx_{table}_{column} ON {table} ({column});"
                            })
    
    # Check for joins and suggest indexes
    if "JOIN" in query.upper():
        # Extract join conditions
        import re
        join_matches = re.findall(r"(\w+)\s+(?:INNER|LEFT|RIGHT|FULL)?\s*JOIN\s+(\w+)\s+ON\s+([^()]+)(?:\s+AND|\s+OR|\s*$)", query, re.IGNORECASE)
        for match in join_matches:
            table1, table2, join_condition = match
            # Extract column names from join condition
            column_matches = re.findall(r"(\w+\.\w+)\s*=\s*(\w+\.\w+)", join_condition)
            for columns in column_matches:
                col1, col2 = columns
                table1_col = col1.split(".")[1]
                table2_col = col2.split(".")[1]
                suggestions.append({
                    "type": "create_index",
                    "description": f"Create indexes on join columns {col1} and {col2} to improve join performance.",
                    "sql": f"CREATE INDEX idx_{table1}_{table1_col} ON {table1} ({table1_col});\nCREATE INDEX idx_{table2}_{table2_col} ON {table2} ({table2_col});"
                })
    
    # Check for LIKE/ILIKE and suggest text indexes
    if "LIKE" in query.upper() or "ILIKE" in query.upper():
        # Extract LIKE conditions
        import re
        like_matches = re.findall(r"(\w+)\s+(?:I)?LIKE\s+'%([^%]+)%'", query, re.IGNORECASE)
        for match in like_matches:
            column, _ = match
            # Extract table name from the query
            from_matches = re.findall(r"FROM\s+(\w+)", query, re.IGNORECASE)
            if from_matches:
                table = from_matches[0]
                suggestions.append({
                    "type": "create_text_index",
                    "description": f"Create a GIN index with pg_trgm extension on {table}.{column} to improve text search performance.",
                    "sql": f"CREATE EXTENSION IF NOT EXISTS pg_trgm;\nCREATE INDEX idx_{table}_{column}_trgm ON {table} USING gin ({column} gin_trgm_ops);"
                })
    
    return suggestions

def main():
    """Main function."""
    logger.info("Starting database query plan analysis...")
    
    # Connect to Supabase
    supabase = connect_to_supabase()
    
    # Get test data
    test_data = get_test_data(supabase)
    logger.info(f"Using test data: {test_data}")
    
    # Analyze queries
    results = []
    
    for query_info in QUERIES:
        logger.info(f"Analyzing query: {query_info['name']}")
        
        # Replace placeholders
        query = replace_placeholders(query_info["query"], test_data)
        
        # Get execution plan
        plan = analyze_query_plan(supabase, query)
        
        # Execute query and measure performance
        result, execution_time = execute_query(supabase, query)
        
        # Analyze plan for issues
        issues = analyze_plan_for_issues(plan)
        
        # Suggest optimizations
        suggestions = suggest_optimizations(query, plan, issues)
        
        # Record results
        results.append({
            "name": query_info["name"],
            "query": query,
            "plan": plan,
            "execution_time_ms": execution_time,
            "row_count": len(result) if result else 0,
            "issues": issues,
            "suggestions": suggestions
        })
    
    # Generate report
    report = {
        "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
        "queries": results
    }
    
    # Create reports directory if it doesn't exist
    os.makedirs("reports", exist_ok=True)
    
    # Save report to file
    report_path = "reports/query_plan_analysis_report.json"
    with open(report_path, "w") as f:
        json.dump(report, f, indent=2)
    
    # Generate text report
    text_report = []
    
    text_report.append("=" * 80)
    text_report.append("CryoProtect v2 Database Query Plan Analysis Report")
    text_report.append("=" * 80)
    text_report.append(f"Timestamp: {report['timestamp']}")
    text_report.append("")
    
    for query_result in report["queries"]:
        text_report.append("-" * 80)
        text_report.append(f"Query: {query_result['name']}")
        text_report.append(f"Execution Time: {query_result['execution_time_ms']:.2f} ms")
        text_report.append(f"Row Count: {query_result['row_count']}")
        text_report.append("")
        text_report.append("SQL:")
        text_report.append(query_result["query"])
        text_report.append("")
        
        if query_result["issues"]:
            text_report.append("Issues:")
            for issue in query_result["issues"]:
                text_report.append(f"  - [{issue['severity'].upper()}] {issue['description']}")
            text_report.append("")
        
        if query_result["suggestions"]:
            text_report.append("Optimization Suggestions:")
            for suggestion in query_result["suggestions"]:
                text_report.append(f"  - {suggestion['description']}")
                text_report.append(f"    SQL: {suggestion['sql']}")
            text_report.append("")
        
        if query_result["plan"]:
            text_report.append("Execution Plan:")
            for line in query_result["plan"]:
                text_report.append(f"  {line}")
            text_report.append("")
    
    text_report.append("=" * 80)
    text_report.append("Summary of Optimization Suggestions:")
    text_report.append("=" * 80)
    
    # Collect unique suggestions
    unique_suggestions = {}
    for query_result in report["queries"]:
        for suggestion in query_result["suggestions"]:
            if suggestion["sql"] not in unique_suggestions:
                unique_suggestions[suggestion["sql"]] = suggestion["description"]
    
    for sql, description in unique_suggestions.items():
        text_report.append(f"- {description}")
        text_report.append(f"  SQL: {sql}")
        text_report.append("")
    
    text_report.append("=" * 80)
    
    # Save text report to file
    text_report_path = "reports/query_plan_analysis_report.txt"
    with open(text_report_path, "w") as f:
        f.write("\n".join(text_report))
    
    logger.info(f"Query plan analysis completed. Reports saved to {report_path} and {text_report_path}")
    
    return 0

if __name__ == "__main__":
    exit(main())