"""
Database Health Check Module for CryoProtect v2

This module provides comprehensive database health checking functionality through
a class-based orchestrator that integrates specialized validator modules:
- SchemaValidator: Validates table structure, column types, constraints
- IntegrityChecker: Checks data integrity, orphaned records, invalid data
- PerformanceAnalyzer: Diagnoses performance issues, slow queries, missing indexes

The DatabaseHealthCheck class serves as the main entry point for running health checks
and generating comprehensive reports in various formats.
"""

import logging
import json
import os
import time
import traceback
from datetime import datetime, timezone
from typing import Dict, List, Any, Optional, Tuple, Set, Union

# Import validator modules
from database.utils.schema_validator import SchemaValidator
from database.utils.integrity_checker import IntegrityChecker
from database.utils.performance_analyzer import PerformanceAnalyzer

# Setup logging
logger = logging.getLogger(__name__)

# Constants
HEALTH_CHECK_REPORT_DIR = os.path.join('reports', 'database')
os.makedirs(HEALTH_CHECK_REPORT_DIR, exist_ok=True)


class DatabaseHealthCheck:
    """
    Orchestrator for the database health check system.
    
    This class integrates specialized validator modules to provide a comprehensive
    health check system for the database. It manages the execution of various checks,
    aggregates results, calculates overall health status, and generates reports in
    different formats.
    
    Attributes:
        connection_pool: Database connection pool
        config: Configuration dictionary for health checks
        validators: Dictionary of registered validator modules
    """
    
    def __init__(self, connection_pool, config=None):
        """
        Initialize the database health check orchestrator.
        
        Args:
            connection_pool: Database connection pool
            config: Optional configuration dictionary for health checks
        """
        self.connection_pool = connection_pool
        self.config = config or {}
        
        # Register validator modules
        self.validators = {
            'schema': SchemaValidator(),
            'integrity': IntegrityChecker(),
            'performance': PerformanceAnalyzer()
        }
        
        logger.info("DatabaseHealthCheck initialized with %d validator modules", len(self.validators))
    
    def run_health_check(self, categories=None) -> Dict[str, Any]:
        """
        Run health checks for selected or all categories.
        
        Args:
            categories: Optional list of categories to check ('schema', 'integrity', 'performance')
                       If None, all checks will be run
        
        Returns:
            Dict with aggregated health check results
        """
        start_time = time.time()
        
        # Initialize results structure
        results = {
            'timestamp': datetime.now(timezone.utc).isoformat(),
            'overall_status': 'passed',
            'categories_checked': [],
            'total_issues_found': 0,
            'results_by_category': {},
            'recommendations': []
        }
        
        try:
            # Determine which categories to check
            categories_to_check = categories or list(self.validators.keys())
            results['categories_checked'] = categories_to_check
            
            # Run checks for each selected category
            for category in categories_to_check:
                if category in self.validators:
                    logger.info("Running %s health checks", category)
                    validator = self.validators[category]
                    
                    # Run the checks for this validator
                    category_results = validator.run_checks(self.connection_pool)
                    
                    # Store results
                    results['results_by_category'][category] = category_results
                    
                    # Update total issues count
                    results['total_issues_found'] += category_results.get('issues_found', 0)
                    
                    # Collect recommendations
                    if 'recommendations' in category_results:
                        results['recommendations'].extend(category_results['recommendations'])
                else:
                    logger.warning("Unknown health check category: %s", category)
            
            # Calculate overall status
            results['overall_status'] = self._calculate_status(results)
            
        except Exception as e:
            logger.error("Error during health check execution: %s", str(e), exc_info=True)
            results['overall_status'] = 'error'
            results['error'] = str(e)
            results['traceback'] = traceback.format_exc()
        
        # Add execution time
        execution_time = time.time() - start_time
        results['execution_time_seconds'] = execution_time
        logger.info("Health check completed in %.2f seconds with status: %s", 
                   execution_time, results['overall_status'])
        
        return results
    
    def _calculate_status(self, results: Dict[str, Any]) -> str:
        """
        Calculate overall health status based on category results.
        
        Args:
            results: Health check results dictionary
            
        Returns:
            Overall status string ('passed', 'warning', 'failed', or 'error')
        """
        # Start with the assumption that everything passed
        overall_status = 'passed'
        
        # Check if any category has error status
        for category, category_results in results['results_by_category'].items():
            category_status = category_results.get('status', 'passed')
            
            # Error is the worst status
            if category_status == 'error':
                return 'error'
            
            # Failed is the next worst status
            elif category_status == 'failed' and overall_status != 'error':
                overall_status = 'failed'
            
            # Warning is the least severe problematic status
            elif category_status == 'warning' and overall_status not in ('error', 'failed'):
                overall_status = 'warning'
        
        return overall_status
    
    def generate_report(self, health_results=None, format="json", output_file=None) -> Union[str, Dict[str, Any]]:
        """
        Generate a health check report in the specified format.
        
        Args:
            health_results: Health check results to include in the report
                           If None, a new health check will be run
            format: Report format ('json', 'markdown', or 'html')
            output_file: Optional file path to save the report
            
        Returns:
            Report content as string (for markdown/html) or dict (for json)
        """
        # Run health check if results not provided
        if health_results is None:
            health_results = self.run_health_check()
        
        # Generate report in the specified format
        if format.lower() == 'json':
            report = health_results
            report_content = json.dumps(report, indent=2)
        
        elif format.lower() == 'markdown':
            report_content = self._generate_markdown_report(health_results)
        
        elif format.lower() == 'html':
            report_content = self._generate_html_report(health_results)
        
        else:
            logger.warning("Unsupported report format: %s. Defaulting to JSON.", format)
            report = health_results
            report_content = json.dumps(report, indent=2)
            format = 'json'
        
        # Save report to file if output_file is specified
        if output_file:
            # Create directory if it doesn't exist
            os.makedirs(os.path.dirname(os.path.abspath(output_file)), exist_ok=True)
            
            with open(output_file, 'w') as f:
                f.write(report_content)
            
            logger.info("Health check report saved to %s", output_file)
        else:
            # If no output file specified, save to default location with timestamp
            timestamp_str = datetime.now().strftime('%Y%m%d_%H%M%S')
            default_filename = f'health_report_{timestamp_str}.{format}'
            default_path = os.path.join(HEALTH_CHECK_REPORT_DIR, default_filename)
            
            with open(default_path, 'w') as f:
                f.write(report_content)
            
            logger.info("Health check report saved to %s", default_path)
        
        # Return the report content
        if format.lower() == 'json':
            return health_results
        else:
            return report_content
    
    def _generate_markdown_report(self, results: Dict[str, Any]) -> str:
        """
        Generate a markdown report from health check results.
        
        Args:
            results: Health check results dictionary
            
        Returns:
            Markdown formatted report as string
        """
        md = []
        
        # Header
        md.append("# Database Health Check Report\n")
        md.append(f"**Generated:** {results['timestamp']}\n")
        md.append(f"**Overall Status:** {results['overall_status'].upper()}\n")
        md.append(f"**Total Issues Found:** {results['total_issues_found']}\n")
        md.append(f"**Categories Checked:** {', '.join(results['categories_checked'])}\n")
        md.append(f"**Execution Time:** {results.get('execution_time_seconds', 0):.2f} seconds\n")
        
        # Summary of issues by category
        md.append("\n## Summary by Category\n")
        for category, category_results in results['results_by_category'].items():
            status = category_results.get('status', 'unknown')
            issues = category_results.get('issues_found', 0)
            md.append(f"- **{category.capitalize()}**: {status.upper()} ({issues} issues)\n")
        
        # Recommendations
        if results['recommendations']:
            md.append("\n## Recommendations\n")
            for i, recommendation in enumerate(results['recommendations'], 1):
                md.append(f"{i}. {recommendation}\n")
        
        # Detailed results by category
        md.append("\n## Detailed Results\n")
        for category, category_results in results['results_by_category'].items():
            md.append(f"\n### {category.capitalize()} Checks\n")
            md.append(f"**Status:** {category_results.get('status', 'unknown').upper()}\n")
            md.append(f"**Issues Found:** {category_results.get('issues_found', 0)}\n")
            
            # Add category-specific details
            details = category_results.get('details', {})
            if details:
                md.append("\n#### Details\n")
                for check_name, check_details in details.items():
                    md.append(f"\n##### {check_name.replace('_', ' ').title()}\n")
                    md.append(f"Status: {check_details.get('status', 'unknown')}\n")
                    
                    # Add more specific details based on the check type
                    if isinstance(check_details, dict):
                        for key, value in check_details.items():
                            if key not in ('status', 'recommendations'):
                                md.append(f"- {key.replace('_', ' ').title()}: ")
                                if isinstance(value, list):
                                    if value:
                                        md.append("\n")
                                        for item in value[:5]:  # Limit to first 5 items
                                            md.append(f"  - {str(item)}\n")
                                        if len(value) > 5:
                                            md.append(f"  - ... and {len(value) - 5} more\n")
                                    else:
                                        md.append("None\n")
                                else:
                                    md.append(f"{str(value)}\n")
        
        # Error information if present
        if 'error' in results:
            md.append("\n## Errors\n")
            md.append(f"```\n{results['error']}\n```\n")
            
            if 'traceback' in results:
                md.append("\n### Traceback\n")
                md.append(f"```\n{results['traceback']}\n```\n")
        
        return ''.join(md)
    
    def _generate_html_report(self, results: Dict[str, Any]) -> str:
        """
        Generate an HTML report from health check results.
        
        Args:
            results: Health check results dictionary
            
        Returns:
            HTML formatted report as string
        """
        # Convert markdown to HTML
        markdown_report = self._generate_markdown_report(results)
        
        # Simple HTML wrapper with some basic styling
        html = [
            "<!DOCTYPE html>",
            "<html>",
            "<head>",
            "    <title>Database Health Check Report</title>",
            "    <style>",
            "        body { font-family: Arial, sans-serif; margin: 40px; line-height: 1.6; }",
            "        h1 { color: #2c3e50; }",
            "        h2 { color: #3498db; border-bottom: 1px solid #eee; padding-bottom: 10px; }",
            "        h3 { color: #2980b9; }",
            "        h4 { color: #27ae60; }",
            "        h5 { color: #16a085; }",
            "        .passed { color: green; }",
            "        .warning { color: orange; }",
            "        .failed { color: red; }",
            "        .error { color: darkred; }",
            "        pre { background-color: #f8f8f8; padding: 10px; border-radius: 5px; overflow-x: auto; }",
            "        table { border-collapse: collapse; width: 100%; }",
            "        th, td { border: 1px solid #ddd; padding: 8px; text-align: left; }",
            "        th { background-color: #f2f2f2; }",
            "    </style>",
            "</head>",
            "<body>"
        ]
        
        # Convert markdown headers to HTML
        markdown_lines = markdown_report.split('\n')
        for line in markdown_lines:
            if line.startswith('# '):
                html.append(f"<h1>{line[2:]}</h1>")
            elif line.startswith('## '):
                html.append(f"<h2>{line[3:]}</h2>")
            elif line.startswith('### '):
                html.append(f"<h3>{line[4:]}</h3>")
            elif line.startswith('#### '):
                html.append(f"<h4>{line[5:]}</h4>")
            elif line.startswith('##### '):
                html.append(f"<h5>{line[6:]}</h5>")
            elif line.startswith('- '):
                html.append(f"<ul><li>{line[2:]}</li></ul>")
            elif line.startswith('```'):
                if html[-1] == '<pre>':
                    html.append('</pre>')
                else:
                    html.append('<pre>')
            elif line.startswith('**Overall Status:**'):
                status = line.split(':')[1].strip().lower()
                html.append(f"<p><strong>Overall Status:</strong> <span class='{status}'>{status.upper()}</span></p>")
            else:
                # Replace markdown bold with HTML bold
                line = line.replace('**', '<strong>', 1)
                if '**' in line:
                    line = line.replace('**', '</strong>', 1)
                html.append(f"<p>{line}</p>")
        
        html.append("</body>")
        html.append("</html>")
        
        return '\n'.join(html)


# Deprecated functions - kept for backward compatibility but will log warnings

def check_schema_validation() -> Dict[str, Any]:
    """
    DEPRECATED: Use DatabaseHealthCheck class instead.
    
    This function is maintained for backward compatibility.
    It will create a DatabaseHealthCheck instance and run schema validation.
    
    Returns:
        Dict with validation results
    """
    from database.utils.connection import supabase_connection
    
    logger.warning("Deprecated function check_schema_validation() called. Use DatabaseHealthCheck class instead.")
    
    with supabase_connection() as conn:
        health_check = DatabaseHealthCheck(conn)
        results = health_check.run_health_check(categories=['schema'])
        return results['results_by_category']['schema']


def check_data_integrity() -> Dict[str, Any]:
    """
    DEPRECATED: Use DatabaseHealthCheck class instead.
    
    This function is maintained for backward compatibility.
    It will create a DatabaseHealthCheck instance and run integrity checks.
    
    Returns:
        Dict with integrity check results
    """
    from database.utils.connection import supabase_connection
    
    logger.warning("Deprecated function check_data_integrity() called. Use DatabaseHealthCheck class instead.")
    
    with supabase_connection() as conn:
        health_check = DatabaseHealthCheck(conn)
        results = health_check.run_health_check(categories=['integrity'])
        return results['results_by_category']['integrity']


def check_performance() -> Dict[str, Any]:
    """
    DEPRECATED: Use DatabaseHealthCheck class instead.
    
    This function is maintained for backward compatibility.
    It will create a DatabaseHealthCheck instance and run performance checks.
    
    Returns:
        Dict with performance check results
    """
    from database.utils.connection import supabase_connection
    
    logger.warning("Deprecated function check_performance() called. Use DatabaseHealthCheck class instead.")
    
    with supabase_connection() as conn:
        health_check = DatabaseHealthCheck(conn)
        results = health_check.run_health_check(categories=['performance'])
        return results['results_by_category']['performance']


def generate_health_report() -> Dict[str, Any]:
    """
    DEPRECATED: Use DatabaseHealthCheck class instead.
    
    This function is maintained for backward compatibility.
    It will create a DatabaseHealthCheck instance and generate a report.
    
    Returns:
        Dict with complete health report
    """
    from database.utils.connection import supabase_connection
    
    logger.warning("Deprecated function generate_health_report() called. Use DatabaseHealthCheck class instead.")
    
    with supabase_connection() as conn:
        health_check = DatabaseHealthCheck(conn)
        return health_check.run_health_check()


def run_scheduled_health_check(interval_hours: int = 24) -> None:
    """
    DEPRECATED: Use DatabaseHealthCheck class instead.
    
    This function is maintained for backward compatibility.
    It will create a DatabaseHealthCheck instance and run scheduled checks.
    
    Args:
        interval_hours: Hours between checks
    """
    import threading
    from database.utils.connection import supabase_connection
    
    logger.warning("Deprecated function run_scheduled_health_check() called. Use DatabaseHealthCheck class instead.")
    
    def _scheduled_run():
        while True:
            try:
                logger.info("Running scheduled database health check")
                with supabase_connection() as conn:
                    health_check = DatabaseHealthCheck(conn)
                    report = health_check.run_health_check()
                    
                    if report['overall_status'] != 'passed':
                        logger.warning(f"Database health check issues found. Status: {report['overall_status']}")
                        
                        # If there are critical issues, could trigger alerts here
                        if report['overall_status'] == 'failed':
                            # Example: send_alert('Database health check failed', report['summary'])
                            pass
                    
                    logger.info(f"Database health check completed. Status: {report['overall_status']}")
            except Exception as e:
                logger.error(f"Error in scheduled health check: {str(e)}", exc_info=True)
            
            # Sleep for the specified interval
            time.sleep(interval_hours * 3600)
    
    # Start the scheduled health check in a separate thread
    thread = threading.Thread(target=_scheduled_run, daemon=True)
    thread.start()
    logger.info(f"Scheduled health check started (interval: {interval_hours} hours)")


def fix_common_issues(fix_type: str = 'all') -> Dict[str, Any]:
    """
    DEPRECATED: Use DatabaseHealthCheck class instead.
    
    This function is maintained for backward compatibility.
    In the future, this functionality should be moved to the validator modules.
    
    Args:
        fix_type: Type of fixes to apply ('orphaned', 'indexes', or 'all')
        
    Returns:
        Dict with fix results
    """
    from database.utils.connection import get_supabase_service_client
    
    logger.warning("Deprecated function fix_common_issues() called. Use DatabaseHealthCheck class instead.")
    
    results = {
        'status': 'success',
        'fixes_applied': 0,
        'details': {}
    }
    
    try:
        # Get the service client for admin operations
        conn = get_supabase_service_client()
        
        # Fix orphaned records
        if fix_type in ('orphaned', 'all'):
            # First generate a health report to identify issues
            integrity_report = check_data_integrity()
            
            if integrity_report['status'] != 'passed':
                for relationship_key, issue in integrity_report['details'].items():
                    if 'orphaned_records' in issue and issue['orphaned_records']:
                        # Extract table name from relationship key (e.g., "parent_to_dependent")
                        parts = relationship_key.split('_to_')
                        if len(parts) == 2:
                            dependent_table = parts[1]
                            
                            # Extract orphaned record IDs
                            orphaned_ids = [record['id'] for record in issue['orphaned_records']]
                            
                            if orphaned_ids:
                                # Delete orphaned records
                                delete_query = f'''
                                    DELETE FROM {dependent_table}
                                    WHERE id IN ({','.join([f"'{id}'" for id in orphaned_ids])});
                                '''
                                
                                conn.rpc('exec_sql', {'query': delete_query}).execute()
                                
                                results['fixes_applied'] += len(orphaned_ids)
                                results['details'][f"fixed_orphaned_{dependent_table}"] = {
                                    'count': len(orphaned_ids),
                                    'ids': orphaned_ids
                                }
        
        # Fix missing indexes
        if fix_type in ('indexes', 'all'):
            performance_report = check_performance()
            
            if 'missing_indexes' in performance_report['details']:
                missing_indexes = performance_report['details']['missing_indexes']
                
                for index_info in missing_indexes:
                    table_name = index_info['table_name']
                    column_name = index_info['column_name']
                    
                    # Create index
                    index_name = f"idx_{table_name}_{column_name}"
                    create_index_query = f'''
                        CREATE INDEX IF NOT EXISTS {index_name}
                        ON {table_name} ({column_name});
                    '''
                    
                    conn.rpc('exec_sql', {'query': create_index_query}).execute()
                    
                    results['fixes_applied'] += 1
                    results['details'][f"created_index_{table_name}_{column_name}"] = {
                        'index_name': index_name
                    }
        
        if results['fixes_applied'] == 0:
            results['status'] = 'no_fixes_needed'
            results['details']['message'] = 'No issues found that require fixing'
    except Exception as e:
        results['status'] = 'failed'
        results['details']['error'] = str(e)
        results['details']['traceback'] = traceback.format_exc()
    
    return results


def get_latest_health_report() -> Dict[str, Any]:
    """
    DEPRECATED: Use DatabaseHealthCheck class instead.
    
    This function is maintained for backward compatibility.
    It will return the most recent health report or generate a new one.
    
    Returns:
        Dict with the most recent health report
    """
    logger.warning("Deprecated function get_latest_health_report() called. Use DatabaseHealthCheck class instead.")
    
    try:
        report_files = [f for f in os.listdir(HEALTH_CHECK_REPORT_DIR) if f.startswith('health_report_') and f.endswith('.json')]
        
        if not report_files:
            # No reports found, generate a new one
            return generate_health_report()
        
        # Sort by timestamp (newest first)
        report_files.sort(reverse=True)
        latest_report_path = os.path.join(HEALTH_CHECK_REPORT_DIR, report_files[0])
        
        with open(latest_report_path, 'r') as f:
            report = json.load(f)
        
        return report
    except Exception as e:
        logger.error(f"Error getting latest health report: {str(e)}", exc_info=True)
        # If error, generate a new report
        return generate_health_report()


# Main function to run when the module is executed directly
if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    
    from database.utils.connection import supabase_connection
    
    with supabase_connection() as conn:
        health_check = DatabaseHealthCheck(conn)
        report = health_check.run_health_check()
        print(json.dumps(report, indent=2))