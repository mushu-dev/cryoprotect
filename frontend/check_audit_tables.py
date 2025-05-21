#!/usr/bin/env python3
"""
Script to check for existing audit and logging tables in the database.
"""

import sys
import json
import logging
from datetime import datetime
import db_utils
from psycopg2.extras import RealDictCursor

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

def check_tables():
    """Check for existing audit and logging tables."""
    try:
        # Get all tables
        query = """
        SELECT 
            table_schema,
            table_name
        FROM 
            information_schema.tables
        WHERE 
            table_schema IN ('public', 'audit', 'logging')
            AND table_type = 'BASE TABLE'
        ORDER BY
            table_schema, table_name;
        """
        
        tables = db_utils.execute_query(query, cursor_factory=RealDictCursor)
        
        # Check for audit schema
        audit_schema_query = """
        SELECT EXISTS (
            SELECT 1 
            FROM information_schema.schemata 
            WHERE schema_name = 'audit'
        );
        """
        audit_schema_exists = db_utils.execute_query(audit_schema_query)[0][0]
        
        # Check for logging schema
        logging_schema_query = """
        SELECT EXISTS (
            SELECT 1 
            FROM information_schema.schemata 
            WHERE schema_name = 'logging'
        );
        """
        logging_schema_exists = db_utils.execute_query(logging_schema_query)[0][0]
        
        # Get all tables that might be related to audit or logging
        related_tables_query = """
        SELECT 
            table_schema,
            table_name
        FROM 
            information_schema.tables
        WHERE 
            table_schema = 'public'
            AND table_type = 'BASE TABLE'
            AND (
                table_name LIKE '%audit%'
                OR table_name LIKE '%log%'
                OR table_name LIKE '%history%'
                OR table_name LIKE '%change%'
                OR table_name LIKE '%track%'
                OR table_name LIKE '%event%'
            )
        ORDER BY
            table_schema, table_name;
        """
        
        related_tables = db_utils.execute_query(related_tables_query, cursor_factory=RealDictCursor)
        
        # Check for existing audit triggers
        audit_triggers_query = """
        SELECT 
            trigger_schema,
            trigger_name,
            event_manipulation,
            event_object_schema,
            event_object_table,
            action_statement
        FROM 
            information_schema.triggers
        WHERE 
            trigger_name LIKE '%audit%'
            OR trigger_name LIKE '%log%'
            OR action_statement LIKE '%audit%'
            OR action_statement LIKE '%log%'
        ORDER BY
            event_object_schema,
            event_object_table,
            trigger_name;
        """
        
        audit_triggers = db_utils.execute_query(audit_triggers_query, cursor_factory=RealDictCursor)
        
        return {
            'tables': tables,
            'audit_schema_exists': audit_schema_exists,
            'logging_schema_exists': logging_schema_exists,
            'related_tables': related_tables,
            'audit_triggers': audit_triggers
        }
    except Exception as e:
        logger.error(f"Error checking tables: {e}")
        return {
            'error': str(e)
        }

def save_report(results, filename_prefix="audit_tables_report"):
    """Save the report to a file."""
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    filename = f"{filename_prefix}_{timestamp}.json"
    
    report = {
        'timestamp': datetime.now().isoformat(),
        'results': results,
        'summary': {
            'total_tables': len(results.get('tables', [])),
            'audit_schema_exists': results.get('audit_schema_exists', False),
            'logging_schema_exists': results.get('logging_schema_exists', False),
            'related_tables_count': len(results.get('related_tables', [])),
            'audit_triggers_count': len(results.get('audit_triggers', []))
        }
    }
    
    try:
        with open(filename, 'w') as f:
            json.dump(report, f, indent=2)
        logger.info(f"Report saved to {filename}")
        return filename
    except Exception as e:
        logger.error(f"Error saving report: {e}")
        return None

def main():
    """Main function."""
    # Check database connection
    if not db_utils.test_connection():
        logger.error("Database connection failed. Exiting.")
        sys.exit(1)
    
    # Check tables
    logger.info("Checking for audit and logging tables...")
    results = check_tables()
    
    if 'error' in results:
        logger.error(f"Error checking tables: {results['error']}")
        sys.exit(1)
    
    # Print summary
    logger.info(f"Total tables: {len(results['tables'])}")
    logger.info(f"Audit schema exists: {results['audit_schema_exists']}")
    logger.info(f"Logging schema exists: {results['logging_schema_exists']}")
    logger.info(f"Related tables: {len(results['related_tables'])}")
    logger.info(f"Audit triggers: {len(results['audit_triggers'])}")
    
    # If any related tables or audit triggers exist, print them
    if results['related_tables']:
        logger.info("\nRelated tables:")
        for table in results['related_tables']:
            logger.info(f"  {table['table_schema']}.{table['table_name']}")
    
    if results['audit_triggers']:
        logger.info("\nAudit triggers:")
        for trigger in results['audit_triggers']:
            logger.info(f"  {trigger['trigger_name']} on {trigger['event_object_schema']}.{trigger['event_object_table']} ({trigger['event_manipulation']})")
    
    # Save report
    report_file = save_report(results)
    if report_file:
        logger.info(f"Report saved to {report_file}")
    
    sys.exit(0)

if __name__ == "__main__":
    main()