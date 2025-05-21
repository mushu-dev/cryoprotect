#!/usr/bin/env python3
"""
Script to implement and test the audit logging system.
This script applies the SQL migration and tests the functionality.
"""

import os
import sys
import json
import logging
from datetime import datetime
import argparse
from psycopg2.extras import RealDictCursor
import db_utils

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

def execute_migration():
    """Execute the migration to implement audit logging."""
    sql_file = os.path.join(
        os.path.dirname(os.path.realpath(__file__)), 
        'migrations', 
        '031_implement_audit_logging.sql'
    )
    
    try:
        if db_utils.execute_sql_file(sql_file):
            logger.info("Successfully executed migration to implement audit logging")
            return True
        else:
            logger.error("Failed to execute migration")
            return False
    except Exception as e:
        logger.error(f"Error executing migration: {e}")
        return False

def check_audit_schema():
    """Check if the audit schema and tables were created successfully."""
    try:
        # Check if audit schema exists
        schema_query = """
        SELECT EXISTS (
            SELECT 1 
            FROM information_schema.schemata 
            WHERE schema_name = 'audit'
        );
        """
        schema_exists = db_utils.execute_query(schema_query)[0][0]
        
        # Check if audit_log table exists
        audit_log_query = """
        SELECT EXISTS (
            SELECT 1 
            FROM information_schema.tables 
            WHERE table_schema = 'audit' 
            AND table_name = 'audit_log'
        );
        """
        audit_log_exists = db_utils.execute_query(audit_log_query)[0][0]
        
        # Check if application_log table exists
        app_log_query = """
        SELECT EXISTS (
            SELECT 1 
            FROM information_schema.tables 
            WHERE table_schema = 'audit' 
            AND table_name = 'application_log'
        );
        """
        app_log_exists = db_utils.execute_query(app_log_query)[0][0]
        
        # Get audit_log columns if it exists
        audit_log_columns = []
        if audit_log_exists:
            columns_query = """
            SELECT column_name, data_type
            FROM information_schema.columns
            WHERE table_schema = 'audit'
            AND table_name = 'audit_log';
            """
            audit_log_columns = db_utils.execute_query(columns_query, cursor_factory=RealDictCursor)
        
        # Get application_log columns if it exists
        app_log_columns = []
        if app_log_exists:
            columns_query = """
            SELECT column_name, data_type
            FROM information_schema.columns
            WHERE table_schema = 'audit'
            AND table_name = 'application_log';
            """
            app_log_columns = db_utils.execute_query(columns_query, cursor_factory=RealDictCursor)
        
        return {
            'audit_schema_exists': schema_exists,
            'audit_log_exists': audit_log_exists,
            'application_log_exists': app_log_exists,
            'audit_log_columns': audit_log_columns,
            'application_log_columns': app_log_columns
        }
    except Exception as e:
        logger.error(f"Error checking audit schema: {e}")
        return {
            'error': str(e)
        }

def check_audit_functions():
    """Check if the audit functions were created successfully."""
    try:
        # List of functions to check
        functions = [
            'audit.audit_trigger_func',
            'audit.enable_table_auditing',
            'audit.enable_schema_auditing',
            'audit.get_record_history',
            'audit.log_application_event',
            'audit.cleanup_old_logs',
            'audit.get_request_context',
            'audit.changed_fields'
        ]
        
        function_results = {}
        
        for function_name in functions:
            schema, name = function_name.split('.')
            query = """
            SELECT EXISTS (
                SELECT 1
                FROM pg_proc p
                JOIN pg_namespace n ON p.pronamespace = n.oid
                WHERE n.nspname = %s
                AND p.proname = %s
            );
            """
            exists = db_utils.execute_query(query, (schema, name))[0][0]
            function_results[function_name] = exists
        
        return function_results
    except Exception as e:
        logger.error(f"Error checking audit functions: {e}")
        return {
            'error': str(e)
        }

def check_audit_triggers():
    """Check if audit triggers were created for the target tables."""
    try:
        # Get all audit triggers
        query = """
        SELECT 
            trigger_schema,
            trigger_name,
            event_object_schema,
            event_object_table,
            action_statement
        FROM 
            information_schema.triggers
        WHERE 
            trigger_name LIKE 'audit_trigger_%'
        ORDER BY
            event_object_schema,
            event_object_table;
        """
        
        triggers = db_utils.execute_query(query, cursor_factory=RealDictCursor)
        
        # Check if triggers exist for core tables
        core_tables = [
            'molecules',
            'mixtures',
            'molecular_properties',
            'teams',
            'team_members',
            'mixture_components'
        ]
        
        core_table_triggers = {}
        
        for table in core_tables:
            has_trigger = False
            for trigger in triggers:
                if trigger['event_object_table'] == table:
                    has_trigger = True
                    break
            core_table_triggers[table] = has_trigger
        
        return {
            'triggers': triggers,
            'core_table_triggers': core_table_triggers
        }
    except Exception as e:
        logger.error(f"Error checking audit triggers: {e}")
        return {
            'error': str(e)
        }

def check_audit_views():
    """Check if the audit views were created successfully."""
    try:
        # List of views to check
        views = [
            'audit.recent_activity',
            'audit.activity_summary',
            'audit.user_activity',
            'audit.recent_app_logs',
            'audit.app_log_summary'
        ]
        
        view_results = {}
        
        for view_name in views:
            schema, name = view_name.split('.')
            query = """
            SELECT EXISTS (
                SELECT 1
                FROM information_schema.views
                WHERE table_schema = %s
                AND table_name = %s
            );
            """
            exists = db_utils.execute_query(query, (schema, name))[0][0]
            view_results[view_name] = exists
        
        return view_results
    except Exception as e:
        logger.error(f"Error checking audit views: {e}")
        return {
            'error': str(e)
        }

def test_audit_logging():
    """Test the audit logging system with a sample operation."""
    try:
        # Create a test record
        test_id = datetime.now().strftime('%Y%m%d%H%M%S')
        create_query = f"""
        INSERT INTO public.teams (name, description, created_by)
        VALUES ('Audit Test Team {test_id}', 'Created to test audit logging', auth.uid())
        RETURNING id;
        """
        
        # Execute the insert
        create_result = db_utils.execute_query(create_query)
        
        if not create_result:
            logger.error("Failed to create test record")
            return {
                'success': False,
                'error': 'Failed to create test record'
            }
        
        team_id = create_result[0][0]
        logger.info(f"Created test team with ID: {team_id}")
        
        # Update the test record
        update_query = f"""
        UPDATE public.teams
        SET description = 'Updated to test audit logging'
        WHERE id = '{team_id}'
        RETURNING id;
        """
        
        # Execute the update
        update_result = db_utils.execute_query(update_query)
        
        if not update_result:
            logger.error("Failed to update test record")
            return {
                'success': False,
                'error': 'Failed to update test record'
            }
        
        logger.info(f"Updated test team with ID: {team_id}")
        
        # Check if audit logs were created
        audit_query = f"""
        SELECT 
            id,
            table_name,
            action,
            row_id,
            changed_fields,
            old_data,
            new_data,
            recorded_at
        FROM 
            audit.audit_log
        WHERE 
            row_id = '{team_id}'
        ORDER BY
            recorded_at;
        """
        
        audit_logs = db_utils.execute_query(audit_query, cursor_factory=RealDictCursor)
        
        if not audit_logs:
            logger.error("No audit logs found for the test operation")
            return {
                'success': False,
                'error': 'No audit logs found for the test operation'
            }
        
        logger.info(f"Found {len(audit_logs)} audit logs for the test operation")
        
        # Test application logging
        app_log_query = """
        SELECT audit.log_application_event(
            'TEST',
            'INFO',
            'Test application logging',
            'AUDIT_TEST',
            '{"test": true}'::jsonb
        );
        """
        
        app_log_result = db_utils.execute_query(app_log_query)
        
        if not app_log_result:
            logger.error("Failed to create application log")
            return {
                'success': False,
                'error': 'Failed to create application log'
            }
        
        app_log_id = app_log_result[0][0]
        logger.info(f"Created application log with ID: {app_log_id}")
        
        # Check if the application log was created
        check_app_log_query = f"""
        SELECT 
            id,
            event_type,
            event_severity,
            message,
            details,
            recorded_at
        FROM 
            audit.application_log
        WHERE 
            id = {app_log_id};
        """
        
        app_logs = db_utils.execute_query(check_app_log_query, cursor_factory=RealDictCursor)
        
        if not app_logs:
            logger.error("No application logs found for the test operation")
            return {
                'success': False,
                'error': 'No application logs found for the test operation'
            }
        
        logger.info(f"Found application log with ID: {app_log_id}")
        
        # Clean up the test record
        delete_query = f"""
        DELETE FROM public.teams
        WHERE id = '{team_id}'
        RETURNING id;
        """
        
        delete_result = db_utils.execute_query(delete_query)
        
        if not delete_result:
            logger.error("Failed to delete test record")
            return {
                'success': False,
                'error': 'Failed to delete test record'
            }
        
        logger.info(f"Deleted test team with ID: {team_id}")
        
        # Check if delete action was audited
        final_audit_query = f"""
        SELECT 
            id,
            action,
            recorded_at
        FROM 
            audit.audit_log
        WHERE 
            row_id = '{team_id}'
        ORDER BY
            recorded_at DESC
        LIMIT 1;
        """
        
        final_audit = db_utils.execute_query(final_audit_query, cursor_factory=RealDictCursor)
        
        delete_audited = False
        if final_audit and final_audit[0]['action'] == 'DELETE':
            delete_audited = True
            logger.info("Delete action was audited correctly")
        else:
            logger.error("Delete action was not audited")
        
        return {
            'success': True,
            'team_id': team_id,
            'audit_logs': audit_logs,
            'app_log_id': app_log_id,
            'app_logs': app_logs,
            'delete_audited': delete_audited
        }
    except Exception as e:
        logger.error(f"Error testing audit logging: {e}")
        return {
            'success': False,
            'error': str(e)
        }

def save_report(schema_check, function_check, trigger_check, view_check, test_results, 
                filename_prefix="audit_logging_implementation_report"):
    """Save the report to a file."""
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    filename = f"{filename_prefix}_{timestamp}.json"
    
    # Calculate success metrics
    schema_success = schema_check.get('audit_schema_exists', False) and \
                     schema_check.get('audit_log_exists', False) and \
                     schema_check.get('application_log_exists', False)
    
    function_success = all(function_check.values()) if function_check and 'error' not in function_check else False
    
    trigger_success = all(trigger_check.get('core_table_triggers', {}).values()) if trigger_check and 'error' not in trigger_check else False
    
    view_success = all(view_check.values()) if view_check and 'error' not in view_check else False
    
    test_success = test_results.get('success', False) if test_results else False
    
    report = {
        'timestamp': datetime.now().isoformat(),
        'schema_check': schema_check,
        'function_check': function_check,
        'trigger_check': trigger_check,
        'view_check': view_check,
        'test_results': test_results,
        'summary': {
            'schema_success': schema_success,
            'function_success': function_success,
            'trigger_success': trigger_success,
            'view_success': view_success,
            'test_success': test_success,
            'overall_success': schema_success and function_success and trigger_success and view_success and test_success
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
    """Main function to run the script."""
    parser = argparse.ArgumentParser(description="Implement and test the audit logging system.")
    parser.add_argument("--test-only", action="store_true", help="Only test the existing audit logging system without applying the migration")
    args = parser.parse_args()
    
    # Check database connection
    if not db_utils.test_connection():
        logger.error("Database connection failed. Exiting.")
        sys.exit(1)
    
    if not args.test_only:
        # Execute the migration
        logger.info("Executing migration to implement audit logging...")
        success = execute_migration()
        
        if not success:
            logger.error("Migration failed. Exiting.")
            sys.exit(1)
    
    # Check audit schema and tables
    logger.info("Checking audit schema and tables...")
    schema_check = check_audit_schema()
    
    if 'error' in schema_check:
        logger.error(f"Error checking audit schema: {schema_check['error']}")
        sys.exit(1)
    
    # Check audit functions
    logger.info("Checking audit functions...")
    function_check = check_audit_functions()
    
    if 'error' in function_check:
        logger.error(f"Error checking audit functions: {function_check['error']}")
        sys.exit(1)
    
    # Check audit triggers
    logger.info("Checking audit triggers...")
    trigger_check = check_audit_triggers()
    
    if 'error' in trigger_check:
        logger.error(f"Error checking audit triggers: {trigger_check['error']}")
        sys.exit(1)
    
    # Check audit views
    logger.info("Checking audit views...")
    view_check = check_audit_views()
    
    if 'error' in view_check:
        logger.error(f"Error checking audit views: {view_check['error']}")
        sys.exit(1)
    
    # Test audit logging
    logger.info("Testing audit logging...")
    test_results = test_audit_logging()
    
    # Save report
    report_file = save_report(
        schema_check, 
        function_check, 
        trigger_check, 
        view_check, 
        test_results
    )
    
    if report_file:
        logger.info(f"Report saved to {report_file}")
    
    # Check overall success
    schema_success = schema_check.get('audit_schema_exists', False) and \
                     schema_check.get('audit_log_exists', False) and \
                     schema_check.get('application_log_exists', False)
    
    function_success = all(function_check.values()) if function_check and 'error' not in function_check else False
    
    trigger_success = all(trigger_check.get('core_table_triggers', {}).values()) if trigger_check and 'error' not in trigger_check else False
    
    view_success = all(view_check.values()) if view_check and 'error' not in view_check else False
    
    test_success = test_results.get('success', False) if test_results else False
    
    overall_success = schema_success and function_success and trigger_success and view_success and test_success
    
    if overall_success:
        logger.info("Audit logging implementation successful!")
        print(f"\n✅ Audit logging implementation successful!")
        print(f"   - Audit schema and tables: {'✅' if schema_success else '❌'}")
        print(f"   - Audit functions: {'✅' if function_success else '❌'}")
        print(f"   - Audit triggers: {'✅' if trigger_success else '❌'}")
        print(f"   - Audit views: {'✅' if view_success else '❌'}")
        print(f"   - Audit logging test: {'✅' if test_success else '❌'}")
        sys.exit(0)
    else:
        logger.warning("Audit logging implementation has issues. Check the report for details.")
        print(f"\n⚠️ Audit logging implementation has issues:")
        print(f"   - Audit schema and tables: {'✅' if schema_success else '❌'}")
        print(f"   - Audit functions: {'✅' if function_success else '❌'}")
        print(f"   - Audit triggers: {'✅' if trigger_success else '❌'}")
        print(f"   - Audit views: {'✅' if view_success else '❌'}")
        print(f"   - Audit logging test: {'✅' if test_success else '❌'}")
        sys.exit(1)

if __name__ == "__main__":
    main()