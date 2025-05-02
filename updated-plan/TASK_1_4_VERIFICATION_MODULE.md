# Task 1.4: Build Verification Module

## Objective
Create a comprehensive database verification module that provides tools for validating database integrity, schema consistency, and data quality.

## Context
The project currently has multiple scripts for verifying database state, but they are scattered and inconsistent. By creating a unified verification module, we'll improve the reliability of database operations, enable automated testing, and provide better tools for diagnosing issues.

## Acceptance Criteria
- A modular verification package with well-defined components
- Support for schema verification, constraint validation, and data quality checks
- CLI interface for running verifications
- Detailed reporting of verification results
- Support for different verification levels (basic, standard, comprehensive)
- Integration with the main database package
- Clear documentation and examples

## Implementation Steps

1. Create the verification module structure:
   ```
   database/verification/
   ├── __init__.py          # Package exports
   ├── runner.py            # Main entry point
   ├── schema.py            # Schema verification
   ├── constraints.py       # Constraint verification
   ├── data.py              # Data quality verification
   ├── reporters/           # Report generation
   │   ├── __init__.py
   │   ├── console.py       # Console reporter
   │   ├── json.py          # JSON reporter
   │   └── html.py          # HTML reporter
   └── utils.py             # Verification utilities
   ```

2. Implement the verification module `__init__.py`:
   ```python
   """
   Database verification module for CryoProtect v2.
   
   This module provides tools for verifying database integrity,
   schema consistency, and data quality.
   """
   
   from database.verification.runner import (
       verify_database,
       verify_schema,
       verify_constraints,
       verify_data_quality
   )
   
   __all__ = [
       'verify_database',
       'verify_schema',
       'verify_constraints',
       'verify_data_quality'
   ]
   ```

3. Implement the main verification runner module:
   ```python
   """
   Main entry point for database verification operations.
   
   This module provides high-level functions for verifying database
   integrity, schema consistency, and data quality.
   """
   
   import logging
   import os
   import json
   from datetime import datetime
   from typing import Dict, List, Optional, Set, Union, Any
   
   from database.utils.connection import create_connection
   from database.verification.schema import run_schema_verification
   from database.verification.constraints import run_constraint_verification
   from database.verification.data import run_data_quality_verification
   
   logger = logging.getLogger(__name__)
   
   # Verification levels
   BASIC = 'basic'        # Fast, essential checks only
   STANDARD = 'standard'  # Default level, balanced
   COMPREHENSIVE = 'comprehensive'  # Thorough checks
   
   def verify_database(
       conn: Optional[Any] = None,
       config: Optional[Dict] = None,
       level: str = STANDARD,
       include_modules: Optional[List[str]] = None,
       exclude_modules: Optional[List[str]] = None,
       report_format: str = 'json',
       output_file: Optional[str] = None
   ) -> Dict:
       """
       Verify database integrity, schema consistency, and data quality.
       
       Args:
           conn: Database connection
           config: Optional configuration dictionary
           level: Verification level (basic, standard, comprehensive)
           include_modules: Optional list of specific modules to include
           exclude_modules: Optional list of specific modules to exclude
           report_format: Report format (console, json, html)
           output_file: Optional path to save report
           
       Returns:
           Verification results dictionary
       """
       logger.info(f"Starting database verification at {level} level")
       
       if conn is None:
           conn = create_connection(config)
           
       # Determine which modules to run
       all_modules = {'schema', 'constraints', 'data_quality'}
       
       if include_modules:
           modules_to_run = set(include_modules) & all_modules
       elif exclude_modules:
           modules_to_run = all_modules - set(exclude_modules)
       else:
           modules_to_run = all_modules
           
       # Run verifications
       results = {
           'timestamp': datetime.now().isoformat(),
           'level': level,
           'modules': list(modules_to_run),
           'results': {}
       }
       
       if 'schema' in modules_to_run:
           results['results']['schema'] = verify_schema(conn, level)
           
       if 'constraints' in modules_to_run:
           results['results']['constraints'] = verify_constraints(conn, level)
           
       if 'data_quality' in modules_to_run:
           results['results']['data_quality'] = verify_data_quality(conn, level)
           
       # Calculate summary
       success = all(
           module_result.get('success', False)
           for module_result in results['results'].values()
       )
       
       results['success'] = success
       
       # Generate report
       if output_file:
           if report_format == 'json':
               with open(output_file, 'w') as f:
                   json.dump(results, f, indent=2)
           elif report_format == 'html':
               from database.verification.reporters.html import generate_html_report
               html_report = generate_html_report(results)
               with open(output_file, 'w') as f:
                   f.write(html_report)
                   
       # Log summary
       logger_method = logger.info if success else logger.error
       logger_method(f"Verification completed with success={success}")
       
       return results
   
   def verify_schema(
       conn: Optional[Any] = None,
       level: str = STANDARD,
       config: Optional[Dict] = None
   ) -> Dict:
       """
       Verify database schema.
       
       Args:
           conn: Database connection
           level: Verification level
           config: Optional configuration dictionary
           
       Returns:
           Schema verification results
       """
       logger.info(f"Verifying database schema at {level} level")
       
       if conn is None:
           conn = create_connection(config)
           
       return run_schema_verification(conn, level)
   
   def verify_constraints(
       conn: Optional[Any] = None,
       level: str = STANDARD,
       config: Optional[Dict] = None
   ) -> Dict:
       """
       Verify database constraints.
       
       Args:
           conn: Database connection
           level: Verification level
           config: Optional configuration dictionary
           
       Returns:
           Constraint verification results
       """
       logger.info(f"Verifying database constraints at {level} level")
       
       if conn is None:
           conn = create_connection(config)
           
       return run_constraint_verification(conn, level)
   
   def verify_data_quality(
       conn: Optional[Any] = None,
       level: str = STANDARD,
       config: Optional[Dict] = None
   ) -> Dict:
       """
       Verify data quality.
       
       Args:
           conn: Database connection
           level: Verification level
           config: Optional configuration dictionary
           
       Returns:
           Data quality verification results
       """
       logger.info(f"Verifying data quality at {level} level")
       
       if conn is None:
           conn = create_connection(config)
           
       return run_data_quality_verification(conn, level)
   
   def main():
       """
       CLI entry point for database verification.
       """
       import argparse
       
       parser = argparse.ArgumentParser(description='Verify CryoProtect database')
       parser.add_argument(
           '--level',
           choices=['basic', 'standard', 'comprehensive'],
           default='standard',
           help='Verification level'
       )
       parser.add_argument(
           '--include',
           nargs='+',
           choices=['schema', 'constraints', 'data_quality'],
           help='Specific verification modules to include'
       )
       parser.add_argument(
           '--exclude',
           nargs='+',
           choices=['schema', 'constraints', 'data_quality'],
           help='Specific verification modules to exclude'
       )
       parser.add_argument(
           '--format',
           choices=['console', 'json', 'html'],
           default='console',
           help='Report format'
       )
       parser.add_argument(
           '--output',
           help='Output file path for report'
       )
       
       args = parser.parse_args()
       
       results = verify_database(
           level=args.level,
           include_modules=args.include,
           exclude_modules=args.exclude,
           report_format=args.format,
           output_file=args.output
       )
       
       if args.format == 'console':
           print(f"Verification Results (Level: {args.level})")
           print("-" * 60)
           for module, module_results in results['results'].items():
               status = "PASSED" if module_results.get('success', False) else "FAILED"
               print(f"{module}: {status}")
               
               if 'issues' in module_results and module_results['issues']:
                   print("  Issues:")
                   for issue in module_results['issues']:
                       print(f"  - {issue['message']}")
               print()
               
           overall = "PASSED" if results['success'] else "FAILED"
           print(f"Overall: {overall}")
       
   if __name__ == '__main__':
       main()
   ```

4. Implement the schema verification module:
   ```python
   """
   Schema verification module.
   
   This module provides functions for verifying database schema consistency.
   """
   
   import logging
   from typing import Any, Dict, List, Optional, Set, Tuple
   
   logger = logging.getLogger(__name__)
   
   def _get_expected_tables() -> Dict[str, List[Dict]]:
       """
       Get expected tables and their columns.
       
       Returns:
           Dictionary mapping table names to column definitions
       """
       # This would ideally load from a schema definition file
       return {
           'molecules': [
               {'name': 'id', 'type': 'uuid', 'nullable': False},
               {'name': 'name', 'type': 'varchar', 'nullable': False},
               {'name': 'formula', 'type': 'varchar', 'nullable': True},
               {'name': 'smiles', 'type': 'text', 'nullable': True},
               {'name': 'created_at', 'type': 'timestamp with time zone', 'nullable': True},
               {'name': 'updated_at', 'type': 'timestamp with time zone', 'nullable': True}
           ],
           'mixtures': [
               {'name': 'id', 'type': 'uuid', 'nullable': False},
               {'name': 'name', 'type': 'varchar', 'nullable': False},
               {'name': 'description', 'type': 'text', 'nullable': True},
               {'name': 'created_at', 'type': 'timestamp with time zone', 'nullable': True},
               {'name': 'updated_at', 'type': 'timestamp with time zone', 'nullable': True}
           ],
           'mixture_components': [
               {'name': 'id', 'type': 'uuid', 'nullable': False},
               {'name': 'mixture_id', 'type': 'uuid', 'nullable': False},
               {'name': 'molecule_id', 'type': 'uuid', 'nullable': False},
               {'name': 'concentration', 'type': 'numeric', 'nullable': True},
               {'name': 'units', 'type': 'varchar', 'nullable': True},
               {'name': 'created_at', 'type': 'timestamp with time zone', 'nullable': True},
               {'name': 'updated_at', 'type': 'timestamp with time zone', 'nullable': True}
           ]
       }
   
   def _get_actual_tables(conn: Any) -> Dict[str, List[Dict]]:
       """
       Get actual tables and their columns from the database.
       
       Args:
           conn: Database connection
           
       Returns:
           Dictionary mapping table names to column definitions
       """
       tables = {}
       
       # Get list of tables
       table_query = """
           SELECT tablename FROM pg_tables
           WHERE schemaname='public'
           ORDER BY tablename
       """
       table_result = conn.sql(table_query).execute()
       
       for table_row in table_result.data:
           table_name = table_row['tablename']
           
           # Get columns for this table
           column_query = """
               SELECT
                   column_name,
                   data_type,
                   is_nullable
               FROM information_schema.columns
               WHERE table_schema='public' AND table_name=$1
               ORDER BY ordinal_position
           """
           column_result = conn.sql(column_query, {"params": [table_name]}).execute()
           
           columns = []
           for column_row in column_result.data:
               columns.append({
                   'name': column_row['column_name'],
                   'type': column_row['data_type'],
                   'nullable': column_row['is_nullable'] == 'YES'
               })
               
           tables[table_name] = columns
           
       return tables
   
   def _verify_table_presence(
       expected_tables: Dict[str, List[Dict]],
       actual_tables: Dict[str, List[Dict]]
   ) -> List[Dict]:
       """
       Verify that all expected tables are present.
       
       Args:
           expected_tables: Expected tables and columns
           actual_tables: Actual tables and columns
           
       Returns:
           List of issues found
       """
       issues = []
       
       # Check for missing tables
       missing_tables = set(expected_tables.keys()) - set(actual_tables.keys())
       for table in missing_tables:
           issues.append({
               'type': 'missing_table',
               'severity': 'error',
               'message': f"Table '{table}' is missing"
           })
           
       return issues
   
   def _verify_table_columns(
       expected_tables: Dict[str, List[Dict]],
       actual_tables: Dict[str, List[Dict]]
   ) -> List[Dict]:
       """
       Verify that tables have expected columns.
       
       Args:
           expected_tables: Expected tables and columns
           actual_tables: Actual tables and columns
           
       Returns:
           List of issues found
       """
       issues = []
       
       for table_name, expected_columns in expected_tables.items():
           if table_name not in actual_tables:
               continue  # Skip missing tables (already reported)
               
           actual_columns = actual_tables[table_name]
           actual_column_names = {col['name'] for col in actual_columns}
           
           # Check for missing columns
           for expected_column in expected_columns:
               column_name = expected_column['name']
               
               if column_name not in actual_column_names:
                   issues.append({
                       'type': 'missing_column',
                       'severity': 'error',
                       'message': f"Column '{column_name}' is missing from table '{table_name}'"
                   })
               else:
                   # Check column properties
                   actual_column = next(col for col in actual_columns if col['name'] == column_name)
                   
                   # Check type
                   if expected_column['type'] not in actual_column['type']:
                       issues.append({
                           'type': 'column_type_mismatch',
                           'severity': 'error',
                           'message': (
                               f"Column '{column_name}' in table '{table_name}' "
                               f"has type '{actual_column['type']}', "
                               f"expected '{expected_column['type']}'"
                           )
                       })
                       
                   # Check nullability
                   if expected_column['nullable'] != actual_column['nullable']:
                       severity = 'warning' if expected_column['nullable'] else 'error'
                       issues.append({
                           'type': 'column_nullability_mismatch',
                           'severity': severity,
                           'message': (
                               f"Column '{column_name}' in table '{table_name}' "
                               f"has nullability={actual_column['nullable']}, "
                               f"expected {expected_column['nullable']}"
                           )
                       })
           
       return issues
   
   def run_schema_verification(conn: Any, level: str = 'standard') -> Dict:
       """
       Run schema verification checks.
       
       Args:
           conn: Database connection
           level: Verification level
           
       Returns:
           Verification results
       """
       logger.info(f"Running schema verification at {level} level")
       
       issues = []
       
       try:
           expected_tables = _get_expected_tables()
           actual_tables = _get_actual_tables(conn)
           
           # Verify table presence
           issues.extend(_verify_table_presence(expected_tables, actual_tables))
           
           # Verify table columns
           issues.extend(_verify_table_columns(expected_tables, actual_tables))
           
           # Additional checks for comprehensive level
           if level == 'comprehensive':
               # Additional checks here...
               pass
               
           success = not any(issue['severity'] == 'error' for issue in issues)
           
           return {
               'success': success,
               'issues': issues
           }
       except Exception as e:
           logger.error(f"Error during schema verification: {str(e)}")
           return {
               'success': False,
               'issues': [
                   {
                       'type': 'verification_error',
                       'severity': 'error',
                       'message': f"Error during schema verification: {str(e)}"
                   }
               ]
           }
   ```

5. Implement the constraint verification module:
   ```python
   """
   Constraint verification module.
   
   This module provides functions for verifying database constraint integrity.
   """
   
   import logging
   from typing import Any, Dict, List, Optional, Set, Tuple
   
   logger = logging.getLogger(__name__)
   
   def _get_expected_foreign_keys() -> List[Dict]:
       """
       Get expected foreign key constraints.
       
       Returns:
           List of expected foreign key definitions
       """
       # This would ideally load from a schema definition file
       return [
           {
               'table': 'mixture_components',
               'column': 'mixture_id',
               'references_table': 'mixtures',
               'references_column': 'id'
           },
           {
               'table': 'mixture_components',
               'column': 'molecule_id',
               'references_table': 'molecules',
               'references_column': 'id'
           }
       ]
   
   def _get_actual_foreign_keys(conn: Any) -> List[Dict]:
       """
       Get actual foreign key constraints from the database.
       
       Args:
           conn: Database connection
           
       Returns:
           List of actual foreign key definitions
       """
       query = """
           SELECT
               tc.table_name,
               kcu.column_name,
               ccu.table_name AS references_table,
               ccu.column_name AS references_column
           FROM information_schema.table_constraints AS tc
           JOIN information_schema.key_column_usage AS kcu
             ON tc.constraint_name = kcu.constraint_name
           JOIN information_schema.constraint_column_usage AS ccu
             ON ccu.constraint_name = tc.constraint_name
           WHERE tc.constraint_type = 'FOREIGN KEY'
           AND tc.table_schema = 'public'
       """
       result = conn.sql(query).execute()
       
       foreign_keys = []
       for row in result.data:
           foreign_keys.append({
               'table': row['table_name'],
               'column': row['column_name'],
               'references_table': row['references_table'],
               'references_column': row['references_column']
           })
           
       return foreign_keys
   
   def _verify_foreign_keys(
       expected_foreign_keys: List[Dict],
       actual_foreign_keys: List[Dict]
   ) -> List[Dict]:
       """
       Verify that expected foreign keys are present.
       
       Args:
           expected_foreign_keys: Expected foreign key constraints
           actual_foreign_keys: Actual foreign key constraints
           
       Returns:
           List of issues found
       """
       issues = []
       
       for expected_fk in expected_foreign_keys:
           found = False
           for actual_fk in actual_foreign_keys:
               if (
                   expected_fk['table'] == actual_fk['table']
                   and expected_fk['column'] == actual_fk['column']
                   and expected_fk['references_table'] == actual_fk['references_table']
                   and expected_fk['references_column'] == actual_fk['references_column']
               ):
                   found = True
                   break
                   
           if not found:
               issues.append({
                   'type': 'missing_foreign_key',
                   'severity': 'error',
                   'message': (
                       f"Foreign key constraint missing: "
                       f"{expected_fk['table']}({expected_fk['column']}) "
                       f"-> {expected_fk['references_table']}({expected_fk['references_column']})"
                   )
               })
               
       return issues
   
   def _verify_integrity_constraints(conn: Any) -> List[Dict]:
       """
       Verify data integrity constraints.
       
       Args:
           conn: Database connection
           
       Returns:
           List of issues found
       """
       issues = []
       
       # Check for orphaned records in mixture_components
       query = """
           SELECT COUNT(*) AS count
           FROM mixture_components
           WHERE mixture_id NOT IN (SELECT id FROM mixtures)
       """
       result = conn.sql(query).execute()
       orphaned_count = result.data[0]['count'] if result.data else 0
       
       if orphaned_count > 0:
           issues.append({
               'type': 'orphaned_records',
               'severity': 'error',
               'message': f"Found {orphaned_count} orphaned records in mixture_components (mixture_id)"
           })
           
       # Check for orphaned records in mixture_components (molecule)
       query = """
           SELECT COUNT(*) AS count
           FROM mixture_components
           WHERE molecule_id NOT IN (SELECT id FROM molecules)
       """
       result = conn.sql(query).execute()
       orphaned_count = result.data[0]['count'] if result.data else 0
       
       if orphaned_count > 0:
           issues.append({
               'type': 'orphaned_records',
               'severity': 'error',
               'message': f"Found {orphaned_count} orphaned records in mixture_components (molecule_id)"
           })
           
       return issues
   
   def run_constraint_verification(conn: Any, level: str = 'standard') -> Dict:
       """
       Run constraint verification checks.
       
       Args:
           conn: Database connection
           level: Verification level
           
       Returns:
           Verification results
       """
       logger.info(f"Running constraint verification at {level} level")
       
       issues = []
       
       try:
           # Verify foreign keys
           expected_foreign_keys = _get_expected_foreign_keys()
           actual_foreign_keys = _get_actual_foreign_keys(conn)
           issues.extend(_verify_foreign_keys(expected_foreign_keys, actual_foreign_keys))
           
           # Verify integrity constraints (standard and comprehensive levels)
           if level in ('standard', 'comprehensive'):
               issues.extend(_verify_integrity_constraints(conn))
           
           # Additional checks for comprehensive level
           if level == 'comprehensive':
               # Additional checks here...
               pass
               
           success = not any(issue['severity'] == 'error' for issue in issues)
           
           return {
               'success': success,
               'issues': issues
           }
       except Exception as e:
           logger.error(f"Error during constraint verification: {str(e)}")
           return {
               'success': False,
               'issues': [
                   {
                       'type': 'verification_error',
                       'severity': 'error',
                       'message': f"Error during constraint verification: {str(e)}"
                   }
               ]
           }
   ```

6. Implement the data quality verification module:
   ```python
   """
   Data quality verification module.
   
   This module provides functions for verifying data quality.
   """
   
   import logging
   from typing import Any, Dict, List, Optional, Set, Tuple
   
   logger = logging.getLogger(__name__)
   
   def _check_required_fields(conn: Any) -> List[Dict]:
       """
       Check for missing required fields.
       
       Args:
           conn: Database connection
           
       Returns:
           List of issues found
       """
       issues = []
       
       # Check for missing required fields in molecules
       query = """
           SELECT COUNT(*) AS count
           FROM molecules
           WHERE name IS NULL OR name = ''
       """
       result = conn.sql(query).execute()
       missing_count = result.data[0]['count'] if result.data else 0
       
       if missing_count > 0:
           issues.append({
               'type': 'missing_required_field',
               'severity': 'error',
               'message': f"Found {missing_count} molecule records with missing name"
           })
           
       # Check for missing required fields in mixtures
       query = """
           SELECT COUNT(*) AS count
           FROM mixtures
           WHERE name IS NULL OR name = ''
       """
       result = conn.sql(query).execute()
       missing_count = result.data[0]['count'] if result.data else 0
       
       if missing_count > 0:
           issues.append({
               'type': 'missing_required_field',
               'severity': 'error',
               'message': f"Found {missing_count} mixture records with missing name"
           })
           
       return issues
   
   def _check_data_consistency(conn: Any) -> List[Dict]:
       """
       Check for data consistency issues.
       
       Args:
           conn: Database connection
           
       Returns:
           List of issues found
       """
       issues = []
       
       # Check for mixtures with no components
       query = """
           SELECT COUNT(*) AS count
           FROM mixtures m
           WHERE NOT EXISTS (
               SELECT 1 FROM mixture_components mc
               WHERE mc.mixture_id = m.id
           )
       """
       result = conn.sql(query).execute()
       empty_count = result.data[0]['count'] if result.data else 0
       
       if empty_count > 0:
           issues.append({
               'type': 'empty_mixture',
               'severity': 'warning',
               'message': f"Found {empty_count} mixtures with no components"
           })
           
       # Check for invalid concentration values
       query = """
           SELECT COUNT(*) AS count
           FROM mixture_components
           WHERE concentration < 0 OR concentration > 100
       """
       result = conn.sql(query).execute()
       invalid_count = result.data[0]['count'] if result.data else 0
       
       if invalid_count > 0:
           issues.append({
               'type': 'invalid_concentration',
               'severity': 'warning',
               'message': f"Found {invalid_count} mixture components with invalid concentration"
           })
           
       return issues
   
   def _check_duplicate_records(conn: Any) -> List[Dict]:
       """
       Check for duplicate records.
       
       Args:
           conn: Database connection
           
       Returns:
           List of issues found
       """
       issues = []
       
       # Check for duplicate molecule names
       query = """
           SELECT name, COUNT(*) AS count
           FROM molecules
           GROUP BY name
           HAVING COUNT(*) > 1
       """
       result = conn.sql(query).execute()
       
       if result.data:
           total_duplicates = sum(row['count'] - 1 for row in result.data)
           issues.append({
               'type': 'duplicate_molecule_name',
               'severity': 'warning',
               'message': f"Found {total_duplicates} duplicate molecule names across {len(result.data)} groups"
           })
           
       # Check for duplicate mixture names
       query = """
           SELECT name, COUNT(*) AS count
           FROM mixtures
           GROUP BY name
           HAVING COUNT(*) > 1
       """
       result = conn.sql(query).execute()
       
       if result.data:
           total_duplicates = sum(row['count'] - 1 for row in result.data)
           issues.append({
               'type': 'duplicate_mixture_name',
               'severity': 'warning',
               'message': f"Found {total_duplicates} duplicate mixture names across {len(result.data)} groups"
           })
           
       return issues
   
   def run_data_quality_verification(conn: Any, level: str = 'standard') -> Dict:
       """
       Run data quality verification checks.
       
       Args:
           conn: Database connection
           level: Verification level
           
       Returns:
           Verification results
       """
       logger.info(f"Running data quality verification at {level} level")
       
       issues = []
       
       try:
           # Basic level checks
           issues.extend(_check_required_fields(conn))
           
           # Standard and comprehensive level checks
           if level in ('standard', 'comprehensive'):
               issues.extend(_check_data_consistency(conn))
           
           # Comprehensive level checks
           if level == 'comprehensive':
               issues.extend(_check_duplicate_records(conn))
               
           success = not any(issue['severity'] == 'error' for issue in issues)
           
           return {
               'success': success,
               'issues': issues
           }
       except Exception as e:
           logger.error(f"Error during data quality verification: {str(e)}")
           return {
               'success': False,
               'issues': [
                   {
                       'type': 'verification_error',
                       'severity': 'error',
                       'message': f"Error during data quality verification: {str(e)}"
                   }
               ]
           }
   ```

7. Implement a simple HTML reporter:
   ```python
   """
   HTML report generator for verification results.
   
   This module provides functions for generating HTML reports.
   """
   
   import logging
   from datetime import datetime
   from typing import Dict
   
   logger = logging.getLogger(__name__)
   
   def generate_html_report(results: Dict) -> str:
       """
       Generate HTML report from verification results.
       
       Args:
           results: Verification results dictionary
           
       Returns:
           HTML report as a string
       """
       timestamp = results.get('timestamp', datetime.now().isoformat())
       level = results.get('level', 'unknown')
       success = results.get('success', False)
       
       html = f"""<!DOCTYPE html>
   <html>
   <head>
       <title>Database Verification Report</title>
       <style>
           body {{
               font-family: Arial, sans-serif;
               line-height: 1.6;
               margin: 0;
               padding: 20px;
               color: #333;
           }}
           h1, h2, h3 {{
               color: #444;
           }}
           .header {{
               margin-bottom: 20px;
               padding-bottom: 10px;
               border-bottom: 1px solid #eee;
           }}
           .timestamp {{
               color: #888;
               font-size: 0.9em;
           }}
           .success {{
               color: green;
           }}
           .failure {{
               color: red;
           }}
           .module {{
               margin-bottom: 20px;
               padding: 15px;
               background-color: #f9f9f9;
               border-radius: 5px;
           }}
           .issue {{
               margin: 10px 0;
               padding: 10px;
               background-color: #fff;
               border-radius: 3px;
               border-left: 4px solid #ddd;
           }}
           .issue.error {{
               border-left-color: #e74c3c;
           }}
           .issue.warning {{
               border-left-color: #f39c12;
           }}
           .summary {{
               margin-top: 20px;
               padding: 15px;
               background-color: #f0f0f0;
               border-radius: 5px;
           }}
       </style>
   </head>
   <body>
       <div class="header">
           <h1>Database Verification Report</h1>
           <div class="timestamp">Generated on {timestamp}</div>
       </div>
       
       <h2>Overview</h2>
       <p>Verification level: <strong>{level}</strong></p>
       <p>Overall status: <span class="{('success' if success else 'failure')}">
           {('PASSED' if success else 'FAILED')}
       </span></p>
       
       <h2>Module Results</h2>
   """
       
       # Add module results
       for module, module_results in results.get('results', {}).items():
           module_success = module_results.get('success', False)
           module_status = 'PASSED' if module_success else 'FAILED'
           status_class = 'success' if module_success else 'failure'
           
           html += f"""
       <div class="module">
           <h3>{module} <span class="{status_class}">{module_status}</span></h3>
   """
           
           # Add issues
           issues = module_results.get('issues', [])
           if issues:
               html += f"<p>Found {len(issues)} issues:</p>"
               
               for issue in issues:
                   severity = issue.get('severity', 'info')
                   message = issue.get('message', 'Unknown issue')
                   
                   html += f"""
           <div class="issue {severity}">
               <strong>{severity.upper()}:</strong> {message}
           </div>
   """
           else:
               html += "<p>No issues found.</p>"
               
           html += "</div>"
           
       # Add summary
       total_issues = sum(
           len(module_results.get('issues', []))
           for module_results in results.get('results', {}).values()
       )
       
       error_count = sum(
           sum(1 for issue in module_results.get('issues', []) if issue.get('severity') == 'error')
           for module_results in results.get('results', {}).values()
       )
       
       warning_count = sum(
           sum(1 for issue in module_results.get('issues', []) if issue.get('severity') == 'warning')
           for module_results in results.get('results', {}).values()
       )
       
       html += f"""
       <div class="summary">
           <h2>Summary</h2>
           <p>Total issues: {total_issues}</p>
           <p>Errors: {error_count}</p>
           <p>Warnings: {warning_count}</p>
       </div>
   </body>
   </html>
   """
       
       return html
   ```

8. Create a README for the verification module:
   ```markdown
   # Database Verification Module

   This module provides tools for verifying database integrity, schema consistency, and data quality for the CryoProtect v2 project.

   ## Usage

   ### Command Line Interface

   The verification module provides a command line interface for running verification checks:

   ```bash
   # Run standard verification
   python -m database.verification

   # Run comprehensive verification
   python -m database.verification --level comprehensive

   # Run specific verification modules
   python -m database.verification --include schema data_quality

   # Generate HTML report
   python -m database.verification --format html --output report.html
   ```

   ### Programmatic Usage

   The verification module can also be used programmatically:

   ```python
   from database.verification import verify_database, verify_schema, verify_constraints, verify_data_quality

   # Run full database verification
   results = verify_database(level='standard')

   # Run specific verification
   schema_results = verify_schema()
   constraint_results = verify_constraints()
   data_quality_results = verify_data_quality()

   # Check success
   if results['success']:
       print("Verification passed!")
   else:
       print("Verification failed!")
       for module, module_results in results['results'].items():
           for issue in module_results.get('issues', []):
               print(f"{module}: {issue['message']}")
   ```

   ## Verification Levels

   The verification module supports three verification levels:

   - **basic**: Fast, essential checks only
   - **standard**: Default level, balanced thoroughness and performance
   - **comprehensive**: Thorough checks, may be slower

   ## Report Formats

   Verification results can be output in several formats:

   - **console**: Print results to the console
   - **json**: Generate a JSON report
   - **html**: Generate an HTML report

   ## Adding Custom Verifications

   To add a custom verification, create a new function in the appropriate module and add it to the verification runner.

   ```python
   # In database/verification/data.py
   def _check_custom_requirement(conn):
       issues = []
       # Implement custom check
       return issues

   # In database/verification/data.py
   def run_data_quality_verification(conn, level):
       issues = []
       # Add custom check to appropriate level
       if level == 'comprehensive':
           issues.extend(_check_custom_requirement(conn))
       # ...
   ```
   ```

9. Update the main database `__init__.py` to include the verification module:
   ```python
   """
   Database package for CryoProtect v2.
   
   This package centralizes all database-related functionality including
   models, operations, migrations, and verification tools.
   """
   
   from database.main import (
       initialize_database,
       populate_database,
       verify_database,
       backup_database,
       restore_database
   )
   
   # Include the population module exports
   from database.population import (
       populate_all,
       populate_specific,
       populate_from_file
   )
   
   # Include the migrations module exports
   from database.migrations import (
       apply_migrations,
       rollback_migrations,
       get_migration_status,
       initialize_migration_tracking
   )
   
   # Include the verification module exports
   from database.verification import (
       verify_database,
       verify_schema,
       verify_constraints,
       verify_data_quality
   )
   
   __all__ = [
       'initialize_database',
       'populate_database',
       'verify_database',
       'backup_database',
       'restore_database',
       'populate_all',
       'populate_specific',
       'populate_from_file',
       'apply_migrations',
       'rollback_migrations',
       'get_migration_status',
       'initialize_migration_tracking',
       'verify_database',
       'verify_schema',
       'verify_constraints',
       'verify_data_quality'
   ]
   ```

10. Update the main entry point to use the new verification module:
    ```python
    """
    Main entry points for the database package.
    
    This module provides high-level functions for common database operations.
    """
    
    import logging
    from typing import Dict, List, Optional, Union
    
    # Import from the population module
    from database.population.runner import populate_all, populate_specific
    
    # Import from the verification module
    from database.verification.runner import verify_database as run_verification
    
    logger = logging.getLogger(__name__)
    
    def initialize_database(config: Dict = None) -> bool:
        """
        Initialize the database with schema and required tables.
        
        Args:
            config: Configuration dictionary with connection details
            
        Returns:
            True if initialization was successful, False otherwise
        """
        logger.info("Initializing database")
        # Implementation to be added
        return True
        
    def populate_database(
        tables: Optional[List[str]] = None,
        environment: str = 'development',
        config: Optional[Dict] = None
    ) -> Dict[str, int]:
        """
        Populate the database with data.
        
        Args:
            tables: List of tables to populate, or None for all tables
            environment: Environment to use ('development', 'staging', 'production')
            config: Optional configuration dictionary
            
        Returns:
            Dictionary with table names as keys and count of inserted records as values
        """
        logger.info(f"Populating database in {environment} environment")
        
        if tables:
            return populate_specific(tables, environment, config)
        else:
            return populate_all(environment, config)
    
    def verify_database(
        level: str = 'standard',
        modules: Optional[List[str]] = None,
        config: Optional[Dict] = None
    ) -> Dict:
        """
        Verify database integrity, schema consistency, and data quality.
        
        Args:
            level: Verification level ('basic', 'standard', 'comprehensive')
            modules: Optional list of specific verification modules
            config: Optional configuration dictionary
            
        Returns:
            Verification results dictionary
        """
        logger.info(f"Verifying database at {level} level")
        
        return run_verification(
            config=config,
            level=level,
            include_modules=modules
        )
    
    # ... rest of the file unchanged
    ```

## Files to Modify

- Create new directory: `/mnt/c/Users/1edwa/documents/cryoprotect v2/database/verification/`
- Create new directory: `/mnt/c/Users/1edwa/documents/cryoprotect v2/database/verification/reporters/`
- Create new file: `/mnt/c/Users/1edwa/documents/cryoprotect v2/database/verification/__init__.py`
- Create new file: `/mnt/c/Users/1edwa/documents/cryoprotect v2/database/verification/runner.py`
- Create new file: `/mnt/c/Users/1edwa/documents/cryoprotect v2/database/verification/schema.py`
- Create new file: `/mnt/c/Users/1edwa/documents/cryoprotect v2/database/verification/constraints.py`
- Create new file: `/mnt/c/Users/1edwa/documents/cryoprotect v2/database/verification/data.py`
- Create new file: `/mnt/c/Users/1edwa/documents/cryoprotect v2/database/verification/reporters/__init__.py`
- Create new file: `/mnt/c/Users/1edwa/documents/cryoprotect v2/database/verification/reporters/html.py`
- Create new file: `/mnt/c/Users/1edwa/documents/cryoprotect v2/database/verification/README.md`
- Update existing file: `/mnt/c/Users/1edwa/documents/cryoprotect v2/database/__init__.py`
- Update existing file: `/mnt/c/Users/1edwa/documents/cryoprotect v2/database/main.py`

## Verification
1. Verify that all modules import properly without errors
2. Check that schema verification works with the database
3. Test constraint verification with sample constraints
4. Verify data quality checks run properly
5. Ensure report generation works in all formats
6. Test the CLI with various options
7. Verify integration with the main database package