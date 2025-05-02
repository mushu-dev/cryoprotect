"""
Schema Validator Module for CryoProtect v2

This module provides a class-based schema validation utility for the database health check system.
It verifies database schema integrity including tables, columns, constraints, views, and RLS policies.
"""

import logging
from typing import Dict, List, Any, Optional, Set, Tuple
import json

logger = logging.getLogger(__name__)

class SchemaValidator:
    """
    Schema validator for database health checks.
    
    This class provides methods to validate database schema integrity:
    - Table existence verification
    - Column definition and type checking
    - Constraint validation (primary keys, foreign keys, unique)
    - View definition verification
    - RLS policy verification
    
    The validator follows a class-based pattern with a run_checks method that
    accepts a connection pool and returns structured results.
    """
    
    # Class attribute for categorization in health check system
    category = "schema"
    
    def __init__(self):
        """Initialize the schema validator."""
        # Define expected tables and their columns
        # This could be loaded from a configuration file in a production system
        self.expected_tables = {
            'molecules': [
                {'name': 'id', 'type': 'uuid', 'nullable': False, 'primary_key': True},
                {'name': 'name', 'type': 'character varying', 'nullable': False},
                {'name': 'smiles', 'type': 'character varying', 'nullable': True},
                {'name': 'created_at', 'type': 'timestamp with time zone', 'nullable': True},
                {'name': 'updated_at', 'type': 'timestamp with time zone', 'nullable': True}
            ],
            'mixtures': [
                {'name': 'id', 'type': 'uuid', 'nullable': False, 'primary_key': True},
                {'name': 'name', 'type': 'character varying', 'nullable': False},
                {'name': 'description', 'type': 'text', 'nullable': True},
                {'name': 'created_at', 'type': 'timestamp with time zone', 'nullable': True},
                {'name': 'updated_at', 'type': 'timestamp with time zone', 'nullable': True},
                {'name': 'created_by', 'type': 'uuid', 'nullable': True}
            ],
            'mixture_components': [
                {'name': 'id', 'type': 'uuid', 'nullable': False, 'primary_key': True},
                {'name': 'mixture_id', 'type': 'uuid', 'nullable': False, 'foreign_key': {'table': 'mixtures', 'column': 'id'}},
                {'name': 'molecule_id', 'type': 'uuid', 'nullable': False, 'foreign_key': {'table': 'molecules', 'column': 'id'}},
                {'name': 'concentration', 'type': 'double precision', 'nullable': True},
                {'name': 'unit', 'type': 'character varying', 'nullable': True},
                {'name': 'created_at', 'type': 'timestamp with time zone', 'nullable': True},
                {'name': 'updated_at', 'type': 'timestamp with time zone', 'nullable': True}
            ],
            'experiments': [
                {'name': 'id', 'type': 'uuid', 'nullable': False, 'primary_key': True},
                {'name': 'mixture_id', 'type': 'uuid', 'nullable': False, 'foreign_key': {'table': 'mixtures', 'column': 'id'}},
                {'name': 'property_type_id', 'type': 'uuid', 'nullable': False},
                {'name': 'numeric_value', 'type': 'double precision', 'nullable': True},
                {'name': 'text_value', 'type': 'text', 'nullable': True},
                {'name': 'boolean_value', 'type': 'boolean', 'nullable': True},
                {'name': 'created_at', 'type': 'timestamp with time zone', 'nullable': True},
                {'name': 'updated_at', 'type': 'timestamp with time zone', 'nullable': True}
            ],
            'lab_verifications': [
                {'name': 'id', 'type': 'uuid', 'nullable': False, 'primary_key': True},
                {'name': 'experiment_id', 'type': 'uuid', 'nullable': False, 'foreign_key': {'table': 'experiments', 'column': 'id'}},
                {'name': 'verification_status', 'type': 'character varying', 'nullable': False},
                {'name': 'verifier', 'type': 'character varying', 'nullable': True},
                {'name': 'equipment_used', 'type': 'character varying', 'nullable': True},
                {'name': 'comments', 'type': 'text', 'nullable': True},
                {'name': 'created_at', 'type': 'timestamp with time zone', 'nullable': True},
                {'name': 'updated_at', 'type': 'timestamp with time zone', 'nullable': True}
            ]
        }
        
        # Define expected views
        self.expected_views = {
            'molecule_details': {
                'base_tables': ['molecules', 'molecular_properties'],
                'columns': ['id', 'name', 'smiles', 'molecular_weight', 'logp']
            },
            'mixture_details': {
                'base_tables': ['mixtures', 'mixture_components', 'molecules'],
                'columns': ['id', 'name', 'description', 'components']
            }
        }
        
        # Define tables that should have RLS policies
        self.tables_with_rls = [
            'molecules',
            'mixtures',
            'mixture_components',
            'experiments',
            'lab_verifications'
        ]
    
    def run_checks(self, connection_pool) -> Dict[str, Any]:
        """
        Run all schema validation checks.
        
        Args:
            connection_pool: Database connection pool
            
        Returns:
            Dict with validation results including status, issues found, and recommendations
        """
        results = {
            'status': 'passed',
            'issues_found': 0,
            'details': {},
            'recommendations': []
        }
        
        try:
            # Get a connection from the pool
            conn = connection_pool
            
            # Run all checks
            table_results = self._check_tables_exist(conn)
            column_results = self._check_column_definitions(conn)
            constraint_results = self._check_constraints(conn)
            view_results = self._check_views(conn)
            rls_results = self._check_rls_policies(conn)
            
            # Combine results
            all_checks = [
                ('tables', table_results),
                ('columns', column_results),
                ('constraints', constraint_results),
                ('views', view_results),
                ('rls_policies', rls_results)
            ]
            
            # Process results
            for check_name, check_result in all_checks:
                results['details'][check_name] = check_result
                
                # Update status based on severity
                if check_result['status'] == 'failed':
                    results['status'] = 'failed'
                elif check_result['status'] == 'warning' and results['status'] != 'failed':
                    results['status'] = 'warning'
                
                # Count issues
                results['issues_found'] += check_result['issues_found']
                
                # Add recommendations
                if 'recommendations' in check_result:
                    results['recommendations'].extend(check_result['recommendations'])
            
        except Exception as e:
            logger.error(f"Error during schema validation: {str(e)}")
            results['status'] = 'error'
            results['details']['error'] = str(e)
            results['recommendations'].append("Fix the error that prevented schema validation from completing")
        
        return results
    
    def _check_column_definitions(self, conn) -> Dict[str, Any]:
        """
        Check column definitions and types for all expected tables.
        
        Args:
            conn: Database connection
            
        Returns:
            Dict with check results
        """
        results = {
            'status': 'passed',
            'issues_found': 0,
            'column_issues': {},
            'recommendations': []
        }
        
        try:
            # For each expected table
            for table_name, expected_columns in self.expected_tables.items():
                # Skip if table doesn't exist (already reported in table check)
                table_exists_query = f"""
                    SELECT EXISTS (
                        SELECT FROM information_schema.tables
                        WHERE table_schema = 'public'
                        AND table_name = '{table_name}'
                    )
                """
                table_exists_response = conn.rpc('exec_sql', {'query': table_exists_query}).execute()
                
                if not (hasattr(table_exists_response, 'data') and
                        table_exists_response.data and
                        table_exists_response.data[0]['exists']):
                    continue
                
                # Get actual columns for this table
                column_query = f"""
                    SELECT
                        column_name,
                        data_type,
                        is_nullable,
                        column_default
                    FROM information_schema.columns
                    WHERE table_schema = 'public'
                    AND table_name = '{table_name}'
                """
                column_response = conn.rpc('exec_sql', {'query': column_query}).execute()
                
                if hasattr(column_response, 'data') and column_response.data:
                    actual_columns = {
                        col['column_name']: {
                            'type': col['data_type'],
                            'nullable': col['is_nullable'] == 'YES',
                            'default': col['column_default']
                        } for col in column_response.data
                    }
                    
                    table_issues = {
                        'missing_columns': [],
                        'type_mismatches': [],
                        'nullable_mismatches': []
                    }
                    
                    # Check each expected column
                    for expected_col in expected_columns:
                        col_name = expected_col['name']
                        
                        # Check if column exists
                        if col_name not in actual_columns:
                            table_issues['missing_columns'].append(col_name)
                            results['issues_found'] += 1
                            results['recommendations'].append(
                                f"Add missing column '{col_name}' to table '{table_name}'"
                            )
                            continue
                        
                        actual_col = actual_columns[col_name]
                        
                        # Check column type
                        if expected_col['type'] not in actual_col['type']:
                            table_issues['type_mismatches'].append({
                                'column': col_name,
                                'expected': expected_col['type'],
                                'actual': actual_col['type']
                            })
                            results['issues_found'] += 1
                            results['recommendations'].append(
                                f"Fix type mismatch for column '{col_name}' in table '{table_name}' "
                                f"(expected: {expected_col['type']}, actual: {actual_col['type']})"
                            )
                        
                        # Check nullability (only if explicitly set in expected schema)
                        if 'nullable' in expected_col and expected_col['nullable'] != actual_col['nullable']:
                            table_issues['nullable_mismatches'].append({
                                'column': col_name,
                                'expected': expected_col['nullable'],
                                'actual': actual_col['nullable']
                            })
                            results['issues_found'] += 1
                            nullable_text = "NOT NULL" if not expected_col['nullable'] else "NULL"
                            results['recommendations'].append(
                                f"Alter column '{col_name}' in table '{table_name}' to be {nullable_text}"
                            )
                    
                    # Add table issues to results if any were found
                    if (table_issues['missing_columns'] or
                        table_issues['type_mismatches'] or
                        table_issues['nullable_mismatches']):
                        results['column_issues'][table_name] = table_issues
                        results['status'] = 'failed'
        
        except Exception as e:
            logger.error(f"Error checking column definitions: {str(e)}")
            results['status'] = 'error'
            results['issues_found'] += 1
            results['error'] = str(e)
        
        return results
    
    def _check_constraints(self, conn) -> Dict[str, Any]:
        """
        Validate constraints (primary keys, foreign keys, unique).
        
        Args:
            conn: Database connection
            
        Returns:
            Dict with check results
        """
        results = {
            'status': 'passed',
            'issues_found': 0,
            'constraint_issues': {},
            'recommendations': []
        }
        
        try:
            # Check primary keys
            pk_issues = self._check_primary_keys(conn)
            if pk_issues:
                results['constraint_issues']['primary_keys'] = pk_issues
                results['issues_found'] += len(pk_issues)
                results['status'] = 'failed'
                
                # Add recommendations
                for table, columns in pk_issues.items():
                    if not columns['actual']:
                        results['recommendations'].append(
                            f"Add primary key constraint on column '{columns['expected']}' for table '{table}'"
                        )
                    else:
                        results['recommendations'].append(
                            f"Fix primary key mismatch for table '{table}' "
                            f"(expected: {columns['expected']}, actual: {columns['actual']})"
                        )
            
            # Check foreign keys
            fk_issues = self._check_foreign_keys(conn)
            if fk_issues:
                results['constraint_issues']['foreign_keys'] = fk_issues
                results['issues_found'] += len(fk_issues)
                results['status'] = 'failed'
                
                # Add recommendations
                for issue in fk_issues:
                    results['recommendations'].append(
                        f"Add missing foreign key constraint on '{issue['table']}.{issue['column']}' "
                        f"referencing '{issue['ref_table']}.{issue['ref_column']}'"
                    )
            
            # Check unique constraints
            unique_issues = self._check_unique_constraints(conn)
            if unique_issues:
                results['constraint_issues']['unique_constraints'] = unique_issues
                results['issues_found'] += len(unique_issues)
                results['status'] = 'failed'
                
                # Add recommendations
                for issue in unique_issues:
                    results['recommendations'].append(
                        f"Add missing unique constraint on column(s) '{', '.join(issue['columns'])}' "
                        f"for table '{issue['table']}'"
                    )
        
        except Exception as e:
            logger.error(f"Error checking constraints: {str(e)}")
            results['status'] = 'error'
            results['issues_found'] += 1
            results['error'] = str(e)
        
        return results
    
    def _check_primary_keys(self, conn) -> Dict[str, Dict[str, Any]]:
        """
        Check primary key constraints.
        
        Args:
            conn: Database connection
            
        Returns:
            Dict of tables with primary key issues
        """
        pk_issues = {}
        
        # Get actual primary keys
        pk_query = """
            SELECT
                tc.table_name,
                kcu.column_name
            FROM information_schema.table_constraints tc
            JOIN information_schema.key_column_usage kcu
                ON tc.constraint_name = kcu.constraint_name
                AND tc.table_schema = kcu.table_schema
            WHERE tc.constraint_type = 'PRIMARY KEY'
            AND tc.table_schema = 'public'
        """
        pk_response = conn.rpc('exec_sql', {'query': pk_query}).execute()
        
        if hasattr(pk_response, 'data') and pk_response.data:
            actual_pks = {row['table_name']: row['column_name'] for row in pk_response.data}
            
            # Check against expected primary keys
            for table_name, columns in self.expected_tables.items():
                expected_pk_cols = [col['name'] for col in columns if col.get('primary_key')]
                
                if expected_pk_cols and table_name not in actual_pks:
                    pk_issues[table_name] = {
                        'expected': expected_pk_cols[0],  # Assuming single-column PK
                        'actual': None
                    }
                elif expected_pk_cols and actual_pks[table_name] != expected_pk_cols[0]:
                    pk_issues[table_name] = {
                        'expected': expected_pk_cols[0],
                        'actual': actual_pks[table_name]
                    }
        
        return pk_issues
    
    def _check_tables_exist(self, conn) -> Dict[str, Any]:
        """
        Verify that all required tables exist in the database.
        
        Args:
            conn: Database connection
            
        Returns:
            Dict with check results
        """
        results = {
            'status': 'passed',
            'issues_found': 0,
            'missing_tables': [],
            'recommendations': []
        }
        
        try:
            # Get list of tables in the database
            query = """
                SELECT tablename 
                FROM pg_tables
                WHERE schemaname = 'public'
            """
            response = conn.rpc('exec_sql', {'query': query}).execute()
            
            if hasattr(response, 'data') and response.data:
                actual_tables = {row['tablename'] for row in response.data}
                
                # Check for missing tables
                expected_table_names = set(self.expected_tables.keys())
                missing_tables = expected_table_names - actual_tables
                
                if missing_tables:
                    results['status'] = 'failed'
                    results['issues_found'] = len(missing_tables)
                    results['missing_tables'] = list(missing_tables)
                    
                    # Add recommendations
                    for table in missing_tables:
                        results['recommendations'].append(f"Create missing table '{table}'")
            else:
                results['status'] = 'error'
                results['issues_found'] = 1
                results['error'] = "Could not retrieve table list"
                results['recommendations'].append("Check database connection and permissions")
        
        except Exception as e:
            logger.error(f"Error checking tables: {str(e)}")
            results['status'] = 'error'
            results['issues_found'] = 1
            results['error'] = str(e)
        
        return results
    
    def _check_column_definitions(self, conn) -> Dict[str, Any]:
        """
        Check column definitions and types for all expected tables.
        
        Args:
            conn: Database connection
            
        Returns:
            Dict with check results
        """
        results = {
            'status': 'passed',
            'issues_found': 0,
            'column_issues': {},
            'recommendations': []
        }
        
        try:
            # For each expected table
            for table_name, expected_columns in self.expected_tables.items():
                # Skip if table doesn't exist (already reported in table check)
                table_exists_query = f"""
                    SELECT EXISTS (
                        SELECT FROM information_schema.tables
                        WHERE table_schema = 'public'
                        AND table_name = '{table_name}'
                    )
                """
                table_exists_response = conn.rpc('exec_sql', {'query': table_exists_query}).execute()
                
                if not (hasattr(table_exists_response, 'data') and
                        table_exists_response.data and
                        table_exists_response.data[0]['exists']):
                    continue
                
                # Get actual columns for this table
                column_query = f"""
                    SELECT
                        column_name,
                        data_type,
                        is_nullable,
                        column_default
                    FROM information_schema.columns
                    WHERE table_schema = 'public'
                    AND table_name = '{table_name}'
                """
                column_response = conn.rpc('exec_sql', {'query': column_query}).execute()
                
                if hasattr(column_response, 'data') and column_response.data:
                    actual_columns = {
                        col['column_name']: {
                            'type': col['data_type'],
                            'nullable': col['is_nullable'] == 'YES',
                            'default': col['column_default']
                        } for col in column_response.data
                    }
                    
                    table_issues = {
                        'missing_columns': [],
                        'type_mismatches': [],
                        'nullable_mismatches': []
                    }
                    
                    # Check each expected column
                    for expected_col in expected_columns:
                        col_name = expected_col['name']
                        
                        # Check if column exists
                        if col_name not in actual_columns:
                            table_issues['missing_columns'].append(col_name)
                            results['issues_found'] += 1
                            results['recommendations'].append(
                                f"Add missing column '{col_name}' to table '{table_name}'"
                            )
                            continue
                        
                        actual_col = actual_columns[col_name]
                        
                        # Check column type
                        if expected_col['type'] not in actual_col['type']:
                            table_issues['type_mismatches'].append({
                                'column': col_name,
                                'expected': expected_col['type'],
                                'actual': actual_col['type']
                            })
                            results['issues_found'] += 1
                            results['recommendations'].append(
                                f"Fix type mismatch for column '{col_name}' in table '{table_name}' "
                                f"(expected: {expected_col['type']}, actual: {actual_col['type']})"
                            )
                        
                        # Check nullability (only if explicitly set in expected schema)
                        if 'nullable' in expected_col and expected_col['nullable'] != actual_col['nullable']:
                            table_issues['nullable_mismatches'].append({
                                'column': col_name,
                                'expected': expected_col['nullable'],
                                'actual': actual_col['nullable']
                            })
                            results['issues_found'] += 1
                            nullable_text = "NOT NULL" if not expected_col['nullable'] else "NULL"
                            results['recommendations'].append(
                                f"Alter column '{col_name}' in table '{table_name}' to be {nullable_text}"
                            )
                    
                    # Add table issues to results if any were found
                    if (table_issues['missing_columns'] or
                        table_issues['type_mismatches'] or
                        table_issues['nullable_mismatches']):
                        results['column_issues'][table_name] = table_issues
                        results['status'] = 'failed'
        
        except Exception as e:
            logger.error(f"Error checking column definitions: {str(e)}")
            results['status'] = 'error'
            results['issues_found'] += 1
            results['error'] = str(e)
        
        return results
    
    def _check_constraints(self, conn) -> Dict[str, Any]:
        """
        Validate constraints (primary keys, foreign keys, unique).
        
        Args:
            conn: Database connection
            
        Returns:
            Dict with check results
        """
        results = {
            'status': 'passed',
            'issues_found': 0,
            'constraint_issues': {},
            'recommendations': []
        }
        
        try:
            # Check primary keys
            pk_issues = self._check_primary_keys(conn)
            if pk_issues:
                results['constraint_issues']['primary_keys'] = pk_issues
                results['issues_found'] += len(pk_issues)
                results['status'] = 'failed'
                
                # Add recommendations
                for table, columns in pk_issues.items():
                    if not columns['actual']:
                        results['recommendations'].append(
                            f"Add primary key constraint on column '{columns['expected']}' for table '{table}'"
                        )
                    else:
                        results['recommendations'].append(
                            f"Fix primary key mismatch for table '{table}' "
                            f"(expected: {columns['expected']}, actual: {columns['actual']})"
                        )
            
            # Check foreign keys
            fk_issues = self._check_foreign_keys(conn)
            if fk_issues:
                results['constraint_issues']['foreign_keys'] = fk_issues
                results['issues_found'] += len(fk_issues)
                results['status'] = 'failed'
                
                # Add recommendations
                for issue in fk_issues:
                    results['recommendations'].append(
                        f"Add missing foreign key constraint on '{issue['table']}.{issue['column']}' "
                        f"referencing '{issue['ref_table']}.{issue['ref_column']}'"
                    )
            
            # Check unique constraints
            unique_issues = self._check_unique_constraints(conn)
            if unique_issues:
                results['constraint_issues']['unique_constraints'] = unique_issues
                results['issues_found'] += len(unique_issues)
                results['status'] = 'failed'
                
                # Add recommendations
                for issue in unique_issues:
                    results['recommendations'].append(
                        f"Add missing unique constraint on column(s) '{', '.join(issue['columns'])}' "
                        f"for table '{issue['table']}'"
                    )
        
        except Exception as e:
            logger.error(f"Error checking constraints: {str(e)}")
            results['status'] = 'error'
            results['issues_found'] += 1
            results['error'] = str(e)
        
        return results
    
    def _check_primary_keys(self, conn) -> Dict[str, Dict[str, Any]]:
        """
        Check primary key constraints.
        
        Args:
            conn: Database connection
            
        Returns:
            Dict of tables with primary key issues
        """
        pk_issues = {}
        
        # Get actual primary keys
        pk_query = """
            SELECT
                tc.table_name,
                kcu.column_name
            FROM information_schema.table_constraints tc
            JOIN information_schema.key_column_usage kcu
                ON tc.constraint_name = kcu.constraint_name
                AND tc.table_schema = kcu.table_schema
            WHERE tc.constraint_type = 'PRIMARY KEY'
            AND tc.table_schema = 'public'
        """
        pk_response = conn.rpc('exec_sql', {'query': pk_query}).execute()
        
        if hasattr(pk_response, 'data') and pk_response.data:
            actual_pks = {row['table_name']: row['column_name'] for row in pk_response.data}
            
            # Check against expected primary keys
            for table_name, columns in self.expected_tables.items():
                expected_pk_cols = [col['name'] for col in columns if col.get('primary_key')]
                
                if expected_pk_cols and table_name not in actual_pks:
                    pk_issues[table_name] = {
                        'expected': expected_pk_cols[0],  # Assuming single-column PK
                        'actual': None
                    }
                elif expected_pk_cols and actual_pks[table_name] != expected_pk_cols[0]:
                    pk_issues[table_name] = {
                        'expected': expected_pk_cols[0],
                        'actual': actual_pks[table_name]
                    }
        
        return pk_issues
    
    def _check_foreign_keys(self, conn) -> List[Dict[str, str]]:
        """
        Check foreign key constraints.
        
        Args:
            conn: Database connection
            
        Returns:
            List of missing foreign key issues
        """
        fk_issues = []
        
        # Get actual foreign keys
        fk_query = """
            SELECT
                tc.table_name,
                kcu.column_name,
                ccu.table_name AS foreign_table_name,
                ccu.column_name AS foreign_column_name
            FROM information_schema.table_constraints tc
            JOIN information_schema.key_column_usage kcu
                ON tc.constraint_name = kcu.constraint_name
                AND tc.table_schema = kcu.table_schema
            JOIN information_schema.constraint_column_usage ccu
                ON ccu.constraint_name = tc.constraint_name
                AND ccu.table_schema = tc.table_schema
            WHERE tc.constraint_type = 'FOREIGN KEY'
            AND tc.table_schema = 'public'
        """
        fk_response = conn.rpc('exec_sql', {'query': fk_query}).execute()
        
        if hasattr(fk_response, 'data') and fk_response.data:
            actual_fks = {
                (row['table_name'], row['column_name']):
                (row['foreign_table_name'], row['foreign_column_name'])
                for row in fk_response.data
            }
            
            # Check against expected foreign keys
            for table_name, columns in self.expected_tables.items():
                for col in columns:
                    if 'foreign_key' in col:
                        fk_def = col['foreign_key']
                        expected_fk = (fk_def['table'], fk_def['column'])
                        
                        if (table_name, col['name']) not in actual_fks:
                            fk_issues.append({
                                'table': table_name,
                                'column': col['name'],
                                'ref_table': fk_def['table'],
                                'ref_column': fk_def['column']
                            })
                        elif actual_fks[(table_name, col['name'])] != expected_fk:
                            fk_issues.append({
                                'table': table_name,
                                'column': col['name'],
                                'ref_table': fk_def['table'],
                                'ref_column': fk_def['column'],
                                'actual_ref': actual_fks[(table_name, col['name'])]
                            })
        
        return fk_issues
    
    def _check_unique_constraints(self, conn) -> List[Dict[str, Any]]:
        """
        Check unique constraints.
        
        Args:
            conn: Database connection
            
        Returns:
            List of missing unique constraint issues
        """
        unique_issues = []
        
        # Get actual unique constraints
        unique_query = """
            SELECT
                tc.table_name,
                kcu.column_name
            FROM information_schema.table_constraints tc
            JOIN information_schema.key_column_usage kcu
                ON tc.constraint_name = kcu.constraint_name
                AND tc.table_schema = kcu.table_schema
            WHERE tc.constraint_type = 'UNIQUE'
            AND tc.table_schema = 'public'
        """
        unique_response = conn.rpc('exec_sql', {'query': unique_query}).execute()
        
        if hasattr(unique_response, 'data') and unique_response.data:
            # Group by table name
            actual_uniques = {}
            for row in unique_response.data:
                table = row['table_name']
                column = row['column_name']
                if table not in actual_uniques:
                    actual_uniques[table] = []
                actual_uniques[table].append(column)
            
            # Check against expected unique constraints
            # For this example, we'll just check if 'name' columns have unique constraints
            # In a real implementation, you would define expected unique constraints in self.expected_tables
            for table_name, columns in self.expected_tables.items():
                name_columns = [col['name'] for col in columns if col['name'] == 'name']
                
                if name_columns and (
                    table_name not in actual_uniques or
                    'name' not in actual_uniques[table_name]
                ):
                    unique_issues.append({
                        'table': table_name,
                        'columns': ['name']
                    })
        
        return unique_issues
    
    def _check_views(self, conn) -> Dict[str, Any]:
        """
        Ensure views are properly defined.
        
        Args:
            conn: Database connection
            
        Returns:
            Dict with check results
        """
        results = {
            'status': 'passed',
            'issues_found': 0,
            'view_issues': {},
            'recommendations': []
        }
        
        try:
            # Get list of views in the database
            view_query = """
                SELECT table_name
                FROM information_schema.views
                WHERE table_schema = 'public'
            """
            view_response = conn.rpc('exec_sql', {'query': view_query}).execute()
            
            if hasattr(view_response, 'data') and view_response.data:
                actual_views = {row['table_name'] for row in view_response.data}
                
                # Check for missing views
                expected_view_names = set(self.expected_views.keys())
                missing_views = expected_view_names - actual_views
                
                if missing_views:
                    results['status'] = 'warning'  # Views are often optional, so warning not failure
                    results['issues_found'] += len(missing_views)
                    results['view_issues']['missing_views'] = list(missing_views)
                    
                    # Add recommendations
                    for view in missing_views:
                        view_def = self.expected_views[view]
                        results['recommendations'].append(
                            f"Create missing view '{view}' based on tables: {', '.join(view_def['base_tables'])}"
                        )
                
                # Check view definitions for existing views
                for view_name in expected_view_names & actual_views:
                    expected_view = self.expected_views[view_name]
                    
                    # Get view definition
                    view_def_query = f"""
                        SELECT column_name
                        FROM information_schema.columns
                        WHERE table_schema = 'public'
                        AND table_name = '{view_name}'
                    """
                    view_def_response = conn.rpc('exec_sql', {'query': view_def_query}).execute()
                    
                    if hasattr(view_def_response, 'data') and view_def_response.data:
                        actual_columns = {row['column_name'] for row in view_def_response.data}
                        expected_columns = set(expected_view['columns'])
                        
                        missing_columns = expected_columns - actual_columns
                        
                        if missing_columns:
                            if 'incorrect_views' not in results['view_issues']:
                                results['view_issues']['incorrect_views'] = {}
                            
                            results['view_issues']['incorrect_views'][view_name] = {
                                'missing_columns': list(missing_columns)
                            }
                            results['issues_found'] += len(missing_columns)
                            results['status'] = 'warning'
                            
                            results['recommendations'].append(
                                f"Update view '{view_name}' to include missing columns: {', '.join(missing_columns)}"
                            )
        
        except Exception as e:
            logger.error(f"Error checking views: {str(e)}")
            results['status'] = 'error'
            results['issues_found'] += 1
            results['error'] = str(e)
        
        return results
    
    def _check_rls_policies(self, conn) -> Dict[str, Any]:
        """
        Verify RLS policies are in place.
        
        Args:
            conn: Database connection
            
        Returns:
            Dict with check results
        """
        results = {
            'status': 'passed',
            'issues_found': 0,
            'rls_issues': {},
            'recommendations': []
        }
        
        try:
            # Check if RLS is enabled for each table
            for table_name in self.tables_with_rls:
                # First check if table exists
                table_exists_query = f"""
                    SELECT EXISTS (
                        SELECT FROM information_schema.tables
                        WHERE table_schema = 'public'
                        AND table_name = '{table_name}'
                    )
                """
                table_exists_response = conn.rpc('exec_sql', {'query': table_exists_query}).execute()
                
                if not (hasattr(table_exists_response, 'data') and
                        table_exists_response.data and
                        table_exists_response.data[0]['exists']):
                    continue
                
                # Check if RLS is enabled
                rls_query = f"""
                    SELECT relrowsecurity
                    FROM pg_class
                    WHERE relname = '{table_name}'
                    AND relnamespace = (SELECT oid FROM pg_namespace WHERE nspname = 'public')
                """
                rls_response = conn.rpc('exec_sql', {'query': rls_query}).execute()
                
                if hasattr(rls_response, 'data') and rls_response.data:
                    rls_enabled = rls_response.data[0]['relrowsecurity']
                    
                    if not rls_enabled:
                        if 'tables_without_rls' not in results['rls_issues']:
                            results['rls_issues']['tables_without_rls'] = []
                        
                        results['rls_issues']['tables_without_rls'].append(table_name)
                        results['issues_found'] += 1
                        results['status'] = 'failed'
                        
                        results['recommendations'].append(
                            f"Enable Row Level Security for table '{table_name}'"
                        )
                
                # Check if policies exist
                policy_query = f"""
                    SELECT polname
                    FROM pg_policy
                    WHERE polrelid = '{table_name}'::regclass
                """
                policy_response = conn.rpc('exec_sql', {'query': policy_query}).execute()
                
                if hasattr(policy_response, 'data'):
                    policies = [row['polname'] for row in policy_response.data] if policy_response.data else []
                    
                    if not policies:
                        if 'tables_without_policies' not in results['rls_issues']:
                            results['rls_issues']['tables_without_policies'] = []
                        
                        results['rls_issues']['tables_without_policies'].append(table_name)
                        results['issues_found'] += 1
                        results['status'] = 'failed'
                        
                        results['recommendations'].append(
                            f"Create RLS policies for table '{table_name}'"
                        )
        
        except Exception as e:
            logger.error(f"Error checking RLS policies: {str(e)}")
            results['status'] = 'error'
            results['issues_found'] += 1
            results['error'] = str(e)
        
        return results


# Example usage:
# from database.utils.connection import initialize_connection_pool
#
# # Initialize connection pool
# initialize_connection_pool()
#
# # Create validator and run checks
# validator = SchemaValidator()
# results = validator.run_checks(connection_pool)
#
# # Process results
# if results['status'] == 'passed':
#     print("Schema validation passed!")
