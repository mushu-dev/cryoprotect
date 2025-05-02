"""
Integrity Checker Module for CryoProtect v2

This module provides a class-based data integrity verification utility for the database health check system.
It verifies database data integrity including foreign key constraints, unique constraints, data consistency,
required field validation, and potential data corruption detection.
"""

import logging
from typing import Dict, List, Any, Optional, Set, Tuple
import json
import traceback

logger = logging.getLogger(__name__)

class IntegrityChecker:
    """
    Data integrity checker for database health checks.
    
    This class provides methods to validate database data integrity:
    - Orphaned records detection (foreign key constraint violations)
    - Unique constraint validation
    - Data consistency across related tables
    - Required field validation
    - Data corruption detection
    
    The checker follows a class-based pattern with a run_checks method that
    accepts a connection pool and returns structured results.
    """
    
    # Class attribute for categorization in health check system
    category = "integrity"
    
    def __init__(self):
        """Initialize the integrity checker."""
        # Define table relationships for integrity checks
        self.table_relationships = {
            'molecules': {
                'dependent_tables': [
                    {'table': 'mixture_components', 'fk_column': 'molecule_id'},
                    {'table': 'molecular_properties', 'fk_column': 'molecule_id'},
                    {'table': 'predictions', 'fk_column': 'molecule_id'}
                ]
            },
            'mixtures': {
                'dependent_tables': [
                    {'table': 'mixture_components', 'fk_column': 'mixture_id'},
                    {'table': 'predictions', 'fk_column': 'mixture_id'},
                    {'table': 'experiments', 'fk_column': 'mixture_id'}
                ]
            },
            'experiments': {
                'dependent_tables': [
                    {'table': 'lab_verifications', 'fk_column': 'experiment_id'}
                ]
            },
            'users': {
                'dependent_tables': [
                    {'table': 'user_team_membership', 'fk_column': 'user_id'},
                    {'table': 'shared_items', 'fk_column': 'owner_id'}
                ]
            }
        }
        
        # Define tables with unique constraints
        self.unique_constraints = {
            'molecules': [
                ['smiles']  # SMILES should be unique
            ],
            'mixtures': [
                ['name', 'created_by']  # Name should be unique per user
            ],
            'users': [
                ['email']  # Email should be unique
            ]
        }
        
        # Define required fields (non-nullable fields that should have values)
        self.required_fields = {
            'molecules': ['name'],
            'mixtures': ['name'],
            'mixture_components': ['mixture_id', 'molecule_id'],
            'experiments': ['mixture_id', 'property_type_id'],
            'lab_verifications': ['experiment_id', 'verification_status']
        }
        
        # Define data consistency rules
        self.consistency_rules = [
            {
                'description': 'Mixture components concentration should be positive',
                'table': 'mixture_components',
                'condition': 'concentration <= 0 OR concentration IS NULL',
                'severity': 'warning'
            },
            {
                'description': 'Experiment numeric values should be within reasonable range',
                'table': 'experiments',
                'condition': 'numeric_value IS NOT NULL AND (numeric_value < -1000 OR numeric_value > 1000)',
                'severity': 'warning'
            }
        ]
    
    def run_checks(self, connection_pool) -> Dict[str, Any]:
        """
        Run all data integrity checks.
        
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
            orphaned_results = self._check_orphaned_records(conn)
            unique_results = self._check_unique_constraints(conn)
            consistency_results = self._check_data_consistency(conn)
            required_results = self._check_required_fields(conn)
            corruption_results = self._check_data_corruption(conn)
            
            # Combine results
            all_checks = [
                ('orphaned_records', orphaned_results),
                ('unique_constraints', unique_results),
                ('data_consistency', consistency_results),
                ('required_fields', required_results),
                ('data_corruption', corruption_results)
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
            logger.error(f"Error during integrity checks: {str(e)}")
            results['status'] = 'error'
            results['details']['error'] = str(e)
            results['details']['traceback'] = traceback.format_exc()
            results['recommendations'].append("Fix the error that prevented integrity checks from completing")
        
        return results
    
    def _check_orphaned_records(self, conn) -> Dict[str, Any]:
        """
        Check for orphaned records that violate foreign key constraints.
        
        Args:
            conn: Database connection
            
        Returns:
            Dict with check results
        """
        results = {
            'status': 'passed',
            'issues_found': 0,
            'orphaned_records': {},
            'recommendations': []
        }
        
        try:
            for parent_table, relationship in self.table_relationships.items():
                for dependent in relationship['dependent_tables']:
                    dependent_table = dependent['table']
                    fk_column = dependent['fk_column']
                    
                    # Check if both tables exist before running the query
                    tables_exist_query = f"""
                        SELECT 
                            (SELECT EXISTS(SELECT 1 FROM information_schema.tables 
                                          WHERE table_schema = 'public' AND table_name = '{parent_table}')) AS parent_exists,
                            (SELECT EXISTS(SELECT 1 FROM information_schema.tables 
                                          WHERE table_schema = 'public' AND table_name = '{dependent_table}')) AS dependent_exists
                    """
                    tables_exist_response = conn.rpc('exec_sql', {'query': tables_exist_query}).execute()
                    
                    if not (hasattr(tables_exist_response, 'data') and 
                            tables_exist_response.data and 
                            tables_exist_response.data[0]['parent_exists'] and 
                            tables_exist_response.data[0]['dependent_exists']):
                        # Skip if either table doesn't exist
                        continue
                    
                    # Check for orphaned records (foreign keys that don't exist in parent table)
                    orphaned_query = f"""
                        SELECT d.id, d.{fk_column}
                        FROM {dependent_table} d
                        LEFT JOIN {parent_table} p ON d.{fk_column} = p.id
                        WHERE p.id IS NULL;
                    """
                    
                    orphaned_response = conn.rpc('exec_sql', {'query': orphaned_query}).execute()
                    
                    if hasattr(orphaned_response, 'data') and orphaned_response.data:
                        orphaned_records = orphaned_response.data
                        if orphaned_records:
                            relation_key = f"{dependent_table}_to_{parent_table}"
                            results['orphaned_records'][relation_key] = orphaned_records
                            results['status'] = 'failed'
                            results['issues_found'] += len(orphaned_records)
                            
                            # Add recommendations
                            results['recommendations'].append(
                                f"Fix orphaned records in '{dependent_table}' that reference non-existent '{parent_table}' records"
                            )
                            
                            # Add specific recommendation for fixing
                            if len(orphaned_records) <= 10:  # Only suggest specific fixes for a reasonable number of records
                                for record in orphaned_records:
                                    results['recommendations'].append(
                                        f"Either delete record with ID '{record['id']}' from '{dependent_table}' or "
                                        f"create a matching record in '{parent_table}' with ID '{record[fk_column]}'"
                                    )
        
        except Exception as e:
            logger.error(f"Error checking orphaned records: {str(e)}")
            results['status'] = 'error'
            results['issues_found'] += 1
            results['error'] = str(e)
        
        return results
    
    def _check_unique_constraints(self, conn) -> Dict[str, Any]:
        """
        Verify that unique constraints are respected.
        
        Args:
            conn: Database connection
            
        Returns:
            Dict with check results
        """
        results = {
            'status': 'passed',
            'issues_found': 0,
            'duplicate_records': {},
            'recommendations': []
        }
        
        try:
            for table, constraints in self.unique_constraints.items():
                # Check if table exists
                table_exists_query = f"""
                    SELECT EXISTS(SELECT 1 FROM information_schema.tables 
                                 WHERE table_schema = 'public' AND table_name = '{table}')
                """
                table_exists_response = conn.rpc('exec_sql', {'query': table_exists_query}).execute()
                
                if not (hasattr(table_exists_response, 'data') and 
                        table_exists_response.data and 
                        table_exists_response.data[0]['exists']):
                    # Skip if table doesn't exist
                    continue
                
                for constraint_columns in constraints:
                    # Check for duplicate values
                    columns_str = ', '.join(constraint_columns)
                    group_by_str = ', '.join(constraint_columns)
                    
                    duplicate_query = f"""
                        SELECT {columns_str}, COUNT(*) as duplicate_count
                        FROM {table}
                        WHERE {' AND '.join([f"{col} IS NOT NULL" for col in constraint_columns])}
                        GROUP BY {group_by_str}
                        HAVING COUNT(*) > 1
                    """
                    
                    duplicate_response = conn.rpc('exec_sql', {'query': duplicate_query}).execute()
                    
                    if hasattr(duplicate_response, 'data') and duplicate_response.data:
                        duplicates = duplicate_response.data
                        if duplicates:
                            constraint_key = f"{table}_{'_'.join(constraint_columns)}"
                            results['duplicate_records'][constraint_key] = duplicates
                            results['status'] = 'failed'
                            results['issues_found'] += len(duplicates)
                            
                            # Add recommendations
                            columns_desc = ', '.join([f"'{col}'" for col in constraint_columns])
                            results['recommendations'].append(
                                f"Fix duplicate values for columns {columns_desc} in table '{table}'"
                            )
                            
                            # Add specific recommendation
                            for duplicate in duplicates:
                                values_desc = ', '.join([f"{col}='{duplicate[col]}'" for col in constraint_columns])
                                results['recommendations'].append(
                                    f"Resolve {duplicate['duplicate_count']} duplicate records with {values_desc} in table '{table}'"
                                )
        
        except Exception as e:
            logger.error(f"Error checking unique constraints: {str(e)}")
            results['status'] = 'error'
            results['issues_found'] += 1
            results['error'] = str(e)
        
        return results
    
    def _check_data_consistency(self, conn) -> Dict[str, Any]:
        """
        Check for data consistency across related tables.
        
        Args:
            conn: Database connection
            
        Returns:
            Dict with check results
        """
        results = {
            'status': 'passed',
            'issues_found': 0,
            'consistency_issues': {},
            'recommendations': []
        }
        
        try:
            # Check predefined consistency rules
            for rule in self.consistency_rules:
                table = rule['table']
                condition = rule['condition']
                description = rule['description']
                severity = rule['severity']
                
                # Check if table exists
                table_exists_query = f"""
                    SELECT EXISTS(SELECT 1 FROM information_schema.tables 
                                 WHERE table_schema = 'public' AND table_name = '{table}')
                """
                table_exists_response = conn.rpc('exec_sql', {'query': table_exists_query}).execute()
                
                if not (hasattr(table_exists_response, 'data') and 
                        table_exists_response.data and 
                        table_exists_response.data[0]['exists']):
                    # Skip if table doesn't exist
                    continue
                
                # Check for records that violate the consistency rule
                violation_query = f"""
                    SELECT id
                    FROM {table}
                    WHERE {condition}
                """
                
                violation_response = conn.rpc('exec_sql', {'query': violation_query}).execute()
                
                if hasattr(violation_response, 'data') and violation_response.data:
                    violations = violation_response.data
                    if violations:
                        rule_key = f"{table}_{description.replace(' ', '_').lower()}"
                        results['consistency_issues'][rule_key] = {
                            'description': description,
                            'violations': violations,
                            'severity': severity
                        }
                        
                        if severity == 'warning':
                            if results['status'] != 'failed':  # Don't downgrade from failed to warning
                                results['status'] = 'warning'
                        else:
                            results['status'] = 'failed'
                            
                        results['issues_found'] += len(violations)
                        
                        # Add recommendations
                        results['recommendations'].append(
                            f"Fix {len(violations)} {severity} consistency issues: {description} in table '{table}'"
                        )
        
        except Exception as e:
            logger.error(f"Error checking data consistency: {str(e)}")
            results['status'] = 'error'
            results['issues_found'] += 1
            results['error'] = str(e)
        
        return results
    
    def _check_required_fields(self, conn) -> Dict[str, Any]:
        """
        Validate that required fields are properly populated.
        
        Args:
            conn: Database connection
            
        Returns:
            Dict with check results
        """
        results = {
            'status': 'passed',
            'issues_found': 0,
            'missing_values': {},
            'recommendations': []
        }
        
        try:
            for table, fields in self.required_fields.items():
                # Check if table exists
                table_exists_query = f"""
                    SELECT EXISTS(SELECT 1 FROM information_schema.tables 
                                 WHERE table_schema = 'public' AND table_name = '{table}')
                """
                table_exists_response = conn.rpc('exec_sql', {'query': table_exists_query}).execute()
                
                if not (hasattr(table_exists_response, 'data') and 
                        table_exists_response.data and 
                        table_exists_response.data[0]['exists']):
                    # Skip if table doesn't exist
                    continue
                
                for field in fields:
                    # Check if field exists in the table
                    field_exists_query = f"""
                        SELECT EXISTS(SELECT 1 FROM information_schema.columns 
                                     WHERE table_schema = 'public' AND table_name = '{table}' 
                                     AND column_name = '{field}')
                    """
                    field_exists_response = conn.rpc('exec_sql', {'query': field_exists_query}).execute()
                    
                    if not (hasattr(field_exists_response, 'data') and 
                            field_exists_response.data and 
                            field_exists_response.data[0]['exists']):
                        # Skip if field doesn't exist
                        continue
                    
                    # Check for NULL or empty values
                    null_query = f"""
                        SELECT id
                        FROM {table}
                        WHERE {field} IS NULL OR {field} = ''
                    """
                    
                    null_response = conn.rpc('exec_sql', {'query': null_query}).execute()
                    
                    if hasattr(null_response, 'data') and null_response.data:
                        null_records = null_response.data
                        if null_records:
                            field_key = f"{table}_{field}"
                            results['missing_values'][field_key] = null_records
                            results['status'] = 'failed'
                            results['issues_found'] += len(null_records)
                            
                            # Add recommendations
                            results['recommendations'].append(
                                f"Populate missing values for required field '{field}' in table '{table}'"
                            )
                            
                            # Add specific recommendation if not too many records
                            if len(null_records) <= 5:
                                for record in null_records:
                                    results['recommendations'].append(
                                        f"Update record with ID '{record['id']}' in table '{table}' to provide a value for '{field}'"
                                    )
        
        except Exception as e:
            logger.error(f"Error checking required fields: {str(e)}")
            results['status'] = 'error'
            results['issues_found'] += 1
            results['error'] = str(e)
        
        return results
    
    def _check_data_corruption(self, conn) -> Dict[str, Any]:
        """
        Identify potential data corruption.
        
        Args:
            conn: Database connection
            
        Returns:
            Dict with check results
        """
        results = {
            'status': 'passed',
            'issues_found': 0,
            'corruption_issues': {},
            'recommendations': []
        }
        
        try:
            # Check for invalid JSON in JSON/JSONB columns
            json_columns_query = """
                SELECT table_name, column_name
                FROM information_schema.columns
                WHERE table_schema = 'public'
                AND data_type IN ('json', 'jsonb')
            """
            
            json_columns_response = conn.rpc('exec_sql', {'query': json_columns_query}).execute()
            
            if hasattr(json_columns_response, 'data') and json_columns_response.data:
                json_columns = json_columns_response.data
                
                for column_info in json_columns:
                    table_name = column_info['table_name']
                    column_name = column_info['column_name']
                    
                    # Check for malformed JSON (this is a PostgreSQL-specific approach)
                    # In PostgreSQL, invalid JSON would not be stored, but we can check for structural issues
                    corruption_query = f"""
                        SELECT id, {column_name}
                        FROM {table_name}
                        WHERE {column_name} IS NOT NULL
                        AND (
                            jsonb_typeof({column_name}::jsonb) = 'null'
                            OR {column_name}::text = '{{}}'
                            OR {column_name}::text = '[]'
                        )
                        LIMIT 10
                    """
                    
                    corruption_response = conn.rpc('exec_sql', {'query': corruption_query}).execute()
                    
                    if hasattr(corruption_response, 'data') and corruption_response.data:
                        corrupt_records = corruption_response.data
                        if corrupt_records:
                            column_key = f"{table_name}_{column_name}"
                            results['corruption_issues'][column_key] = corrupt_records
                            results['status'] = 'warning'  # JSON might be empty by design, so just a warning
                            results['issues_found'] += len(corrupt_records)
                            
                            # Add recommendations
                            results['recommendations'].append(
                                f"Review potentially empty or null JSON data in column '{column_name}' of table '{table_name}'"
                            )
            
            # Check for text columns with control characters or invalid UTF-8
            text_columns_query = """
                SELECT table_name, column_name
                FROM information_schema.columns
                WHERE table_schema = 'public'
                AND data_type IN ('text', 'character varying')
            """
            
            text_columns_response = conn.rpc('exec_sql', {'query': text_columns_query}).execute()
            
            if hasattr(text_columns_response, 'data') and text_columns_response.data:
                text_columns = text_columns_response.data
                
                for column_info in text_columns:
                    table_name = column_info['table_name']
                    column_name = column_info['column_name']
                    
                    # Check for control characters (ASCII < 32 except tabs, newlines, etc.)
                    corruption_query = f"""
                        SELECT id, {column_name}
                        FROM {table_name}
                        WHERE {column_name} ~ '[\\x00-\\x08\\x0B\\x0C\\x0E-\\x1F]'
                        LIMIT 10
                    """
                    
                    corruption_response = conn.rpc('exec_sql', {'query': corruption_query}).execute()
                    
                    if hasattr(corruption_response, 'data') and corruption_response.data:
                        corrupt_records = corruption_response.data
                        if corrupt_records:
                            column_key = f"{table_name}_{column_name}_control_chars"
                            results['corruption_issues'][column_key] = corrupt_records
                            results['status'] = 'warning'
                            results['issues_found'] += len(corrupt_records)
                            
                            # Add recommendations
                            results['recommendations'].append(
                                f"Clean control characters from text in column '{column_name}' of table '{table_name}'"
                            )
        
        except Exception as e:
            logger.error(f"Error checking data corruption: {str(e)}")
            results['status'] = 'error'
            results['issues_found'] += 1
            results['error'] = str(e)
        
        return results