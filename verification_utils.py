#!/usr/bin/env python3
"""
Verification Utilities for CryoProtect v2 Database Operations

This module provides data quality verification capabilities for database operations:
1. Schema validation
2. Data integrity checks
3. Count verification
4. Consistency checks
5. Relationship validation

Based on specifications in DATABASE_POPULATION_ISSUES.md (Section 5.3)
"""

import os
import time
import json
import logging
import threading
import psycopg2
from typing import Dict, Any, List, Optional, Union, Callable, Tuple, Set
from datetime import datetime, timedelta
from contextlib import contextmanager

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger('verification_utils')

# Try to import required modules, handle import errors gracefully
try:
    from db_connection_utils import ConnectionManager, safe_transaction
    from monitoring_utils import PerformanceMetrics, timed_operation
except ImportError as e:
    logger.warning(f"Could not import required modules: {str(e)}")
    ConnectionManager = None
    safe_transaction = None
    PerformanceMetrics = None
    timed_operation = None


class DataVerifier:
    """
    Verifies data quality and integrity in the database.
    
    This class provides methods to verify data quality, integrity, and
    consistency in the database. It can be used to validate imported data
    and ensure it meets the expected standards.
    """
    
    def __init__(self, verification_file: Optional[str] = None):
        """
        Initialize data verifier.
        
        Args:
            verification_file: Optional path to file for storing verification results
        """
        self.verification_file = verification_file
        self.lock = threading.RLock()
        self.metrics = PerformanceMetrics("data_verification", "monitoring/data_verification.json") if PerformanceMetrics else None
        
        # Initialize verification results
        self.verification_results = {
            "last_verification": None,
            "schema_validation": {},
            "count_verification": {},
            "integrity_checks": {},
            "consistency_checks": {},
            "relationship_validation": {},
            "overall_status": "Unknown"
        }
        if verification_file and os.path.exists(verification_file):
            try:
                with open(verification_file, 'r') as f:
                    saved_results = json.load(f)
                    self.verification_results.update(saved_results)
                    logger.info(f"Loaded verification results from {verification_file}")
            except Exception as e:
                logger.warning(f"Could not load verification results from {verification_file}: {str(e)}")
        
        logger.info("DataVerifier initialized")
    
    def _save_results(self) -> None:
        """Save verification results to file."""
        if not self.verification_file:
            return
            
        try:
            # Create directory if it doesn't exist
            verification_dir = os.path.dirname(self.verification_file)
            if verification_dir and not os.path.exists(verification_dir):
                os.makedirs(verification_dir)
                
            with open(self.verification_file, 'w') as f:
                json.dump(self.verification_results, f, indent=2)
                logger.info(f"Saved verification results to {self.verification_file}")
        except Exception as e:
            logger.warning(f"Could not save verification results to {self.verification_file}: {str(e)}")
    
    def verify_schema(self, expected_tables: List[Dict[str, Any]]) -> Dict[str, Any]:
        """
        Verify database schema against expected tables.
        
        Args:
            expected_tables: List of expected tables with their columns
                Each table should be a dictionary with:
                - name: Table name
                - columns: List of column dictionaries with name, type, nullable
                
        Returns:
            Dictionary with schema validation results
        """
        if not ConnectionManager:
            logger.error("ConnectionManager not available, cannot verify schema")
            return {"status": "Error", "message": "ConnectionManager not available"}
        
        start_time = time.time()
        connection = None
        results = {
            "status": "Unknown",
            "timestamp": datetime.now().isoformat(),
            "duration": 0,
            "tables_checked": 0,
            "tables_missing": [],
            "tables_extra": [],
            "column_issues": [],
            "details": {}
        }
        
        try:
            # Get database connection
            connection = ConnectionManager.get_instance().get_connection()
            if not connection:
                raise Exception("Could not get database connection")
            
            # Get actual tables and columns
            actual_tables = {}
            
            # Query to get all tables
            tables_query = """
            SELECT table_name
            FROM information_schema.tables
            WHERE table_schema = 'public'
            """
            tables_result = connection.execute_query(tables_query)
            
            # Check each table
            for table_row in tables_result:
                table_name = table_row['table_name']
                
                # Query to get columns for this table
                columns_query = """
                SELECT column_name, data_type, is_nullable
                FROM information_schema.columns
                WHERE table_schema = 'public' AND table_name = %s
                """
                columns_result = connection.execute_query(columns_query, [table_name])
                
                # Store columns
                actual_tables[table_name] = {
                    "columns": [
                        {
                            "name": col['column_name'],
                            "type": col['data_type'],
                            "nullable": col['is_nullable'] == 'YES'
                        }
                        for col in columns_result
                    ]
                }
            
            # Convert expected tables to dictionary for easier comparison
            expected_tables_dict = {table['name']: table for table in expected_tables}
            
            # Check for missing and extra tables
            expected_table_names = set(expected_tables_dict.keys())
            actual_table_names = set(actual_tables.keys())
            
            missing_tables = expected_table_names - actual_table_names
            extra_tables = actual_table_names - expected_table_names
            
            results["tables_missing"] = list(missing_tables)
            results["tables_extra"] = list(extra_tables)
            results["tables_checked"] = len(expected_tables)
            
            # Check columns for each expected table
            for table_name, expected_table in expected_tables_dict.items():
                if table_name in actual_tables:
                    actual_columns = actual_tables[table_name]["columns"]
                    expected_columns = expected_table["columns"]
                    
                    # Convert to dictionaries for easier comparison
                    actual_columns_dict = {col['name']: col for col in actual_columns}
                    expected_columns_dict = {col['name']: col for col in expected_columns}
                    
                    # Check for missing and extra columns
                    expected_column_names = set(expected_columns_dict.keys())
                    actual_column_names = set(actual_columns_dict.keys())
                    
                    missing_columns = expected_column_names - actual_column_names
                    extra_columns = actual_column_names - expected_column_names
                    
                    # Check column types and nullability
                    type_mismatches = []
                    nullable_mismatches = []
                    
                    for col_name in expected_column_names & actual_column_names:
                        expected_col = expected_columns_dict[col_name]
                        actual_col = actual_columns_dict[col_name]
                        
                        # Check type (simplified, might need more sophisticated comparison)
                        if expected_col['type'].lower() != actual_col['type'].lower():
                            type_mismatches.append({
                                "column": col_name,
                                "expected_type": expected_col['type'],
                                "actual_type": actual_col['type']
                            })
                        
                        # Check nullability
                        if expected_col['nullable'] != actual_col['nullable']:
                            nullable_mismatches.append({
                                "column": col_name,
                                "expected_nullable": expected_col['nullable'],
                                "actual_nullable": actual_col['nullable']
                            })
                    
                    # Store results for this table
                    table_results = {
                        "missing_columns": list(missing_columns),
                        "extra_columns": list(extra_columns),
                        "type_mismatches": type_mismatches,
                        "nullable_mismatches": nullable_mismatches,
                        "status": "OK"
                    }
                    
                    # Determine table status
                    if missing_columns or type_mismatches or nullable_mismatches:
                        table_results["status"] = "Error"
                        
                        # Add to column issues
                        if missing_columns:
                            results["column_issues"].append({
                                "table": table_name,
                                "issue": "missing_columns",
                                "columns": list(missing_columns)
                            })
                        
                        if type_mismatches:
                            results["column_issues"].append({
                                "table": table_name,
                                "issue": "type_mismatches",
                                "mismatches": type_mismatches
                            })
                        
                        if nullable_mismatches:
                            results["column_issues"].append({
                                "table": table_name,
                                "issue": "nullable_mismatches",
                                "mismatches": nullable_mismatches
                            })
                    
                    results["details"][table_name] = table_results
            
            # Determine overall status
            if missing_tables or results["column_issues"]:
                results["status"] = "Error"
            else:
                results["status"] = "OK"
            
        except Exception as e:
            logger.error(f"Error verifying schema: {str(e)}")
            results["status"] = "Error"
            results["message"] = str(e)
        finally:
            if connection and hasattr(connection, 'close'):
                connection.close()
            
            # Record duration
            results["duration"] = time.time() - start_time
            
            # Update verification results
            with self.lock:
                self.verification_results["schema_validation"] = results
                self.verification_results["last_verification"] = datetime.now().isoformat()
                self._update_overall_status()
                self._save_results()
            
            # Record metrics
            if self.metrics:
                self.metrics.record_operation(
                    operation_type="schema_verification",
                    success=results["status"] == "OK",
                    execution_time=results["duration"],
                    items_processed=results["tables_checked"]
                )
        
        return results
        self.verification_file = verification_file
        self.lock = threading.RLock()
        self.metrics = PerformanceMetrics("data_verification", "monitoring/data_verification.json") if PerformanceMetrics else None
        
        # Initialize verification results
        self.verification_results = {
            "last_verification": None,
            "schema_validation": {},
            "count_verification": {},
            "integrity_checks": {},
            "consistency_checks": {},
            "relationship_validation": {},
            "overall_status": "Unknown"
        }
def verify_counts(self, expected_counts: Dict[str, int]) -> Dict[str, Any]:
        """
        Verify record counts in tables against expected counts.
        
        Args:
            expected_counts: Dictionary mapping table names to expected record counts
                
        Returns:
            Dictionary with count verification results
        """
        if not ConnectionManager:
            logger.error("ConnectionManager not available, cannot verify counts")
            return {"status": "Error", "message": "ConnectionManager not available"}
        
        start_time = time.time()
        connection = None
        results = {
            "status": "Unknown",
            "timestamp": datetime.now().isoformat(),
            "duration": 0,
            "tables_checked": 0,
            "count_mismatches": [],
            "details": {}
        }
        
        try:
            # Get database connection
            connection = ConnectionManager.get_instance().get_connection()
            if not connection:
                raise Exception("Could not get database connection")
            
            # Check count for each table
            for table_name, expected_count in expected_counts.items():
                # Query to get actual count
                count_query = f"SELECT COUNT(*) as count FROM {table_name}"
                try:
                    count_result = connection.execute_query(count_query)
                    actual_count = count_result[0]['count'] if count_result else 0
                    
                    # Store results for this table
                    table_results = {
                        "expected_count": expected_count,
                        "actual_count": actual_count,
                        "difference": actual_count - expected_count,
                        "status": "OK"
                    }
                    
                    # Check if count matches
                    if actual_count != expected_count:
                        table_results["status"] = "Error"
                        results["count_mismatches"].append({
                            "table": table_name,
                            "expected": expected_count,
                            "actual": actual_count,
                            "difference": actual_count - expected_count
                        })
                    
                    results["details"][table_name] = table_results
                    results["tables_checked"] += 1
                    
                except Exception as e:
                    logger.warning(f"Error checking count for table {table_name}: {str(e)}")
                    results["details"][table_name] = {
                        "expected_count": expected_count,
                        "actual_count": None,
                        "error": str(e),
                        "status": "Error"
                    }
                    results["count_mismatches"].append({
                        "table": table_name,
                        "error": str(e)
                    })
            
            # Determine overall status
            if results["count_mismatches"]:
                results["status"] = "Error"
            else:
                results["status"] = "OK"
            
        except Exception as e:
            logger.error(f"Error verifying counts: {str(e)}")
            results["status"] = "Error"
            results["message"] = str(e)
        finally:
            if connection and hasattr(connection, 'close'):
                connection.close()
            
            # Record duration
            results["duration"] = time.time() - start_time
            
            # Update verification results
            with self.lock:
                self.verification_results["count_verification"] = results
                self.verification_results["last_verification"] = datetime.now().isoformat()
                self._update_overall_status()
                self._save_results()
            
            # Record metrics
            if self.metrics:
                self.metrics.record_operation(
                    operation_type="count_verification",
                    success=results["status"] == "OK",
                    execution_time=results["duration"],
                    items_processed=results["tables_checked"]
                )
        
        return results
def verify_integrity(self, integrity_checks: List[Dict[str, Any]]) -> Dict[str, Any]:
        """
        Verify data integrity using custom SQL checks.
        
        Args:
            integrity_checks: List of integrity check dictionaries with:
                - name: Name of the check
                - description: Description of what the check verifies
                - query: SQL query that should return 0 rows if data is valid
                - severity: 'critical', 'warning', or 'info'
                
        Returns:
            Dictionary with integrity check results
        """
        if not ConnectionManager:
            logger.error("ConnectionManager not available, cannot verify integrity")
            return {"status": "Error", "message": "ConnectionManager not available"}
        
        start_time = time.time()
        connection = None
        results = {
            "status": "Unknown",
            "timestamp": datetime.now().isoformat(),
            "duration": 0,
            "checks_performed": 0,
            "checks_failed": 0,
            "critical_failures": 0,
            "details": []
        }
        
        try:
            # Get database connection
            connection = ConnectionManager.get_instance().get_connection()
            if not connection:
                raise Exception("Could not get database connection")
            
            # Run each integrity check
            for check in integrity_checks:
                check_name = check.get('name', 'Unnamed check')
                check_description = check.get('description', '')
                check_query = check.get('query', '')
                check_severity = check.get('severity', 'warning')
                
                if not check_query:
                    logger.warning(f"Skipping integrity check '{check_name}' with empty query")
                    continue
                
                check_result = {
                    "name": check_name,
                    "description": check_description,
                    "severity": check_severity,
                    "status": "Unknown",
                    "violations": 0,
                    "sample_violations": []
                }
                
                try:
                    # Execute the check query
                    query_result = connection.execute_query(check_query)
                    violations = len(query_result) if query_result else 0
                    
                    # Store results
                    check_result["violations"] = violations
                    check_result["status"] = "OK" if violations == 0 else "Error"
                    
                    # Store sample violations (up to 10)
                    if violations > 0:
                        check_result["sample_violations"] = query_result[:10]
                        results["checks_failed"] += 1
                        
                        if check_severity == 'critical':
                            results["critical_failures"] += 1
                    
                except Exception as e:
                    logger.warning(f"Error running integrity check '{check_name}': {str(e)}")
                    check_result["status"] = "Error"
                    check_result["error"] = str(e)
                    results["checks_failed"] += 1
                
                results["details"].append(check_result)
                results["checks_performed"] += 1
            
            # Determine overall status
            if results["critical_failures"] > 0:
                results["status"] = "Critical"
            elif results["checks_failed"] > 0:
                results["status"] = "Warning"
            else:
                results["status"] = "OK"
            
        except Exception as e:
            logger.error(f"Error verifying integrity: {str(e)}")
            results["status"] = "Error"
            results["message"] = str(e)
        finally:
            if connection and hasattr(connection, 'close'):
                connection.close()
            
            # Record duration
            results["duration"] = time.time() - start_time
            
            # Update verification results
            with self.lock:
                self.verification_results["integrity_checks"] = results
                self.verification_results["last_verification"] = datetime.now().isoformat()
                self._update_overall_status()
                self._save_results()
            
            # Record metrics
            if self.metrics:
                self.metrics.record_operation(
                    operation_type="integrity_verification",
                    success=results["status"] in ["OK", "Warning"],
                    execution_time=results["duration"],
                    items_processed=results["checks_performed"]
                )
        
        return results
def verify_consistency(self, consistency_checks: List[Dict[str, Any]]) -> Dict[str, Any]:
        """
        Verify data consistency across tables.
        
        Args:
            consistency_checks: List of consistency check dictionaries with:
                - name: Name of the check
                - description: Description of what the check verifies
                - source_query: SQL query for source data
                - target_query: SQL query for target data
                - comparison_key: Column name to use for joining results
                - comparison_columns: List of columns to compare
                - severity: 'critical', 'warning', or 'info'
                
        Returns:
            Dictionary with consistency check results
        """
        if not ConnectionManager:
            logger.error("ConnectionManager not available, cannot verify consistency")
            return {"status": "Error", "message": "ConnectionManager not available"}
        
        start_time = time.time()
        connection = None
        results = {
            "status": "Unknown",
            "timestamp": datetime.now().isoformat(),
            "duration": 0,
            "checks_performed": 0,
            "checks_failed": 0,
            "critical_failures": 0,
            "details": []
        }
        
        try:
            # Get database connection
            connection = ConnectionManager.get_instance().get_connection()
            if not connection:
                raise Exception("Could not get database connection")
            
            # Run each consistency check
            for check in consistency_checks:
                check_name = check.get('name', 'Unnamed check')
                check_description = check.get('description', '')
                source_query = check.get('source_query', '')
                target_query = check.get('target_query', '')
                comparison_key = check.get('comparison_key', '')
                comparison_columns = check.get('comparison_columns', [])
                check_severity = check.get('severity', 'warning')
                
                if not source_query or not target_query or not comparison_key:
                    logger.warning(f"Skipping consistency check '{check_name}' with missing parameters")
                    continue
                
                check_result = {
                    "name": check_name,
                    "description": check_description,
                    "severity": check_severity,
                    "status": "Unknown",
                    "inconsistencies": 0,
                    "sample_inconsistencies": []
                }
                
                try:
                    # Execute source and target queries
                    source_data = connection.execute_query(source_query)
                    target_data = connection.execute_query(target_query)
                    
                    # Convert to dictionaries keyed by comparison_key
                    source_dict = {row[comparison_key]: row for row in source_data} if source_data else {}
                    target_dict = {row[comparison_key]: row for row in target_data} if target_data else {}
                    
                    # Find inconsistencies
                    inconsistencies = []
                    
                    for key, source_row in source_dict.items():
                        if key in target_dict:
                            target_row = target_dict[key]
                            
                            # Compare specified columns
                            for column in comparison_columns:
                                if column in source_row and column in target_row:
                                    if source_row[column] != target_row[column]:
                                        inconsistencies.append({
                                            "key": key,
                                            "column": column,
                                            "source_value": source_row[column],
                                            "target_value": target_row[column]
                                        })
                        else:
                            # Key exists in source but not in target
                            inconsistencies.append({
                                "key": key,
                                "issue": "missing_in_target"
                            })
                    
                    # Check for keys in target but not in source
                    for key in target_dict:
                        if key not in source_dict:
                            inconsistencies.append({
                                "key": key,
                                "issue": "missing_in_source"
                            })
                    
                    # Store results
                    check_result["inconsistencies"] = len(inconsistencies)
                    check_result["status"] = "OK" if len(inconsistencies) == 0 else "Error"
                    
                    # Store sample inconsistencies (up to 10)
                    if inconsistencies:
                        check_result["sample_inconsistencies"] = inconsistencies[:10]
                        results["checks_failed"] += 1
                        
                        if check_severity == 'critical':
                            results["critical_failures"] += 1
                    
                except Exception as e:
                    logger.warning(f"Error running consistency check '{check_name}': {str(e)}")
                    check_result["status"] = "Error"
                    check_result["error"] = str(e)
                    results["checks_failed"] += 1
                
                results["details"].append(check_result)
                results["checks_performed"] += 1
            
            # Determine overall status
            if results["critical_failures"] > 0:
                results["status"] = "Critical"
            elif results["checks_failed"] > 0:
                results["status"] = "Warning"
            else:
                results["status"] = "OK"
            
        except Exception as e:
            logger.error(f"Error verifying consistency: {str(e)}")
            results["status"] = "Error"
            results["message"] = str(e)
        finally:
            if connection and hasattr(connection, 'close'):
                connection.close()
            
            # Record duration
            results["duration"] = time.time() - start_time
            
            # Update verification results
            with self.lock:
                self.verification_results["consistency_checks"] = results
                self.verification_results["last_verification"] = datetime.now().isoformat()
                self._update_overall_status()
                self._save_results()
            
            # Record metrics
            if self.metrics:
                self.metrics.record_operation(
                    operation_type="consistency_verification",
                    success=results["status"] in ["OK", "Warning"],
                    execution_time=results["duration"],
                    items_processed=results["checks_performed"]
                )
        
        return results
def verify_relationships(self, relationships: List[Dict[str, Any]]) -> Dict[str, Any]:
        """
        Verify foreign key relationships.
        
        Args:
            relationships: List of relationship dictionaries with:
                - name: Name of the relationship
                - source_table: Source table name
                - source_column: Source column name
                - target_table: Target table name
                - target_column: Target column name
                - severity: 'critical', 'warning', or 'info'
                
        Returns:
            Dictionary with relationship validation results
        """
        if not ConnectionManager:
            logger.error("ConnectionManager not available, cannot verify relationships")
            return {"status": "Error", "message": "ConnectionManager not available"}
        
        start_time = time.time()
        connection = None
        results = {
            "status": "Unknown",
            "timestamp": datetime.now().isoformat(),
            "duration": 0,
            "relationships_checked": 0,
            "relationships_invalid": 0,
            "critical_failures": 0,
            "details": []
        }
        
        try:
            # Get database connection
            connection = ConnectionManager.get_instance().get_connection()
            if not connection:
                raise Exception("Could not get database connection")
            
            # Check each relationship
            for rel in relationships:
                rel_name = rel.get('name', 'Unnamed relationship')
                source_table = rel.get('source_table', '')
                source_column = rel.get('source_column', '')
                target_table = rel.get('target_table', '')
                target_column = rel.get('target_column', '')
                rel_severity = rel.get('severity', 'critical')
                
                if not source_table or not source_column or not target_table or not target_column:
                    logger.warning(f"Skipping relationship check '{rel_name}' with missing parameters")
                    continue
                
                rel_result = {
                    "name": rel_name,
                    "source_table": source_table,
                    "source_column": source_column,
                    "target_table": target_table,
                    "target_column": target_column,
                    "severity": rel_severity,
                    "status": "Unknown",
                    "orphaned_records": 0,
                    "sample_orphans": []
                }
                
                try:
                    # Query to find orphaned records
                    orphan_query = f"""
                    SELECT s.{source_column} as orphaned_key
                    FROM {source_table} s
                    LEFT JOIN {target_table} t ON s.{source_column} = t.{target_column}
                    WHERE s.{source_column} IS NOT NULL
                    AND t.{target_column} IS NULL
                    LIMIT 11
                    """
                    
                    orphan_result = connection.execute_query(orphan_query)
                    orphaned_count = len(orphan_result) if orphan_result else 0
                    
                    # Store results
                    rel_result["orphaned_records"] = orphaned_count
                    rel_result["status"] = "OK" if orphaned_count == 0 else "Error"
                    
                    # Store sample orphans (up to 10)
                    if orphaned_count > 0:
                        rel_result["sample_orphans"] = [row['orphaned_key'] for row in orphan_result[:10]]
                        results["relationships_invalid"] += 1
                        
                        if rel_severity == 'critical':
                            results["critical_failures"] += 1
                    
                except Exception as e:
                    logger.warning(f"Error checking relationship '{rel_name}': {str(e)}")
                    rel_result["status"] = "Error"
                    rel_result["error"] = str(e)
                    results["relationships_invalid"] += 1
                
                results["details"].append(rel_result)
                results["relationships_checked"] += 1
            
            # Determine overall status
            if results["critical_failures"] > 0:
                results["status"] = "Critical"
            elif results["relationships_invalid"] > 0:
                results["status"] = "Warning"
            else:
                results["status"] = "OK"
            
        except Exception as e:
            logger.error(f"Error verifying relationships: {str(e)}")
            results["status"] = "Error"
            results["message"] = str(e)
        finally:
            if connection and hasattr(connection, 'close'):
                connection.close()
            
            # Record duration
            results["duration"] = time.time() - start_time
            
            # Update verification results
            with self.lock:
                self.verification_results["relationship_validation"] = results
                self.verification_results["last_verification"] = datetime.now().isoformat()
                self._update_overall_status()
                self._save_results()
            
            # Record metrics
            if self.metrics:
                self.metrics.record_operation(
                    operation_type="relationship_verification",
                    success=results["status"] in ["OK", "Warning"],
                    execution_time=results["duration"],
                    items_processed=results["relationships_checked"]
                )
        
        return results
    
   def _update_overall_status(self) -> None:
        """Update overall verification status based on individual check results."""
        status_priority = {
            "Critical": 1,
            "Error": 2,
            "Warning": 3,
            "OK": 4,
            "Unknown": 5
        }
        
        # Get status of each verification type
        statuses = []
        
        if "schema_validation" in self.verification_results and "status" in self.verification_results["schema_validation"]:
            statuses.append(self.verification_results["schema_validation"]["status"])
        
        if "count_verification" in self.verification_results and "status" in self.verification_results["count_verification"]:
            statuses.append(self.verification_results["count_verification"]["status"])
        
        if "integrity_checks" in self.verification_results and "status" in self.verification_results["integrity_checks"]:
            statuses.append(self.verification_results["integrity_checks"]["status"])
        
        if "consistency_checks" in self.verification_results and "status" in self.verification_results["consistency_checks"]:
            statuses.append(self.verification_results["consistency_checks"]["status"])
        
        if "relationship_validation" in self.verification_results and "status" in self.verification_results["relationship_validation"]:
            statuses.append(self.verification_results["relationship_validation"]["status"])
        
        # Determine overall status (highest priority status)
        if statuses:
            # Sort by priority and take the first (highest priority) status
            sorted_statuses = sorted(statuses, key=lambda s: status_priority.get(s, 999))
            self.verification_results["overall_status"] = sorted_statuses[0]
        else:
            self.verification_results["overall_status"] = "Unknown"