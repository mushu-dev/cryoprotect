"""
Performance Analyzer Module for CryoProtect v2

This module provides a class-based performance analysis utility for the database health check system.
It analyzes database performance including query execution times, index usage, slow queries,
connection pool utilization, and table statistics.
"""

import logging
import time
import json
import traceback
from typing import Dict, List, Any, Optional, Tuple
from datetime import datetime, timedelta

logger = logging.getLogger(__name__)

class PerformanceAnalyzer:
    """
    Database performance analyzer for health checks.
    
    This class provides methods to analyze database performance:
    - Measure query execution times for common operations
    - Check index usage and effectiveness
    - Identify slow queries and bottlenecks
    - Monitor connection pool utilization
    - Analyze table statistics and growth patterns
    """
    
    # Class attribute for categorization in health check system
    category = "performance"
    
    def __init__(self):
        """Initialize the performance analyzer."""
        # Define common queries to test performance
        self.test_queries = [
            {
                'name': 'get_molecules_with_properties',
                'query': '''
                    SELECT m.id, m.name, m.smiles, mp.property_name, mp.numeric_value
                    FROM molecules m
                    JOIN molecular_properties mp ON m.id = mp.molecule_id
                    WHERE mp.property_name = 'Molecular Weight'
                    LIMIT 10;
                ''',
                'threshold_ms': 100
            },
            {
                'name': 'get_mixtures_with_components',
                'query': '''
                    SELECT m.id, m.name, mc.molecule_id, mc.concentration
                    FROM mixtures m
                    JOIN mixture_components mc ON m.id = mc.mixture_id
                    LIMIT 10;
                ''',
                'threshold_ms': 100
            },
            {
                'name': 'get_experiments_with_verifications',
                'query': '''
                    SELECT e.id, e.mixture_id, e.numeric_value, lv.verification_status
                    FROM experiments e
                    LEFT JOIN lab_verifications lv ON e.id = lv.experiment_id
                    LIMIT 10;
                ''',
                'threshold_ms': 100
            }
        ]
        
        # Define important tables to analyze
        self.important_tables = [
            'molecules', 'mixtures', 'mixture_components', 'molecular_properties',
            'experiments', 'lab_verifications', 'predictions'
        ]
        
        # Define expected indexes for optimal performance
        self.expected_indexes = {
            'molecules': ['id', 'name'],
            'mixtures': ['id', 'created_by'],
            'mixture_components': ['mixture_id', 'molecule_id'],
            'molecular_properties': ['molecule_id', 'property_name'],
            'experiments': ['mixture_id', 'property_type_id'],
            'lab_verifications': ['experiment_id']
        }
    
    def run_checks(self, connection_pool) -> Dict[str, Any]:
        """
        Run all performance checks.
        
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
            query_time_results = self._check_query_execution_times(conn)
            index_usage_results = self._check_index_usage(conn)
            slow_query_results = self._check_slow_queries(conn)
            connection_pool_results = self._check_connection_pool(conn)
            table_stats_results = self._check_table_statistics(conn)
            
            # Combine results
            all_checks = [
                ('query_execution_times', query_time_results),
                ('index_usage', index_usage_results),
                ('slow_queries', slow_query_results),
                ('connection_pool', connection_pool_results),
                ('table_statistics', table_stats_results)
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
            logger.error(f"Error during performance checks: {str(e)}")
            results['status'] = 'error'
            results['details']['error'] = str(e)
            results['details']['traceback'] = traceback.format_exc()
            results['recommendations'].append("Fix the error that prevented performance checks from completing")
        
        return results
    
    def _check_query_execution_times(self, conn) -> Dict[str, Any]:
        """Measure execution times for common database operations."""
        results = {
            'status': 'passed',
            'issues_found': 0,
            'query_times': {},
            'recommendations': []
        }
        
        try:
            for query_info in self.test_queries:
                query_name = query_info['name']
                query = query_info['query']
                threshold_ms = query_info['threshold_ms']
                
                # Measure execution time
                start_time = time.time()
                query_response = conn.rpc('exec_sql', {'query': query}).execute()
                execution_time_ms = (time.time() - start_time) * 1000
                
                # Store the result
                results['query_times'][query_name] = {
                    'execution_time_ms': execution_time_ms,
                    'threshold_ms': threshold_ms,
                    'status': 'ok' if execution_time_ms <= threshold_ms else 'slow'
                }
                
                # Check if execution time exceeds threshold
                if execution_time_ms > threshold_ms:
                    results['status'] = 'warning'
                    results['issues_found'] += 1
                    results['recommendations'].append(
                        f"Optimize query '{query_name}' which took {execution_time_ms:.2f}ms (threshold: {threshold_ms}ms)"
                    )
                    
                    # Get query plan for slow queries
                    try:
                        plan_response = conn.rpc('exec_sql', {'query': f'EXPLAIN ANALYZE {query}'}).execute()
                        if hasattr(plan_response, 'data') and plan_response.data:
                            results['query_times'][query_name]['query_plan'] = plan_response.data
                    except Exception as e:
                        logger.warning(f"Could not get query plan for {query_name}: {str(e)}")
        
        except Exception as e:
            logger.error(f"Error checking query execution times: {str(e)}")
            results['status'] = 'error'
            results['issues_found'] += 1
            results['error'] = str(e)
        
        return results
    
    def _check_index_usage(self, conn) -> Dict[str, Any]:
        """Check index usage and effectiveness."""
        results = {
            'status': 'passed',
            'issues_found': 0,
            'missing_indexes': [],
            'recommendations': []
        }
        
        try:
            # Check for missing indexes on foreign keys
            missing_indexes_query = '''
                SELECT
                    t.relname AS table_name,
                    a.attname AS column_name
                FROM pg_constraint c
                JOIN pg_class t ON t.oid = c.conrelid
                JOIN pg_attribute a ON a.attrelid = t.oid AND a.attnum = ANY(c.conkey)
                WHERE c.contype = 'f'
                AND t.relnamespace = (SELECT oid FROM pg_namespace WHERE nspname = 'public')
                AND NOT EXISTS (
                    SELECT 1
                    FROM pg_index i
                    JOIN pg_attribute a2 ON a2.attrelid = t.oid AND a2.attnum = ANY(i.indkey)
                    WHERE i.indrelid = t.oid AND a2.attname = a.attname
                );
            '''
            
            missing_response = conn.rpc('exec_sql', {'query': missing_indexes_query}).execute()
            
            if hasattr(missing_response, 'data') and missing_response.data:
                missing_indexes = missing_response.data
                if missing_indexes:
                    results['missing_indexes'] = missing_indexes
                    results['status'] = 'warning'
                    results['issues_found'] += len(missing_indexes)
                    
                    for idx in missing_indexes:
                        results['recommendations'].append(
                            f"Create index on {idx['table_name']}({idx['column_name']}) to improve join performance"
                        )
        
        except Exception as e:
            logger.error(f"Error checking index usage: {str(e)}")
            results['status'] = 'error'
            results['issues_found'] += 1
            results['error'] = str(e)
        
        return results
    
    def _check_slow_queries(self, conn) -> Dict[str, Any]:
        """Identify slow queries and bottlenecks."""
        results = {
            'status': 'passed',
            'issues_found': 0,
            'slow_queries': [],
            'recommendations': []
        }
        
        try:
            # Check for long-running active queries
            active_queries_query = '''
                SELECT
                    pid,
                    now() - pg_stat_activity.query_start AS duration,
                    substring(query, 1, 100) AS query_snippet
                FROM pg_stat_activity
                WHERE state = 'active'
                AND now() - pg_stat_activity.query_start > interval '5 seconds'
                ORDER BY duration DESC;
            '''
            
            active_response = conn.rpc('exec_sql', {'query': active_queries_query}).execute()
            
            if hasattr(active_response, 'data') and active_response.data:
                active_queries = active_response.data
                if active_queries:
                    results['slow_queries'] = active_queries
                    results['status'] = 'warning'
                    results['issues_found'] += len(active_queries)
                    
                    for query in active_queries:
                        results['recommendations'].append(
                            f"Long-running query (PID {query['pid']}, duration: {query['duration']}): "
                            f"'{query['query_snippet']}...'"
                        )
        
        except Exception as e:
            logger.error(f"Error checking slow queries: {str(e)}")
            results['status'] = 'error'
            results['issues_found'] += 1
            results['error'] = str(e)
        
        return results
    
    def _check_connection_pool(self, conn) -> Dict[str, Any]:
        """Monitor connection pool utilization."""
        results = {
            'status': 'passed',
            'issues_found': 0,
            'pool_stats': {},
            'recommendations': []
        }
        
        try:
            # Check active connections
            connections_query = '''
                SELECT
                    count(*) AS total_connections,
                    count(*) FILTER (WHERE state = 'active') AS active_connections,
                    count(*) FILTER (WHERE state = 'idle') AS idle_connections,
                    count(*) FILTER (WHERE state = 'idle in transaction') AS idle_in_transaction
                FROM pg_stat_activity
                WHERE backend_type = 'client backend';
            '''
            
            connections_response = conn.rpc('exec_sql', {'query': connections_query}).execute()
            
            if hasattr(connections_response, 'data') and connections_response.data:
                conn_stats = connections_response.data[0]
                results['pool_stats']['connections'] = conn_stats
                
                # Check for connection pool issues
                if conn_stats.get('idle_in_transaction', 0) > 5:
                    results['status'] = 'warning'
                    results['issues_found'] += 1
                    results['recommendations'].append(
                        f"High number of idle-in-transaction connections ({conn_stats['idle_in_transaction']}). "
                        "Check for transactions not being properly committed or rolled back."
                    )
        
        except Exception as e:
            logger.error(f"Error checking connection pool: {str(e)}")
            results['status'] = 'error'
            results['issues_found'] += 1
            results['error'] = str(e)
        
        return results
    
    def _check_table_statistics(self, conn) -> Dict[str, Any]:
        """Analyze table statistics and growth patterns."""
        results = {
            'status': 'passed',
            'issues_found': 0,
            'table_stats': {},
            'recommendations': []
        }
        
        try:
            # Get table sizes and row counts
            tables_query = '''
                SELECT
                    relname AS table_name,
                    pg_size_pretty(pg_total_relation_size(c.oid)) AS total_size,
                    pg_size_pretty(pg_relation_size(c.oid)) AS table_size,
                    pg_size_pretty(pg_total_relation_size(c.oid) - pg_relation_size(c.oid)) AS index_size,
                    reltuples::bigint AS row_estimate
                FROM pg_class c
                JOIN pg_namespace n ON n.oid = c.relnamespace
                WHERE relkind = 'r'
                AND n.nspname = 'public'
                AND c.relname IN (
                    'molecules', 'mixtures', 'mixture_components', 'molecular_properties',
                    'experiments', 'lab_verifications', 'predictions'
                )
                ORDER BY pg_total_relation_size(c.oid) DESC;
            '''
            
            tables_response = conn.rpc('exec_sql', {'query': tables_query}).execute()
            
            if hasattr(tables_response, 'data') and tables_response.data:
                table_stats = tables_response.data
                results['table_stats']['sizes'] = table_stats
                
                # Check for tables that might need optimization
                for table in table_stats:
                    table_name = table['table_name']
                    row_estimate = table['row_estimate']
                    
                    # Check if the table has many rows but no recent analyze
                    if row_estimate > 10000:
                        # Check last analyze time
                        analyze_query = f'''
                            SELECT
                                last_analyze,
                                last_autoanalyze
                            FROM pg_stat_user_tables
                            WHERE schemaname = 'public'
                            AND relname = '{table_name}';
                        '''
                        
                        analyze_response = conn.rpc('exec_sql', {'query': analyze_query}).execute()
                        
                        if hasattr(analyze_response, 'data') and analyze_response.data:
                            analyze_info = analyze_response.data[0]
                            last_analyze = analyze_info.get('last_analyze')
                            last_autoanalyze = analyze_info.get('last_autoanalyze')
                            
                            # If neither has happened recently or at all
                            if (not last_analyze and not last_autoanalyze) or \
                               (last_analyze is None and last_autoanalyze is None):
                                results['status'] = 'warning'
                                results['issues_found'] += 1
                                results['recommendations'].append(
                                    f"Table '{table_name}' has {row_estimate} rows but has never been analyzed. "
                                    f"Run ANALYZE {table_name};"
                                )
        
        except Exception as e:
            logger.error(f"Error checking table statistics: {str(e)}")
            results['status'] = 'error'
            results['issues_found'] += 1
            results['error'] = str(e)
        
        return results
