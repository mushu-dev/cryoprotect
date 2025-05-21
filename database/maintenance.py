#!/usr/bin/env python3
"""
Database maintenance utilities for CryoProtect.

This module provides functions for database maintenance tasks like:
- Vacuum tables to reclaim space
- Reindex tables to improve performance
- Analyze tables to update statistics
- Check for long-running queries
- Clean up temporary data
- Verify data integrity
"""

import logging
import time
import sys
from typing import Dict, List, Any, Optional, Tuple
from datetime import datetime, timedelta

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Maintenance configurations
DEFAULT_VACUUM_THRESHOLD = 20  # Vacuum if > 20% of tuples are dead
DEFAULT_ANALYZE_THRESHOLD = 10  # Analyze if > 10% of tuples are modified
DEFAULT_MAX_RUNTIME = 60 * 60  # Maximum runtime in seconds (1 hour)
DEFAULT_QUERY_TIMEOUT = 60 * 10  # Timeout for long-running queries (10 minutes)

def get_table_stats() -> List[Dict[str, Any]]:
    """
    Get statistics for all tables in the database.
    
    Returns:
        List of table statistics dictionaries
    """
    try:
        from database.db import execute_query
        
        query = """
            SELECT
                schemaname,
                relname AS table_name,
                n_live_tup AS live_tuples,
                n_dead_tup AS dead_tuples,
                CASE WHEN n_live_tup > 0 
                    THEN ROUND((n_dead_tup::float / n_live_tup::float) * 100, 2)
                    ELSE 0
                END AS dead_tuple_pct,
                last_vacuum,
                last_autovacuum,
                last_analyze,
                last_autoanalyze
            FROM pg_stat_user_tables
            ORDER BY dead_tuple_pct DESC;
        """
        
        result = execute_query(query)
        logger.info(f"Retrieved statistics for {len(result)} tables")
        return result
    except Exception as e:
        logger.error(f"Error getting table statistics: {e}")
        return []

def get_index_stats() -> List[Dict[str, Any]]:
    """
    Get statistics for all indexes in the database.
    
    Returns:
        List of index statistics dictionaries
    """
    try:
        from database.db import execute_query
        
        query = """
            SELECT
                schemaname,
                relname AS table_name,
                indexrelname AS index_name,
                idx_scan AS index_scans,
                idx_tup_read AS tuples_read,
                idx_tup_fetch AS tuples_fetched,
                CASE WHEN idx_scan > 0 
                    THEN ROUND((idx_tup_fetch::float / idx_scan::float), 2)
                    ELSE 0
                END AS efficiency_ratio
            FROM pg_stat_user_indexes
            ORDER BY idx_scan DESC;
        """
        
        result = execute_query(query)
        logger.info(f"Retrieved statistics for {len(result)} indexes")
        return result
    except Exception as e:
        logger.error(f"Error getting index statistics: {e}")
        return []

def get_long_running_queries() -> List[Dict[str, Any]]:
    """
    Get information about long-running queries.
    
    Returns:
        List of query information dictionaries
    """
    try:
        from database.db import execute_query
        
        query = """
            SELECT
                pid,
                application_name,
                client_addr,
                state,
                EXTRACT(EPOCH FROM now() - query_start) AS duration_seconds,
                EXTRACT(EPOCH FROM now() - state_change) AS state_change_seconds,
                wait_event_type,
                wait_event,
                query
            FROM pg_stat_activity
            WHERE state != 'idle'
              AND query != '<IDLE>'
              AND query NOT ILIKE '%pg_stat_activity%'
              AND EXTRACT(EPOCH FROM now() - query_start) > 5
            ORDER BY duration_seconds DESC;
        """
        
        result = execute_query(query)
        logger.info(f"Found {len(result)} active queries")
        return result
    except Exception as e:
        logger.error(f"Error getting long-running queries: {e}")
        return []

def vacuum_table(table_name: str, full: bool = False) -> bool:
    """
    Vacuum a table to reclaim space and update statistics.
    
    Args:
        table_name: The table to vacuum
        full: Whether to perform a FULL vacuum (slower but more thorough)
        
    Returns:
        True if successful, False otherwise
    """
    try:
        from database.db import execute_query
        
        start_time = time.time()
        vacuum_type = "FULL" if full else ""
        query = f"VACUUM {vacuum_type} {table_name};"
        
        execute_query(query)
        
        elapsed = time.time() - start_time
        logger.info(f"Vacuumed table {table_name} in {elapsed:.2f} seconds")
        return True
    except Exception as e:
        logger.error(f"Error vacuuming table {table_name}: {e}")
        return False

def analyze_table(table_name: str) -> bool:
    """
    Analyze a table to update statistics for the query planner.
    
    Args:
        table_name: The table to analyze
        
    Returns:
        True if successful, False otherwise
    """
    try:
        from database.db import execute_query
        
        start_time = time.time()
        query = f"ANALYZE {table_name};"
        
        execute_query(query)
        
        elapsed = time.time() - start_time
        logger.info(f"Analyzed table {table_name} in {elapsed:.2f} seconds")
        return True
    except Exception as e:
        logger.error(f"Error analyzing table {table_name}: {e}")
        return False

def reindex_table(table_name: str) -> bool:
    """
    Reindex a table to rebuild all its indexes.
    
    Args:
        table_name: The table to reindex
        
    Returns:
        True if successful, False otherwise
    """
    try:
        from database.db import execute_query
        
        start_time = time.time()
        query = f"REINDEX TABLE {table_name};"
        
        execute_query(query)
        
        elapsed = time.time() - start_time
        logger.info(f"Reindexed table {table_name} in {elapsed:.2f} seconds")
        return True
    except Exception as e:
        logger.error(f"Error reindexing table {table_name}: {e}")
        return False

def cancel_long_query(pid: int) -> bool:
    """
    Cancel a long-running query by process ID.
    
    Args:
        pid: The process ID of the query to cancel
        
    Returns:
        True if successful, False otherwise
    """
    try:
        from database.db import execute_query
        
        # First try to cancel the query nicely
        query = f"SELECT pg_cancel_backend({pid});"
        result = execute_query(query)
        
        # Check if the cancellation worked
        time.sleep(2)
        check_query = f"""
            SELECT pid FROM pg_stat_activity 
            WHERE pid = {pid} AND state != 'idle';
        """
        check_result = execute_query(check_query)
        
        if not check_result:
            logger.info(f"Successfully cancelled query with PID {pid}")
            return True
        
        # If not, try to terminate the backend
        terminate_query = f"SELECT pg_terminate_backend({pid});"
        terminate_result = execute_query(terminate_query)
        
        logger.warning(f"Had to terminate backend for query with PID {pid}")
        return True
    except Exception as e:
        logger.error(f"Error cancelling query with PID {pid}: {e}")
        return False

def clean_invalidation_events(days: int = 7) -> int:
    """
    Clean up old cache invalidation events.
    
    Args:
        days: Delete events older than this many days
        
    Returns:
        Number of events deleted
    """
    try:
        from database.db import execute_query
        
        query = """
            DELETE FROM cache_invalidation_events 
            WHERE processed_at IS NOT NULL 
              AND created_at < NOW() - INTERVAL '%s days'
            RETURNING COUNT(*) as count
        """
        result = execute_query(query, (days,))
        
        if result and 'count' in result[0]:
            count = result[0]['count']
            logger.info(f"Deleted {count} old cache invalidation events")
            return count
        
        return 0
    except Exception as e:
        logger.error(f"Error cleaning up old invalidation events: {e}")
        return 0

def check_database_size() -> Dict[str, Any]:
    """
    Get the size of the database and its largest tables.
    
    Returns:
        Dictionary with database size information
    """
    try:
        from database.db import execute_query
        
        # Get database size
        db_size_query = """
            SELECT
                pg_database_size(current_database()) as size_bytes,
                pg_size_pretty(pg_database_size(current_database())) as size_pretty;
        """
        db_size = execute_query(db_size_query)[0]
        
        # Get table sizes
        table_size_query = """
            SELECT
                relname AS table_name,
                pg_total_relation_size(relid) AS total_bytes,
                pg_size_pretty(pg_total_relation_size(relid)) AS total_size,
                pg_relation_size(relid) AS table_bytes,
                pg_size_pretty(pg_relation_size(relid)) AS table_size,
                pg_total_relation_size(relid) - pg_relation_size(relid) AS index_bytes,
                pg_size_pretty(pg_total_relation_size(relid) - pg_relation_size(relid)) AS index_size
            FROM pg_catalog.pg_statio_user_tables
            ORDER BY pg_total_relation_size(relid) DESC
            LIMIT 10;
        """
        table_sizes = execute_query(table_size_query)
        
        return {
            'database_size': db_size,
            'largest_tables': table_sizes
        }
    except Exception as e:
        logger.error(f"Error checking database size: {e}")
        return {}

def check_table_bloat() -> List[Dict[str, Any]]:
    """
    Check for table bloat (tables that could benefit from vacuuming).
    
    Returns:
        List of tables with bloat information
    """
    try:
        from database.db import execute_query
        
        query = """
            WITH constants AS (
                SELECT current_setting('block_size')::numeric AS bs,
                       23 AS hdr,
                       8 AS ma
            ),
            bloat_info AS (
                SELECT
                    ma, bs, schemaname, tablename, 
                    (datawidth + (hdr + ma - (case when hdr%ma=0 then ma else hdr%ma end)))::numeric AS datahdr,
                    (maxfracsum * (bs - hdr))::numeric AS maxchunk
                FROM (
                    SELECT 
                        schemaname, tablename, hdr, ma, bs,
                        SUM((1 - null_frac) * avg_width) AS datawidth,
                        MAX(null_frac) AS maxfracsum,
                        hdr + (
                            SELECT 1 + count(*) / 8
                            FROM pg_stats s2
                            WHERE null_frac != 0 AND s2.schemaname = s.schemaname AND s2.tablename = s.tablename
                        ) AS nullhdr
                    FROM pg_stats s, constants
                    GROUP BY 1, 2, 3, 4, 5
                ) AS foo
            ),
            table_bloat AS (
                SELECT
                    schemaname, tablename, 
                    cc.relpages, bs,
                    CEIL((cc.reltuples * ((datahdr + ma - 
                        (CASE WHEN datahdr%ma=0 THEN ma ELSE datahdr%ma END)) + nullhdr + 4)) / 
                        (bs - 20)) AS otta,
                    COALESCE(c2.relname,'?') AS iname, COALESCE(c2.reltuples,0) AS ituples, COALESCE(c2.relpages,0) AS ipages,
                    COALESCE(CEIL((c2.reltuples * (datahdr - 12)) / (bs - 20)),0) AS iotta
                FROM bloat_info 
                JOIN pg_class cc ON cc.relname = bloat_info.tablename
                JOIN pg_namespace nn ON cc.relnamespace = nn.oid AND nn.nspname = bloat_info.schemaname
                LEFT JOIN pg_index i ON i.indrelid = cc.oid
                LEFT JOIN pg_class c2 ON c2.oid = i.indexrelid
            )
            SELECT
                schemaname, 
                tablename, 
                ROUND((CASE WHEN otta=0 THEN 0.0 ELSE relpages::numeric/otta END - 1) * 100, 1) AS bloat_pct,
                CASE WHEN relpages <= otta THEN 0 ELSE relpages::BIGINT - otta END AS bloat_pages,
                CASE WHEN relpages <= otta THEN 0 ELSE (bs * (relpages - otta))::BIGINT END AS bloat_bytes,
                pg_size_pretty((CASE WHEN relpages <= otta THEN 0 ELSE (bs * (relpages - otta))::BIGINT END)) AS bloat_size,
                relpages, 
                otta
            FROM table_bloat
            WHERE schemaname NOT IN ('pg_catalog', 'information_schema')
            ORDER BY bloat_pct DESC
            LIMIT 20;
        """
        
        result = execute_query(query)
        logger.info(f"Retrieved bloat information for {len(result)} tables")
        return result
    except Exception as e:
        logger.error(f"Error checking table bloat: {e}")
        return []

def auto_vacuum_analyze() -> Dict[str, List[str]]:
    """
    Automatically vacuum and analyze tables based on statistics.
    
    Returns:
        Dictionary with lists of tables that were vacuumed and analyzed
    """
    try:
        # Get table statistics
        tables = get_table_stats()
        
        # Track which tables were vacuumed and analyzed
        vacuumed_tables = []
        analyzed_tables = []
        
        for table in tables:
            table_name = table['table_name']
            schema_name = table['schemaname']
            full_table_name = f"{schema_name}.{table_name}"
            
            # Skip system tables
            if schema_name in ['pg_catalog', 'information_schema']:
                continue
                
            # Check if table needs vacuum
            if table['dead_tuple_pct'] >= DEFAULT_VACUUM_THRESHOLD:
                logger.info(f"Table {full_table_name} has {table['dead_tuple_pct']}% dead tuples, vacuuming...")
                if vacuum_table(full_table_name):
                    vacuumed_tables.append(full_table_name)
            
            # Check if table needs analyze
            last_analyze = table.get('last_analyze') or table.get('last_autoanalyze')
            if not last_analyze or (datetime.now() - last_analyze).days > 7:
                logger.info(f"Table {full_table_name} hasn't been analyzed recently, analyzing...")
                if analyze_table(full_table_name):
                    analyzed_tables.append(full_table_name)
        
        return {
            'vacuumed': vacuumed_tables,
            'analyzed': analyzed_tables
        }
    except Exception as e:
        logger.error(f"Error in auto_vacuum_analyze: {e}")
        return {
            'vacuumed': [],
            'analyzed': []
        }

def check_unused_indexes() -> List[Dict[str, Any]]:
    """
    Find unused or rarely used indexes that could be candidates for removal.
    
    Returns:
        List of rarely used indexes
    """
    try:
        from database.db import execute_query
        
        query = """
            SELECT
                schemaname,
                relname AS table_name,
                indexrelname AS index_name,
                idx_scan AS index_scans,
                pg_size_pretty(pg_relation_size(i.indexrelid)) AS index_size,
                pg_relation_size(i.indexrelid) AS index_bytes
            FROM pg_stat_user_indexes ui
            JOIN pg_index i ON ui.indexrelid = i.indexrelid
            WHERE NOT indisunique
              AND idx_scan < 100
              AND pg_relation_size(i.indexrelid) > 1024 * 1024
            ORDER BY index_scans, pg_relation_size(i.indexrelid) DESC;
        """
        
        result = execute_query(query)
        logger.info(f"Found {len(result)} unused or rarely used indexes")
        return result
    except Exception as e:
        logger.error(f"Error checking unused indexes: {e}")
        return []

def check_duplicate_indexes() -> List[Dict[str, Any]]:
    """
    Find duplicate indexes that could be candidates for removal.
    
    Returns:
        List of duplicate indexes
    """
    try:
        from database.db import execute_query
        
        query = """
            WITH index_cols AS (
                SELECT
                    i.indrelid,
                    i.indexrelid,
                    array_to_string(array_agg(a.attname ORDER BY x.ordinality), ', ') as cols,
                    array_to_string(array_agg(a.attname::text || ' ' || 
                        CASE WHEN i.indoption[x.ordinality-1] & 1 = 0 THEN 'ASC' ELSE 'DESC' END
                        ORDER BY x.ordinality), ', ') as cols_with_order,
                    array_agg(a.attname ORDER BY x.ordinality) as col_array,
                    array_agg(a.attnum ORDER BY x.ordinality) as colnum_array
                FROM
                    pg_index i
                JOIN pg_class c ON i.indrelid = c.oid
                JOIN pg_namespace n ON c.relnamespace = n.oid
                CROSS JOIN LATERAL unnest(i.indkey) WITH ORDINALITY as x(attnum, ordinality)
                JOIN pg_attribute a ON a.attrelid = i.indrelid AND a.attnum = x.attnum
                WHERE
                    n.nspname NOT IN ('pg_catalog', 'pg_toast')
                    AND a.attnum > 0
                    AND NOT i.indisprimary
                GROUP BY i.indrelid, i.indexrelid
            ),
            duplicate_candidates AS (
                SELECT
                    ic.indrelid,
                    array_agg(ic.indexrelid) as indexes,
                    ic.cols
                FROM
                    index_cols ic
                GROUP BY
                    ic.indrelid, ic.cols
                HAVING
                    count(*) > 1
            )
            SELECT
                c.relname as table_name,
                n.nspname as schema_name,
                dc.cols as columns,
                array_agg(ci.relname) as index_names,
                array_agg(ci.reltuples) as row_estimates,
                array_agg(pg_size_pretty(pg_relation_size(ci.oid))) as index_sizes,
                array_agg(s.idx_scan) as index_scans
            FROM
                duplicate_candidates dc
            JOIN pg_class c ON dc.indrelid = c.oid
            JOIN pg_namespace n ON c.relnamespace = n.oid
            JOIN pg_class ci ON ci.oid = ANY(dc.indexes)
            LEFT JOIN pg_stat_user_indexes s ON s.indexrelid = ci.oid
            GROUP BY
                c.relname, n.nspname, dc.cols
            ORDER BY
                n.nspname, c.relname, dc.cols;
        """
        
        result = execute_query(query)
        logger.info(f"Found {len(result)} sets of duplicate indexes")
        return result
    except Exception as e:
        logger.error(f"Error checking duplicate indexes: {e}")
        return []

def check_foreign_keys_without_indexes() -> List[Dict[str, Any]]:
    """
    Find foreign keys without indexes, which can lead to performance issues.
    
    Returns:
        List of foreign keys without indexes
    """
    try:
        from database.db import execute_query
        
        query = """
            SELECT c.conrelid::regclass AS table_name,
                   att.attname AS column_name,
                   c.conname AS foreign_key_name,
                   c.confrelid::regclass AS referenced_table
            FROM pg_constraint c
            JOIN pg_attribute att ON att.attrelid = c.conrelid AND att.attnum = ANY(c.conkey)
            LEFT JOIN pg_index i ON i.indrelid = c.conrelid AND 
                 (att.attnum = ANY(i.indkey) AND array_position(i.indkey, att.attnum) = 0)
            WHERE c.contype = 'f'
              AND i.indexrelid IS NULL
            ORDER BY 1, 2;
        """
        
        result = execute_query(query)
        logger.info(f"Found {len(result)} foreign keys without indexes")
        return result
    except Exception as e:
        logger.error(f"Error checking foreign keys without indexes: {e}")
        return []

def clean_old_sessions(days: int = 30) -> int:
    """
    Clean up old user sessions.
    
    Args:
        days: Delete sessions older than this many days
        
    Returns:
        Number of sessions deleted
    """
    try:
        from database.db import execute_query
        
        # Check if sessions table exists
        tables_query = """
            SELECT tablename FROM pg_tables 
            WHERE schemaname = 'public' AND tablename = 'user_sessions';
        """
        tables = execute_query(tables_query)
        
        if not tables:
            logger.info("No user_sessions table found, skipping cleanup")
            return 0
        
        # Delete old sessions
        query = """
            DELETE FROM user_sessions
            WHERE last_activity < NOW() - INTERVAL '%s days'
            RETURNING COUNT(*) as count;
        """
        
        result = execute_query(query, (days,))
        
        if result and 'count' in result[0]:
            count = result[0]['count']
            logger.info(f"Deleted {count} old user sessions")
            return count
        
        return 0
    except Exception as e:
        logger.error(f"Error cleaning up old sessions: {e}")
        return 0

def clean_old_audit_logs(days: int = 90) -> int:
    """
    Clean up old audit logs.
    
    Args:
        days: Delete logs older than this many days
        
    Returns:
        Number of logs deleted
    """
    try:
        from database.db import execute_query
        
        # Check if audit_log table exists
        tables_query = """
            SELECT tablename FROM pg_tables 
            WHERE schemaname = 'public' AND tablename = 'audit_log';
        """
        tables = execute_query(tables_query)
        
        if not tables:
            logger.info("No audit_log table found, skipping cleanup")
            return 0
        
        # Delete old logs
        query = """
            DELETE FROM audit_log
            WHERE created_at < NOW() - INTERVAL '%s days'
            RETURNING COUNT(*) as count;
        """
        
        result = execute_query(query, (days,))
        
        if result and 'count' in result[0]:
            count = result[0]['count']
            logger.info(f"Deleted {count} old audit logs")
            return count
        
        return 0
    except Exception as e:
        logger.error(f"Error cleaning up old audit logs: {e}")
        return 0

def run_maintenance(full: bool = False) -> Dict[str, Any]:
    """
    Run a comprehensive database maintenance routine.
    
    Args:
        full: Whether to perform a full maintenance (slower but more thorough)
        
    Returns:
        Dictionary with maintenance results
    """
    start_time = time.time()
    results = {
        'start_time': datetime.now().isoformat(),
        'tables_vacuumed': [],
        'tables_analyzed': [],
        'indexes_reindexed': [],
        'audit_logs_cleaned': 0,
        'sessions_cleaned': 0,
        'invalidation_events_cleaned': 0,
        'long_queries_cancelled': 0,
        'errors': []
    }
    
    try:
        logger.info(f"Starting {'full' if full else 'standard'} database maintenance...")
        
        # Check for long-running queries
        long_queries = get_long_running_queries()
        cancelled_queries = 0
        
        for query in long_queries:
            if query['duration_seconds'] > DEFAULT_QUERY_TIMEOUT:
                pid = query['pid']
                logger.warning(f"Cancelling long-running query (PID: {pid}, Duration: {query['duration_seconds']}s)")
                if cancel_long_query(pid):
                    cancelled_queries += 1
        
        results['long_queries_cancelled'] = cancelled_queries
        
        # Auto-vacuum and analyze tables based on statistics
        vacuum_results = auto_vacuum_analyze()
        results['tables_vacuumed'] = vacuum_results['vacuumed']
        results['tables_analyzed'] = vacuum_results['analyzed']
        
        # Clean up old data
        results['invalidation_events_cleaned'] = clean_invalidation_events(7)
        results['sessions_cleaned'] = clean_old_sessions(30)
        results['audit_logs_cleaned'] = clean_old_audit_logs(90)
        
        # If full maintenance, also reindex important tables
        if full:
            important_tables = ['molecules', 'molecular_properties', 'cryoprotection_scores']
            
            for table in important_tables:
                logger.info(f"Reindexing table {table}...")
                if reindex_table(table):
                    results['indexes_reindexed'].append(table)
        
        # Check for database issues
        results['unused_indexes'] = len(check_unused_indexes())
        results['duplicate_indexes'] = len(check_duplicate_indexes())
        results['foreign_keys_without_indexes'] = len(check_foreign_keys_without_indexes())
        
        # Finish up
        end_time = time.time()
        results['end_time'] = datetime.now().isoformat()
        results['duration_seconds'] = end_time - start_time
        
        logger.info(f"Maintenance completed in {results['duration_seconds']:.2f} seconds")
        return results
    except Exception as e:
        logger.error(f"Error during maintenance: {e}")
        results['errors'].append(str(e))
        results['end_time'] = datetime.now().isoformat()
        results['duration_seconds'] = time.time() - start_time
        return results

def generate_maintenance_report(results: Dict[str, Any]) -> str:
    """
    Generate a human-readable maintenance report.
    
    Args:
        results: Maintenance results dictionary
        
    Returns:
        Text report
    """
    # Generate report header
    report = [
        "Database Maintenance Report",
        "==========================",
        f"Start Time: {results.get('start_time', 'Unknown')}",
        f"End Time: {results.get('end_time', 'Unknown')}",
        f"Duration: {results.get('duration_seconds', 0):.2f} seconds",
        "",
        "Actions Performed:",
        f"- Tables Vacuumed: {len(results.get('tables_vacuumed', []))}",
        f"- Tables Analyzed: {len(results.get('tables_analyzed', []))}",
        f"- Indexes Reindexed: {len(results.get('indexes_reindexed', []))}",
        f"- Long Queries Cancelled: {results.get('long_queries_cancelled', 0)}",
        "",
        "Cleanup Results:",
        f"- Cache Invalidation Events Cleaned: {results.get('invalidation_events_cleaned', 0)}",
        f"- User Sessions Cleaned: {results.get('sessions_cleaned', 0)}",
        f"- Audit Logs Cleaned: {results.get('audit_logs_cleaned', 0)}",
        "",
        "Issues Detected:",
        f"- Unused Indexes: {results.get('unused_indexes', 0)}",
        f"- Duplicate Indexes: {results.get('duplicate_indexes', 0)}",
        f"- Foreign Keys Without Indexes: {results.get('foreign_keys_without_indexes', 0)}",
        "",
    ]
    
    # Add details for vacuumed tables
    if results.get('tables_vacuumed'):
        report.append("Tables Vacuumed:")
        for table in results.get('tables_vacuumed', []):
            report.append(f"- {table}")
        report.append("")
    
    # Add details for analyzed tables
    if results.get('tables_analyzed'):
        report.append("Tables Analyzed:")
        for table in results.get('tables_analyzed', []):
            report.append(f"- {table}")
        report.append("")
    
    # Add errors if any
    if results.get('errors'):
        report.append("Errors:")
        for error in results.get('errors', []):
            report.append(f"- {error}")
        report.append("")
    
    return "\n".join(report)

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description='Database maintenance utility')
    parser.add_argument('--full', action='store_true', help='Run full maintenance (slower but more thorough)')
    parser.add_argument('--report', action='store_true', help='Generate a detailed report after maintenance')
    parser.add_argument('--report-file', type=str, help='Save report to file')
    parser.add_argument('--verbose', action='store_true', help='Enable verbose logging')
    
    args = parser.parse_args()
    
    # Configure logging
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    
    # Run maintenance
    results = run_maintenance(full=args.full)
    
    # Generate report if requested
    if args.report or args.report_file:
        report = generate_maintenance_report(results)
        
        if args.report_file:
            try:
                with open(args.report_file, 'w') as f:
                    f.write(report)
                print(f"Maintenance report saved to {args.report_file}")
            except Exception as e:
                print(f"Error saving report to file: {e}")
                print(report)
        else:
            print(report)