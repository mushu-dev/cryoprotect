#!/usr/bin/env python3
"""
PostgreSQL Health Check Tool for CryoProtect.

This script checks the health and performance of the PostgreSQL database
and provides optimization recommendations based on the results.
"""

import os
import sys
import psycopg2
import psycopg2.extras
import json
import subprocess
from datetime import datetime
from dotenv import load_dotenv

# Load environment variables
load_dotenv()

# Database connection parameters
DB_HOST = os.getenv("SUPABASE_DB_HOST", "localhost")
DB_PORT = os.getenv("SUPABASE_DB_PORT", "5432")
DB_NAME = os.getenv("SUPABASE_DB_NAME", "postgres")
DB_USER = os.getenv("SUPABASE_DB_USER", "postgres")
DB_PASSWORD = os.getenv("SUPABASE_DB_PASSWORD", "postgres")

class PostgresHealthCheck:
    """PostgreSQL health check tool."""
    
    def __init__(self):
        """Initialize the health check tool."""
        self.conn = None
        self.results = {
            "timestamp": datetime.now().isoformat(),
            "connection": {"status": "unknown"},
            "version": "unknown",
            "settings": {},
            "statistics": {},
            "issues": [],
            "recommendations": []
        }
        
    def connect(self):
        """Connect to the PostgreSQL database."""
        try:
            self.conn = psycopg2.connect(
                host=DB_HOST,
                port=DB_PORT,
                dbname=DB_NAME,
                user=DB_USER,
                password=DB_PASSWORD
            )
            self.results["connection"]["status"] = "ok"
            self.results["connection"]["details"] = {
                "host": DB_HOST,
                "port": DB_PORT,
                "dbname": DB_NAME,
                "user": DB_USER
            }
            return True
        except Exception as e:
            self.results["connection"]["status"] = "failed"
            self.results["connection"]["error"] = str(e)
            return False
            
    def check_version(self):
        """Check the PostgreSQL version."""
        try:
            cursor = self.conn.cursor()
            cursor.execute("SELECT version()")
            version = cursor.fetchone()[0]
            self.results["version"] = version
            cursor.close()
        except Exception as e:
            self.results["issues"].append({
                "severity": "error",
                "message": f"Failed to get PostgreSQL version: {e}"
            })
            
    def check_settings(self):
        """Check important PostgreSQL settings."""
        try:
            cursor = self.conn.cursor(cursor_factory=psycopg2.extras.RealDictCursor)
            
            # Get the most important settings
            settings_query = """
                SELECT name, setting, unit, context, short_desc
                FROM pg_settings
                WHERE name IN (
                    'shared_buffers', 'work_mem', 'maintenance_work_mem',
                    'effective_cache_size', 'max_connections',
                    'random_page_cost', 'effective_io_concurrency',
                    'max_worker_processes', 'max_parallel_workers',
                    'max_parallel_workers_per_gather', 'autovacuum',
                    'default_statistics_target', 'synchronous_commit',
                    'wal_buffers', 'checkpoint_completion_target',
                    'max_wal_size', 'min_wal_size'
                )
                ORDER BY name;
            """
            cursor.execute(settings_query)
            settings = cursor.fetchall()
            
            for setting in settings:
                self.results["settings"][setting["name"]] = {
                    "value": setting["setting"],
                    "unit": setting["unit"],
                    "context": setting["context"],
                    "description": setting["short_desc"]
                }
                
            cursor.close()
        except Exception as e:
            self.results["issues"].append({
                "severity": "error",
                "message": f"Failed to check PostgreSQL settings: {e}"
            })
            
    def check_statistics(self):
        """Check database statistics."""
        try:
            cursor = self.conn.cursor(cursor_factory=psycopg2.extras.RealDictCursor)
            
            # Database size
            cursor.execute("SELECT pg_size_pretty(pg_database_size(current_database())) as size")
            db_size = cursor.fetchone()
            self.results["statistics"]["database_size"] = db_size["size"]
            
            # Table counts
            cursor.execute("SELECT COUNT(*) as count FROM information_schema.tables WHERE table_schema = 'public'")
            table_count = cursor.fetchone()
            self.results["statistics"]["table_count"] = table_count["count"]
            
            # Connection statistics
            cursor.execute("""
                SELECT 
                    count(*) AS active_connections,
                    sum(CASE WHEN state = 'active' THEN 1 ELSE 0 END) AS active_queries,
                    sum(CASE WHEN state = 'idle' THEN 1 ELSE 0 END) AS idle_connections
                FROM pg_stat_activity
            """)
            connections = cursor.fetchone()
            self.results["statistics"]["connections"] = connections
            
            # Cache hit ratio
            cursor.execute("""
                SELECT 
                    sum(heap_blks_read) as heap_read,
                    sum(heap_blks_hit) as heap_hit,
                    sum(heap_blks_hit) / (sum(heap_blks_hit) + sum(heap_blks_read)) as ratio
                FROM pg_statio_user_tables
            """)
            cache = cursor.fetchone()
            self.results["statistics"]["cache_hit_ratio"] = cache["ratio"] if cache["ratio"] is not None else 0
            
            # Get largest tables
            cursor.execute("""
                SELECT
                    tablename AS table_name,
                    pg_size_pretty(pg_total_relation_size('"' || tablename || '"')) AS size
                FROM pg_tables
                WHERE schemaname = 'public'
                ORDER BY pg_total_relation_size('"' || tablename || '"') DESC
                LIMIT 5
            """)
            largest_tables = cursor.fetchall()
            self.results["statistics"]["largest_tables"] = largest_tables
            
            cursor.close()
        except Exception as e:
            self.results["issues"].append({
                "severity": "error",
                "message": f"Failed to check database statistics: {e}"
            })
            
    def analyze_performance(self):
        """Analyze database performance and suggest optimizations."""
        try:
            # Check cache hit ratio
            if "cache_hit_ratio" in self.results["statistics"]:
                ratio = self.results["statistics"]["cache_hit_ratio"]
                if ratio < 0.95:
                    self.results["issues"].append({
                        "severity": "warning",
                        "message": f"Cache hit ratio is low: {ratio:.2%}",
                        "recommendation": "Consider increasing shared_buffers"
                    })
                    self.results["recommendations"].append(
                        "Increase shared_buffers to improve cache hit ratio"
                    )
                    
            # Check connections
            if "connections" in self.results["statistics"]:
                conn_data = self.results["statistics"]["connections"]
                max_conn = int(self.results["settings"].get("max_connections", {}).get("value", 100))
                conn_percent = conn_data["active_connections"] / max_conn
                
                if conn_percent > 0.7:
                    self.results["issues"].append({
                        "severity": "warning",
                        "message": f"High connection usage: {conn_percent:.2%} of max_connections",
                        "recommendation": "Consider implementing connection pooling"
                    })
                    self.results["recommendations"].append(
                        "Implement connection pooling to reduce the number of active connections"
                    )
                    
            # Check for missing indexes on foreign keys
            cursor = self.conn.cursor(cursor_factory=psycopg2.extras.RealDictCursor)
            cursor.execute("""
                SELECT
                    conrelid::regclass AS table_from,
                    a.attname AS column,
                    confrelid::regclass AS table_to,
                    af.attname AS ftable_column,
                    (SELECT COUNT(1) FROM pg_index i WHERE indrelid=conrelid AND (a.attnum = ANY(i.indkey))) AS has_index
                FROM
                    pg_constraint c
                    JOIN pg_attribute a ON a.attnum = ANY(c.conkey) AND a.attrelid = c.conrelid
                    JOIN pg_attribute af ON af.attnum = ANY(c.confkey) AND af.attrelid = c.confrelid
                WHERE
                    c.contype = 'f'
                    AND (SELECT COUNT(1) FROM pg_index i WHERE indrelid=conrelid AND (a.attnum = ANY(i.indkey))) = 0
                ORDER BY
                    conrelid::regclass::text,
                    a.attname;
            """)
            missing_indexes = cursor.fetchall()
            
            if missing_indexes:
                self.results["issues"].append({
                    "severity": "warning",
                    "message": f"Found {len(missing_indexes)} foreign keys without indexes",
                    "details": [f"{row['table_from']}.{row['column']} -> {row['table_to']}.{row['ftable_column']}" 
                               for row in missing_indexes[:5]]
                })
                self.results["recommendations"].append(
                    "Add indexes to foreign key columns to improve join performance"
                )
            
            cursor.close()
            
        except Exception as e:
            self.results["issues"].append({
                "severity": "error",
                "message": f"Failed to analyze performance: {e}"
            })
            
    def run_checks(self):
        """Run all health checks."""
        if not self.connect():
            return
            
        self.check_version()
        self.check_settings()
        self.check_statistics()
        self.analyze_performance()
        
        # Close connection
        if self.conn:
            self.conn.close()
            
    def print_report(self):
        """Print the health check report to the console."""
        print("=" * 70)
        print("PostgreSQL Health Check Report")
        print("=" * 70)
        print(f"Timestamp: {self.results['timestamp']}")
        print(f"Connection: {self.results['connection']['status']}")
        print(f"Version: {self.results['version']}")
        print()
        
        print("Database Statistics:")
        print(f"  Database Size: {self.results['statistics'].get('database_size', 'Unknown')}")
        print(f"  Table Count: {self.results['statistics'].get('table_count', 'Unknown')}")
        if "connections" in self.results["statistics"]:
            conn = self.results["statistics"]["connections"]
            print(f"  Active Connections: {conn.get('active_connections', 'Unknown')}")
            print(f"  Active Queries: {conn.get('active_queries', 'Unknown')}")
        if "cache_hit_ratio" in self.results["statistics"]:
            print(f"  Cache Hit Ratio: {self.results['statistics']['cache_hit_ratio']:.2%}")
        print()
        
        print("Key Settings:")
        for name, details in self.results["settings"].items():
            value = details["value"]
            unit = details["unit"] if details["unit"] else ""
            print(f"  {name}: {value}{unit}")
        print()
        
        if self.results["issues"]:
            print("Issues Found:")
            for issue in self.results["issues"]:
                severity = issue["severity"].upper()
                message = issue["message"]
                print(f"  [{severity}] {message}")
                if "details" in issue:
                    for detail in issue["details"]:
                        print(f"    - {detail}")
            print()
            
        if self.results["recommendations"]:
            print("Recommendations:")
            for i, rec in enumerate(self.results["recommendations"], 1):
                print(f"  {i}. {rec}")
            print()
        
        print("=" * 70)
        
    def save_report(self, filename):
        """Save the health check report to a file."""
        with open(filename, 'w') as f:
            json.dump(self.results, f, indent=2)
        print(f"Report saved to {filename}")

def main():
    """Main function."""
    health_check = PostgresHealthCheck()
    health_check.run_checks()
    health_check.print_report()
    
    # Save report to file
    report_dir = "reports"
    os.makedirs(report_dir, exist_ok=True)
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    filename = f"{report_dir}/postgres_health_check_{timestamp}.json"
    health_check.save_report(filename)
    
    # Create symlink to latest report
    latest_link = f"{report_dir}/postgres_health_latest.json"
    if os.path.exists(latest_link):
        os.remove(latest_link)
    os.symlink(filename, latest_link)
    
    # Return exit code based on issues
    critical_issues = sum(1 for issue in health_check.results["issues"] if issue["severity"] == "error")
    return 1 if critical_issues > 0 else 0

if __name__ == "__main__":
    sys.exit(main())
