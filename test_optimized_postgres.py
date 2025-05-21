#!/usr/bin/env python3
"""
Test PostgreSQL Connection with Optimized Settings

This script tests the connection to the optimized PostgreSQL database
and performs basic benchmarks to verify performance improvements.
"""

import os
import sys
import time
import statistics
import psycopg2
import psycopg2.extras
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

# Default benchmark parameters
DEFAULT_ITERATIONS = 10
DEFAULT_QUERIES = [
    ("Simple SELECT", "SELECT 1"),
    ("Version Query", "SELECT version()"),
    ("Count Tables", "SELECT COUNT(*) FROM information_schema.tables WHERE table_schema = 'public'"),
    ("System Catalog", "SELECT * FROM pg_settings LIMIT 100"),
]

# Color constants
GREEN = '\033[0;32m'
YELLOW = '\033[0;33m'
RED = '\033[0;31m'
BLUE = '\033[0;34m'
NC = '\033[0m'  # No Color


def connect_to_database():
    """Connect to the PostgreSQL database."""
    try:
        print(f"Connecting to PostgreSQL at {DB_HOST}:{DB_PORT}...")
        conn = psycopg2.connect(
            host=DB_HOST,
            port=DB_PORT,
            dbname=DB_NAME,
            user=DB_USER,
            password=DB_PASSWORD
        )
        print(f"{GREEN}✓ Connected successfully!{NC}")
        return conn
    except Exception as e:
        print(f"{RED}✗ Connection failed: {str(e)}{NC}")
        return None


def run_benchmark(conn, iterations=DEFAULT_ITERATIONS, queries=DEFAULT_QUERIES, custom_queries=None):
    """Run a set of benchmark queries."""
    if not conn:
        return {}
        
    print(f"\n{BLUE}Running PostgreSQL Performance Benchmark{NC}")
    print(f"Iterations per query: {iterations}")
    
    # Combine default and custom queries
    benchmark_queries = list(queries)
    if custom_queries:
        benchmark_queries.extend(custom_queries)
    
    results = {}
    
    # Run each query multiple times
    for name, query in benchmark_queries:
        print(f"\nRunning query: {name}")
        print(f"SQL: {query}")
        
        times = []
        cursor = conn.cursor()
        
        for i in range(iterations):
            sys.stdout.write(f"\rIteration {i+1}/{iterations}...")
            sys.stdout.flush()
            
            start_time = time.time()
            try:
                cursor.execute(query)
                cursor.fetchall()  # Fetch all to complete the query
            except Exception as e:
                print(f"\n{RED}Query failed: {str(e)}{NC}")
                break
                
            end_time = time.time()
            times.append(end_time - start_time)
            
        if times:
            avg_time = statistics.mean(times)
            min_time = min(times)
            max_time = max(times)
            
            # Calculate standard deviation if we have enough samples
            if len(times) > 1:
                stdev = statistics.stdev(times)
            else:
                stdev = 0
                
            sys.stdout.write("\r" + " " * 50 + "\r")  # Clear the line
            print(f"Results for {name}:")
            print(f"  Average time: {avg_time:.6f} seconds")
            print(f"  Min time: {min_time:.6f} seconds")
            print(f"  Max time: {max_time:.6f} seconds")
            print(f"  Standard deviation: {stdev:.6f} seconds")
            
            results[name] = {
                "avg_time": avg_time,
                "min_time": min_time,
                "max_time": max_time,
                "stdev": stdev,
                "iterations": len(times)
            }
        
        cursor.close()
    
    return results


def test_molecule_queries(conn):
    """Run molecule-specific benchmark queries."""
    if not conn:
        return {}
        
    print(f"\n{BLUE}Running Molecule-Specific Queries{NC}")
    
    # Check if the molecules table exists
    cursor = conn.cursor()
    cursor.execute("SELECT EXISTS (SELECT FROM information_schema.tables WHERE table_schema = 'public' AND table_name = 'molecules')")
    table_exists = cursor.fetchone()[0]
    
    if not table_exists:
        print(f"{YELLOW}The molecules table does not exist, skipping molecule queries{NC}")
        cursor.close()
        return {}
    
    # Count molecules
    print("Counting molecules...")
    cursor.execute("SELECT COUNT(*) FROM public.molecules")
    molecule_count = cursor.fetchone()[0]
    print(f"Found {molecule_count} molecules")
    
    # Skip further tests if no molecules
    if molecule_count == 0:
        print(f"{YELLOW}No molecules found, skipping molecule queries{NC}")
        cursor.close()
        return {}
    
    custom_queries = [
        ("Count Molecules", "SELECT COUNT(*) FROM public.molecules"),
        ("Molecule Properties Count", "SELECT COUNT(*) FROM public.molecular_properties"),
    ]
    
    # Add a complex join query if we have enough data
    if molecule_count > 0:
        custom_queries.append((
            "Molecule with Properties Join",
            """
            SELECT m.id, m.name, m.smiles, mp.property_type, mp.property_name, mp.property_value 
            FROM public.molecules m
            LEFT JOIN public.molecular_properties mp ON m.id = mp.molecule_id
            LIMIT 100
            """
        ))
        
        # Add a text search query
        custom_queries.append((
            "Molecule Name Search",
            "SELECT * FROM public.molecules WHERE name ILIKE '%methyl%' LIMIT 100"
        ))
    
    cursor.close()
    return run_benchmark(conn, iterations=5, queries=custom_queries)


def check_connection_pool(conn):
    """Check if connection pooling is configured."""
    if not conn:
        return {}
        
    print(f"\n{BLUE}Checking Connection Pool Configuration{NC}")
    
    cursor = conn.cursor(cursor_factory=psycopg2.extras.RealDictCursor)
    
    # Check for connection pooler like PgBouncer
    pooler_check_query = """
    SELECT datname, usename, application_name, count(*), state
    FROM pg_stat_activity
    GROUP BY datname, usename, application_name, state
    ORDER BY count(*) DESC;
    """
    
    cursor.execute(pooler_check_query)
    activity = cursor.fetchall()
    
    # Check connection pool settings
    settings_query = """
    SELECT name, setting, unit 
    FROM pg_settings 
    WHERE name IN (
        'max_connections', 
        'superuser_reserved_connections',
        'idle_in_transaction_session_timeout'
    )
    """
    
    cursor.execute(settings_query)
    pool_settings = {row['name']: row for row in cursor.fetchall()}
    
    # Check if idle sessions exist
    cursor.execute("SELECT COUNT(*) AS count FROM pg_stat_activity WHERE state = 'idle'")
    idle_count = cursor.fetchone()['count']
    
    # Calculate connection usage
    cursor.execute("SELECT COUNT(*) AS count FROM pg_stat_activity")
    total_connections = cursor.fetchone()['count']
    
    max_connections = int(pool_settings.get('max_connections', {}).get('setting', '100'))
    connection_usage = total_connections / max_connections
    
    # Print connection pool information
    print("Connection Pool Information:")
    print(f"  Max connections: {max_connections}")
    print(f"  Current connections: {total_connections}")
    print(f"  Idle connections: {idle_count}")
    print(f"  Connection usage: {connection_usage:.2%}")
    
    # Detect pooler
    pooler_detected = False
    common_poolers = ["pgbouncer", "pooler", "pgpool"]
    
    for row in activity:
        app_name = row.get('application_name', '').lower()
        if any(pooler in app_name for pooler in common_poolers):
            pooler_detected = True
            print(f"  Detected connection pooler: {row['application_name']}")
    
    if not pooler_detected:
        print(f"{YELLOW}  No connection pooler detected{NC}")
        if connection_usage > 0.7:
            print(f"{YELLOW}  High connection usage ({connection_usage:.2%}). Consider implementing connection pooling.{NC}")
    
    cursor.close()
    
    return {
        "max_connections": max_connections,
        "current_connections": total_connections,
        "idle_connections": idle_count,
        "connection_usage": connection_usage,
        "pooler_detected": pooler_detected
    }


def format_duration(seconds):
    """Format a duration in seconds to a human-readable format."""
    if seconds < 0.001:
        return f"{seconds * 1000000:.2f} μs"
    elif seconds < 1:
        return f"{seconds * 1000:.2f} ms"
    else:
        return f"{seconds:.4f} s"


def print_summary(results, pool_info):
    """Print a summary of all benchmark results."""
    print(f"\n{BLUE}{'='*50}{NC}")
    print(f"{BLUE}PostgreSQL Performance Summary{NC}")
    print(f"{BLUE}{'='*50}{NC}")
    
    # Print connection pool summary
    if pool_info:
        print("\nConnection Pool:")
        usage_color = GREEN if pool_info["connection_usage"] < 0.7 else (YELLOW if pool_info["connection_usage"] < 0.9 else RED)
        print(f"  Usage: {usage_color}{pool_info['connection_usage']:.2%}{NC} ({pool_info['current_connections']}/{pool_info['max_connections']})")
        pooler_status = f"{GREEN}Detected{NC}" if pool_info.get("pooler_detected", False) else f"{YELLOW}Not detected{NC}"
        print(f"  Connection Pooler: {pooler_status}")
    
    # Print query performance summary
    if results:
        print("\nQuery Performance:")
        for name, data in results.items():
            avg_time = data["avg_time"]
            # Color code based on performance
            if avg_time < 0.01:
                time_color = GREEN
            elif avg_time < 0.1:
                time_color = BLUE
            elif avg_time < 1:
                time_color = YELLOW
            else:
                time_color = RED
                
            print(f"  {name}: {time_color}{format_duration(avg_time)}{NC}")
    
    print(f"\n{BLUE}Overall Assessment:{NC}")
    
    # Determine if there are any performance concerns
    concerns = []
    
    if pool_info.get("connection_usage", 0) > 0.7 and not pool_info.get("pooler_detected", False):
        concerns.append(f"{YELLOW}- High connection usage without a connection pooler{NC}")
    
    if results:
        slow_queries = [name for name, data in results.items() if data["avg_time"] > 0.5]
        if slow_queries:
            concerns.append(f"{YELLOW}- Slow queries detected: {', '.join(slow_queries)}{NC}")
    
    if concerns:
        print("\nPerformance concerns:")
        for concern in concerns:
            print(f"  {concern}")
        print("\nRecommendations:")
        if pool_info.get("connection_usage", 0) > 0.7 and not pool_info.get("pooler_detected", False):
            print("  - Implement connection pooling using PgBouncer")
        print("  - Apply optimized PostgreSQL settings for Fedora")
        print("  - Use the ./check_postgres_health.py tool for detailed analysis")
    else:
        print(f"  {GREEN}No critical performance issues detected{NC}")
    
    print(f"\n{BLUE}{'='*50}{NC}")


def main():
    """Main function."""
    # Connect to the database
    conn = connect_to_database()
    if not conn:
        return 1
    
    try:
        # Run baseline benchmarks
        benchmark_results = run_benchmark(conn)
        
        # Run molecule-specific queries
        molecule_results = test_molecule_queries(conn)
        
        # Check connection pool
        pool_info = check_connection_pool(conn)
        
        # Combine results
        all_results = {**benchmark_results, **molecule_results}
        
        # Print summary
        print_summary(all_results, pool_info)
        
        # Close connection
        conn.close()
        print("\nConnection closed")
        
        return 0
    except Exception as e:
        print(f"{RED}Error: {str(e)}{NC}")
        if conn:
            conn.close()
        return 1


if __name__ == "__main__":
    sys.exit(main())