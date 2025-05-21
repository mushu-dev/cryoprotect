#!/usr/bin/env python3
"""
CryoProtect v2 - Direct PostgreSQL Connection to Supabase

This script provides direct PostgreSQL database connections to Supabase,
bypassing DNS resolution issues by using IP addresses directly.
"""

import os
import sys
import socket
import json
import psycopg2
import psycopg2.extras
from dotenv import load_dotenv

# Load environment variables
load_dotenv()

# ANSI color codes for better readability
class Colors:
    HEADER = '\033[95m'
    BLUE = '\033[94m'
    CYAN = '\033[96m'
    GREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

print(f"{Colors.HEADER}{Colors.BOLD}CryoProtect v2 - Direct PostgreSQL Connection to Supabase{Colors.ENDC}")
print(f"{Colors.CYAN}This script provides direct database connections, bypassing DNS issues.{Colors.ENDC}")
print()

# Extract Supabase configuration
SUPABASE_URL = os.getenv("SUPABASE_URL")
SUPABASE_PROJECT_ID = os.getenv("SUPABASE_PROJECT_ID")
SUPABASE_DB_HOST = os.getenv("SUPABASE_DB_HOST")
SUPABASE_DB_PORT = os.getenv("SUPABASE_DB_PORT", "5432")
SUPABASE_DB_NAME = os.getenv("SUPABASE_DB_NAME", "postgres")
SUPABASE_DB_USER = os.getenv("SUPABASE_DB_USER", "postgres")
SUPABASE_DB_PASSWORD = os.getenv("SUPABASE_DB_PASSWORD", "postgres")

print(f"{Colors.BOLD}Current Supabase Configuration:{Colors.ENDC}")
print(f"  URL: {SUPABASE_URL}")
print(f"  Project ID: {SUPABASE_PROJECT_ID}")
print(f"  DB Host: {SUPABASE_DB_HOST}")
print(f"  DB Port: {SUPABASE_DB_PORT}")
print(f"  DB Name: {SUPABASE_DB_NAME}")
print(f"  DB User: {SUPABASE_DB_USER}")
print()

# Try to find IP address of DB host
try:
    print(f"{Colors.BOLD}Attempting to resolve DB hostname...{Colors.ENDC}")
    db_ip = socket.gethostbyname(SUPABASE_DB_HOST)
    print(f"  {Colors.GREEN}DB hostname resolved: {SUPABASE_DB_HOST} -> {db_ip}{Colors.ENDC}")
    db_host_resolved = True
except socket.gaierror as e:
    print(f"  {Colors.FAIL}Failed to resolve DB hostname: {SUPABASE_DB_HOST}{Colors.ENDC}")
    print(f"  Error: {str(e)}")
    db_host_resolved = False
    
    # Try to infer DB IP from the API hostname
    print(f"\n{Colors.BOLD}Attempting to infer DB IP from API hostname...{Colors.ENDC}")
    try:
        # Extract hostname from URL
        import urllib.parse
        parsed_url = urllib.parse.urlparse(SUPABASE_URL)
        api_hostname = parsed_url.netloc
        
        # Try to resolve API hostname
        api_ip = socket.gethostbyname(api_hostname)
        print(f"  API hostname resolved: {api_hostname} -> {api_ip}")
        
        # Common pattern: db IPs are in the same subnet as API IPs
        # Typically for Supabase: 172.64.x.y -> 172.64.x.247
        ip_parts = api_ip.split('.')
        if len(ip_parts) == 4:
            if ip_parts[0] == '172' and ip_parts[1] == '64':
                db_ip = f"{ip_parts[0]}.{ip_parts[1]}.{ip_parts[2]}.247"
                print(f"  {Colors.WARNING}Inferred potential DB IP: {db_ip}{Colors.ENDC}")
            else:
                db_ip = api_ip
                print(f"  {Colors.WARNING}Using API IP as fallback for DB: {db_ip}{Colors.ENDC}")
        else:
            db_ip = "172.64.149.247"  # Hardcode based on previous findings
            print(f"  {Colors.WARNING}Using hardcoded DB IP: {db_ip}{Colors.ENDC}")
    except Exception as e:
        print(f"  {Colors.FAIL}Failed to infer DB IP: {str(e)}{Colors.ENDC}")
        db_ip = "172.64.149.247"  # Hardcode based on previous findings
        print(f"  {Colors.WARNING}Using hardcoded DB IP as last resort: {db_ip}{Colors.ENDC}")

# Ask user about Supabase dashboard
print(f"\n{Colors.BOLD}Checking project configuration...{Colors.ENDC}")
print(f"  Please verify these settings in your Supabase dashboard:")
print(f"  1. Your project is active (not paused)")
print(f"  2. The database password in your .env file is correct")
print(f"  3. Your IP address has been allowed in the Database settings")
print()
input("Press Enter to continue...")

# Try connection with either the resolved or inferred IP
print(f"\n{Colors.BOLD}Testing PostgreSQL connection...{Colors.ENDC}")
try:
    # Create connection string
    if db_host_resolved:
        conn_str = f"host={SUPABASE_DB_HOST} port={SUPABASE_DB_PORT} dbname={SUPABASE_DB_NAME} user={SUPABASE_DB_USER} password={SUPABASE_DB_PASSWORD}"
        print(f"  Using resolved hostname: {SUPABASE_DB_HOST}")
    else:
        conn_str = f"host={db_ip} port={SUPABASE_DB_PORT} dbname={SUPABASE_DB_NAME} user={SUPABASE_DB_USER} password={SUPABASE_DB_PASSWORD}"
        print(f"  Using direct IP address: {db_ip}")
    
    # Attempt connection
    print("  Connecting to database...")
    conn = psycopg2.connect(conn_str)
    print(f"  {Colors.GREEN}Connection successful!{Colors.ENDC}")
    
    # Test a simple query
    print("  Testing simple query...")
    cur = conn.cursor(cursor_factory=psycopg2.extras.RealDictCursor)
    cur.execute("SELECT current_database() as db, current_user as user")
    result = cur.fetchone()
    print(f"  {Colors.GREEN}Query successful!{Colors.ENDC}")
    print(f"  Connected to database: {result['db']} as user: {result['user']}")
    
    # Test table query
    print("  Testing molecule table query...")
    cur.execute("SELECT COUNT(*) as count FROM public.molecules")
    result = cur.fetchone()
    print(f"  {Colors.GREEN}Molecule count: {result['count']}{Colors.ENDC}")
    
    # Close connection
    cur.close()
    conn.close()
    connection_ok = True
except Exception as e:
    print(f"  {Colors.FAIL}Connection failed: {str(e)}{Colors.ENDC}")
    connection_ok = False
    
    # Provide tips based on error type
    error_str = str(e).lower()
    if "password authentication failed" in error_str:
        print(f"\n{Colors.WARNING}Invalid database password. Please check your .env file.{Colors.ENDC}")
    elif "connection refused" in error_str:
        print(f"\n{Colors.WARNING}Connection refused. The server may be down or your IP might be blocked.{Colors.ENDC}")
    elif "timeout" in error_str:
        print(f"\n{Colors.WARNING}Connection timeout. The server may be unreachable or blocked by a firewall.{Colors.ENDC}")
    elif "no route to host" in error_str:
        print(f"\n{Colors.WARNING}No route to host. Network connectivity issue or incorrect IP address.{Colors.ENDC}")

# Create helper class for database connections
if connection_ok:
    print(f"\n{Colors.BOLD}Creating database connection helper class...{Colors.ENDC}")
    helper_code = """#!/usr/bin/env python3
\"\"\"
CryoProtect v2 - Supabase PostgreSQL Connection Helper

This module provides direct PostgreSQL database connections to Supabase,
bypassing DNS resolution issues with connection pooling for better performance.
\"\"\"

import os
import threading
import logging
from typing import Dict, List, Any, Optional, Union, Tuple
import psycopg2
import psycopg2.pool
import psycopg2.extras
from dotenv import load_dotenv

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Load environment variables
load_dotenv()

class SupabasePostgresConnection:
    \"\"\"
    Singleton class for managing direct PostgreSQL connections to Supabase.
    Provides thread-safe connection pooling and SQL execution functions.
    \"\"\"
    _instance = None
    _lock = threading.RLock()
    
    def __init__(self, host=None, port=None, dbname=None, user=None, password=None):
        \"\"\"
        Initialize the connection pool. This should not be called directly.
        Use get_instance() instead.
        
        Args:
            host: Database hostname or IP (optional, defaults to env vars)
            port: Database port (optional, defaults to env vars)
            dbname: Database name (optional, defaults to env vars)
            user: Database user (optional, defaults to env vars)
            password: Database password (optional, defaults to env vars)
        \"\"\"
        if SupabasePostgresConnection._instance is not None:
            raise RuntimeError("Use SupabasePostgresConnection.get_instance() instead")
        
        # Read connection parameters from environment variables or use provided values
        self.db_host = host or os.environ.get('SUPABASE_DB_HOST', '""" + (db_ip if not db_host_resolved else SUPABASE_DB_HOST) + """')
        self.db_port = port or int(os.environ.get('SUPABASE_DB_PORT', 5432))
        self.db_name = dbname or os.environ.get('SUPABASE_DB_NAME', 'postgres')
        self.db_user = user or os.environ.get('SUPABASE_DB_USER', 'postgres')
        self.db_password = password or os.environ.get('SUPABASE_DB_PASSWORD')
        
        # Connection pool settings
        self.min_connections = int(os.environ.get('SUPABASE_DB_MIN_CONNECTIONS', 1))
        self.max_connections = int(os.environ.get('SUPABASE_DB_MAX_CONNECTIONS', 10))
        
        # Validate required parameters
        if not self.db_password:
            raise ValueError("Database password is required")
        
        # Initialize connection pool
        self._initialize_pool()
        
        # Statistics
        self.query_count = 0
        self.error_count = 0
        self.last_error = None
        
        logger.info(f"Initialized SupabasePostgresConnection with pool size {{self.min_connections}}-{{self.max_connections}}")
    
    @classmethod
    def get_instance(cls, **kwargs) -> 'SupabasePostgresConnection':
        \"\"\"
        Get the singleton instance of SupabasePostgresConnection.
        
        Args:
            **kwargs: Connection parameters to override env vars
            
        Returns:
            The singleton instance
        \"\"\"
        with cls._lock:
            if cls._instance is None:
                cls._instance = SupabasePostgresConnection(**kwargs)
            return cls._instance
    
    def _initialize_pool(self):
        \"\"\"Initialize the connection pool.\"\"\"
        try:
            self.pool = psycopg2.pool.ThreadedConnectionPool(
                minconn=self.min_connections,
                maxconn=self.max_connections,
                host=self.db_host,
                port=self.db_port,
                dbname=self.db_name,
                user=self.db_user,
                password=self.db_password
            )
            logger.info("Connection pool initialized successfully")
        except Exception as e:
            logger.error(f"Error initializing connection pool: {{str(e)}}")
            self.last_error = str(e)
            self.error_count += 1
            raise
    
    def get_connection(self):
        \"\"\"
        Get a connection from the pool.
        
        Returns:
            A database connection
            
        Raises:
            psycopg2.pool.PoolError: If unable to get a connection from the pool
        \"\"\"
        try:
            conn = self.pool.getconn()
            return conn
        except Exception as e:
            logger.error(f"Error getting connection from pool: {{str(e)}}")
            self.last_error = str(e)
            self.error_count += 1
            raise
    
    def release_connection(self, conn):
        \"\"\"
        Release a connection back to the pool.
        
        Args:
            conn: The connection to release
        \"\"\"
        try:
            self.pool.putconn(conn)
        except Exception as e:
            logger.error(f"Error releasing connection to pool: {{str(e)}}")
            # If we can't return it to the pool, try to close it
            try:
                conn.close()
            except:
                pass
    
    def close_all(self):
        \"\"\"Close all connections in the pool.\"\"\"
        try:
            self.pool.closeall()
            logger.info("All connections closed")
        except Exception as e:
            logger.error(f"Error closing connections: {{str(e)}}")
    
    def execute_query(self, query: str, params: Optional[Dict[str, Any]] = None) -> Optional[List[Dict[str, Any]]]:
        \"\"\"
        Execute a SQL query and return the results as a list of dictionaries.
        
        Args:
            query: The SQL query to execute
            params: Optional parameters for the query
            
        Returns:
            A list of dictionaries representing the query results, or None for non-SELECT queries
            
        Raises:
            Exception: If an error occurs during query execution
        \"\"\"
        conn = None
        
        try:
            conn = self.get_connection()
            with conn.cursor(cursor_factory=psycopg2.extras.RealDictCursor) as cursor:
                cursor.execute(query, params)
                
                # For SELECT queries, return the results
                if cursor.description:
                    result = cursor.fetchall()
                    # Convert from RealDictRow to regular dict
                    result = [dict(row) for row in result]
                else:
                    result = None
                
                conn.commit()
                self.query_count += 1
                return result
                
        except Exception as e:
            if conn:
                conn.rollback()
            self.error_count += 1
            self.last_error = str(e)
            logger.error(f"Error executing query: {{str(e)}}")
            raise
        finally:
            if conn:
                self.release_connection(conn)
    
    def execute_batch(self, queries: List[str], transaction: bool = True) -> bool:
        \"\"\"
        Execute a batch of SQL queries, optionally in a transaction.
        
        Args:
            queries: A list of SQL queries to execute
            transaction: Whether to execute the queries in a transaction
            
        Returns:
            True if all queries were executed successfully, False otherwise
            
        Raises:
            Exception: If an error occurs during query execution and transaction=True
        \"\"\"
        if not queries:
            return True
            
        conn = None
        
        try:
            conn = self.get_connection()
            with conn.cursor() as cursor:
                for query in queries:
                    cursor.execute(query)
                
                if transaction:
                    conn.commit()
                
                self.query_count += len(queries)
                return True
                
        except Exception as e:
            if conn and transaction:
                conn.rollback()
            self.error_count += 1
            self.last_error = str(e)
            logger.error(f"Error executing batch queries: {{str(e)}}")
            if transaction:
                raise
            return False
        finally:
            if conn:
                self.release_connection(conn)
    
    def get_stats(self) -> Dict[str, Any]:
        \"\"\"
        Get statistics about the connection pool and query execution.
        
        Returns:
            A dictionary with statistics
        \"\"\"
        return {{
            "min_connections": self.min_connections,
            "max_connections": self.max_connections,
            "query_count": self.query_count,
            "error_count": self.error_count,
            "last_error": self.last_error
        }}
    
    def test_connection(self) -> bool:
        \"\"\"
        Test the connection to the database.
        
        Returns:
            True if the connection is successful, False otherwise
        \"\"\"
        try:
            result = self.execute_query("SELECT 1 AS test")
            return result and len(result) > 0 and result[0]['test'] == 1
        except Exception:
            return False

# Example usage
if __name__ == "__main__":
    # Get instance
    try:
        db = SupabasePostgresConnection.get_instance()
        
        # Test connection
        if db.test_connection():
            print("Connection test PASSED")
            
            # Execute query
            result = db.execute_query("SELECT COUNT(*) as count FROM public.molecules")
            print(f"Molecule count: {{result[0]['count']}}")
            
            # Close connections when done
            db.close_all()
        else:
            print("Connection test FAILED")
    except Exception as e:
        print(f"Error: {{str(e)}}")
"""
    
    try:
        with open("supabase_postgres_helper.py", "w") as f:
            f.write(helper_code)
        print(f"  {Colors.GREEN}Successfully wrote helper class to supabase_postgres_helper.py{Colors.ENDC}")
    except Exception as e:
        print(f"  {Colors.FAIL}Error writing helper class: {str(e)}{Colors.ENDC}")

# Create test script
if connection_ok:
    print(f"\n{Colors.BOLD}Creating test script...{Colors.ENDC}")
    test_script = """#!/usr/bin/env python3
\"\"\"
CryoProtect v2 - Test PostgreSQL Connection

This script tests the direct PostgreSQL connection to Supabase.
\"\"\"

from supabase_postgres_helper import SupabasePostgresConnection

print("CryoProtect v2 - Test PostgreSQL Connection")
print("Testing the direct connection to Supabase PostgreSQL database")
print()

# Get instance
try:
    print("Initializing connection...")
    db = SupabasePostgresConnection.get_instance()
    
    # Test connection
    print("\\nTesting connection...")
    if db.test_connection():
        print("SUCCESS: Connection test passed!")
        
        # Test molecule query
        print("\\nTesting molecule query...")
        result = db.execute_query("SELECT COUNT(*) as count FROM public.molecules")
        print(f"SUCCESS: Found {result[0]['count']} molecules in the database")
        
        # Test molecule properties query
        print("\\nTesting property query...")
        result = db.execute_query("SELECT COUNT(*) as count FROM public.molecular_properties")
        print(f"SUCCESS: Found {result[0]['count']} molecular properties in the database")
        
        # Test getting one molecule
        print("\\nTesting single molecule retrieval...")
        result = db.execute_query("SELECT * FROM public.molecules LIMIT 1")
        if result and len(result) > 0:
            molecule = result[0]
            print(f"SUCCESS: Retrieved molecule with ID: {molecule.get('id')}")
            print(f"  Name: {molecule.get('name')}")
            print(f"  SMILES: {molecule.get('smiles', 'N/A')}")
            print(f"  Formula: {molecule.get('formula', 'N/A')}")
        else:
            print("WARNING: No molecules found in the database")
        
        # Close connections when done
        db.close_all()
        print("\\nClosed all database connections")
    else:
        print("ERROR: Connection test failed")
except Exception as e:
    print(f"ERROR: {str(e)}")

print("\\nTest complete")
"""
    
    try:
        with open("test_postgres_connection.py", "w") as f:
            f.write(test_script)
        print(f"  {Colors.GREEN}Successfully wrote test script to test_postgres_connection.py{Colors.ENDC}")
    except Exception as e:
        print(f"  {Colors.FAIL}Error writing test script: {str(e)}{Colors.ENDC}")

# Create a README for postgres connection
if connection_ok:
    print(f"\n{Colors.BOLD}Creating README...{Colors.ENDC}")
    readme = """# CryoProtect v2 - Direct PostgreSQL Connection to Supabase

This directory contains a solution for DNS resolution issues with Supabase
by providing direct PostgreSQL connections to the database.

## Issue

The issue is that the database hostname (`db.tsdlmynydfuypiugmkev.supabase.co`) cannot be resolved
through DNS. This causes connections to fail when trying to connect to the database through the
Supabase API or directly.

## Solution

The solution is to connect directly to the PostgreSQL database using the IP address instead of the
hostname. This bypasses the DNS resolution issues entirely.

## Files

- `supabase_postgres_connection.py`: Script that sets up and tests the connection
- `supabase_postgres_helper.py`: Helper class for connecting to the database
- `test_postgres_connection.py`: Script for testing the connection

## Usage

To use the direct PostgreSQL connection in your code:

```python
from supabase_postgres_helper import SupabasePostgresConnection

# Get instance
db = SupabasePostgresConnection.get_instance()

# Execute query
result = db.execute_query("SELECT * FROM public.molecules LIMIT 10")
for molecule in result:
    print(f"Molecule: {molecule['name']}")

# Close connections when done
db.close_all()
```

## Important Notes

- This is a temporary workaround until DNS resolution issues are fixed
- You need to allow your IP address in the Supabase database settings
- The database password needs to be set correctly in the .env file
- The connection information is read from the following environment variables:
  - `SUPABASE_DB_HOST`: Database hostname or IP
  - `SUPABASE_DB_PORT`: Database port (default: 5432)
  - `SUPABASE_DB_NAME`: Database name (default: postgres)
  - `SUPABASE_DB_USER`: Database user (default: postgres)
  - `SUPABASE_DB_PASSWORD`: Database password
"""
    
    try:
        with open("POSTGRES_CONNECTION_README.md", "w") as f:
            f.write(readme)
        print(f"  {Colors.GREEN}Successfully wrote README to POSTGRES_CONNECTION_README.md{Colors.ENDC}")
    except Exception as e:
        print(f"  {Colors.FAIL}Error writing README: {str(e)}{Colors.ENDC}")

# Final summary
print("\n" + "="*80)
print("SUMMARY")
print("="*80)
print(f"DB Hostname resolution: {'SUCCESS' if db_host_resolved else 'FAILED'}")
print(f"DB Connection: {'SUCCESS' if connection_ok else 'FAILED'}")

if connection_ok:
    print("\nFiles created:")
    print("- supabase_postgres_helper.py: Helper class for connecting to the database")
    print("- test_postgres_connection.py: Script for testing the connection")
    print("- POSTGRES_CONNECTION_README.md: Documentation for the solution")
    
    print("\nNext steps:")
    print("1. Run the test script: python test_postgres_connection.py")
    print("2. Update your code to use the SupabasePostgresConnection class")
    
    print(f"\n{Colors.GREEN}The PostgreSQL connection has been successfully set up!{Colors.ENDC}")
else:
    print("\nTroubleshooting steps:")
    print("1. Verify your database password in the .env file")
    print("2. Check if your project is paused in the Supabase dashboard")
    print("3. Allow your IP address in the Database settings")
    print("4. Try a different IP address for the database")
    
    print(f"\n{Colors.FAIL}The PostgreSQL connection could not be set up.{Colors.ENDC}")
    print("Please follow the troubleshooting steps and try again.")
