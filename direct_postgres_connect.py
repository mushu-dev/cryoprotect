#!/usr/bin/env python3
"""
CryoProtect v2 - Direct PostgreSQL Connection to Supabase

This is a simplified version of the PostgreSQL connection script
that avoids complex string formatting issues.
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

print("CryoProtect v2 - Direct PostgreSQL Connection to Supabase")
print("This script provides direct database connections, bypassing DNS issues.")
print()

# Extract Supabase configuration
SUPABASE_URL = os.getenv("SUPABASE_URL")
SUPABASE_PROJECT_ID = os.getenv("SUPABASE_PROJECT_ID")
SUPABASE_DB_HOST = os.getenv("SUPABASE_DB_HOST")
SUPABASE_DB_PORT = os.getenv("SUPABASE_DB_PORT", "5432")
SUPABASE_DB_NAME = os.getenv("SUPABASE_DB_NAME", "postgres")
SUPABASE_DB_USER = os.getenv("SUPABASE_DB_USER", "postgres")
SUPABASE_DB_PASSWORD = os.getenv("SUPABASE_DB_PASSWORD", "postgres")

print("Current Supabase Configuration:")
print(f"  URL: {SUPABASE_URL}")
print(f"  Project ID: {SUPABASE_PROJECT_ID}")
print(f"  DB Host: {SUPABASE_DB_HOST}")
print(f"  DB Port: {SUPABASE_DB_PORT}")
print(f"  DB Name: {SUPABASE_DB_NAME}")
print(f"  DB User: {SUPABASE_DB_USER}")
print()

# Try to find IP address of DB host
try:
    print("Attempting to resolve DB hostname...")
    db_ip = socket.gethostbyname(SUPABASE_DB_HOST)
    print(f"DB hostname resolved: {SUPABASE_DB_HOST} -> {db_ip}")
    db_host_resolved = True
except socket.gaierror as e:
    print(f"Failed to resolve DB hostname: {SUPABASE_DB_HOST}")
    print(f"Error: {str(e)}")
    db_host_resolved = False
    
    # Try to infer DB IP from the API hostname
    print("\nAttempting to infer DB IP from API hostname...")
    try:
        # Extract hostname from URL
        import urllib.parse
        parsed_url = urllib.parse.urlparse(SUPABASE_URL)
        api_hostname = parsed_url.netloc
        
        # Try to resolve API hostname
        api_ip = socket.gethostbyname(api_hostname)
        print(f"API hostname resolved: {api_hostname} -> {api_ip}")
        
        # Common pattern for Supabase: db IPs are in the same subnet as API IPs with last octet as 247
        ip_parts = api_ip.split('.')
        if len(ip_parts) == 4:
            if ip_parts[0] == '172' and ip_parts[1] == '64':
                db_ip = f"{ip_parts[0]}.{ip_parts[1]}.{ip_parts[2]}.247"
                print(f"Inferred potential DB IP: {db_ip}")
            else:
                db_ip = api_ip
                print(f"Using API IP as fallback for DB: {db_ip}")
        else:
            db_ip = "172.64.149.247"  # Hardcode based on previous findings
            print(f"Using hardcoded DB IP: {db_ip}")
    except Exception as e:
        print(f"Failed to infer DB IP: {str(e)}")
        db_ip = "172.64.149.247"  # Hardcode based on previous findings
        print(f"Using hardcoded DB IP as last resort: {db_ip}")

# Ask user about Supabase dashboard
print("\nChecking project configuration...")
print("  Please verify these settings in your Supabase dashboard:")
print("  1. Your project is active (not paused)")
print("  2. The database password in your .env file is correct")
print("  3. Your IP address has been allowed in the Database settings")
print()
input("Press Enter to continue...")

# Try connection with either the resolved or inferred IP
print("\nTesting PostgreSQL connection...")
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
    print(f"  Connection successful!")
    
    # Test a simple query
    print("  Testing simple query...")
    cur = conn.cursor(cursor_factory=psycopg2.extras.RealDictCursor)
    cur.execute("SELECT current_database() as db, current_user as user")
    result = cur.fetchone()
    print(f"  Query successful!")
    print(f"  Connected to database: {result['db']} as user: {result['user']}")
    
    # Test table query
    print("  Testing molecule table query...")
    cur.execute("SELECT COUNT(*) as count FROM public.molecules")
    result = cur.fetchone()
    print(f"  Molecule count: {result['count']}")
    
    # Close connection
    cur.close()
    conn.close()
    connection_ok = True
except Exception as e:
    print(f"  Connection failed: {str(e)}")
    connection_ok = False
    
    # Provide tips based on error type
    error_str = str(e).lower()
    if "password authentication failed" in error_str:
        print("\nInvalid database password. Please check your .env file.")
    elif "connection refused" in error_str:
        print("\nConnection refused. The server may be down or your IP might be blocked.")
    elif "timeout" in error_str:
        print("\nConnection timeout. The server may be unreachable or blocked by a firewall.")
    elif "no route to host" in error_str:
        print("\nNo route to host. Network connectivity issue or incorrect IP address.")

# Create a helper class file if connection was successful
if connection_ok:
    print("\nCreating database helper class...")
    
    # This is the DB host to use - either resolved or inferred
    effective_db_host = SUPABASE_DB_HOST if db_host_resolved else db_ip
    
    # Write the helper class to a file
    helper_file_content = """#!/usr/bin/env python3
\"\"\"
CryoProtect v2 - PostgreSQL Connection Helper

This module provides direct PostgreSQL database connections to Supabase,
bypassing DNS resolution issues by using IP addresses directly.
\"\"\"

import os
import threading
import logging
from typing import Dict, List, Any, Optional, Union
import psycopg2
import psycopg2.pool
import psycopg2.extras
from dotenv import load_dotenv

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Load environment variables
load_dotenv()

class PostgresConnection:
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
        if PostgresConnection._instance is not None:
            raise RuntimeError("Use PostgresConnection.get_instance() instead")
        
        # Read connection parameters from environment variables or use provided values
        self.db_host = host or os.environ.get('SUPABASE_DB_HOST', '""" + effective_db_host + """')
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
        
        logger.info(f"Initialized PostgresConnection with pool size {self.min_connections}-{self.max_connections}")
    
    @classmethod
    def get_instance(cls, **kwargs) -> 'PostgresConnection':
        \"\"\"
        Get the singleton instance of PostgresConnection.
        
        Args:
            **kwargs: Connection parameters to override env vars
            
        Returns:
            The singleton instance
        \"\"\"
        with cls._lock:
            if cls._instance is None:
                cls._instance = PostgresConnection(**kwargs)
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
            logger.error(f"Error initializing connection pool: {str(e)}")
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
            logger.error(f"Error getting connection from pool: {str(e)}")
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
            logger.error(f"Error releasing connection to pool: {str(e)}")
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
            logger.error(f"Error closing connections: {str(e)}")
    
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
            logger.error(f"Error executing query: {str(e)}")
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
            logger.error(f"Error executing batch queries: {str(e)}")
            if transaction:
                raise
            return False
        finally:
            if conn:
                self.release_connection(conn)
    
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
        db = PostgresConnection.get_instance()
        
        # Test connection
        if db.test_connection():
            print("Connection test PASSED")
            
            # Execute query
            result = db.execute_query("SELECT COUNT(*) as count FROM public.molecules")
            print(f"Molecule count: {result[0]['count']}")
            
            # Close connections when done
            db.close_all()
        else:
            print("Connection test FAILED")
    except Exception as e:
        print(f"Error: {str(e)}")
"""
    
    # Write to file
    with open("postgres_helper.py", "w") as f:
        f.write(helper_file_content)
    print("Database helper class written to postgres_helper.py")
    
    # Create a test script
    test_script = """#!/usr/bin/env python3
\"\"\"
CryoProtect v2 - Test PostgreSQL Connection

This script tests the direct PostgreSQL connection to Supabase.
\"\"\"

from postgres_helper import PostgresConnection

print("CryoProtect v2 - Test PostgreSQL Connection")
print("Testing the direct connection to Supabase PostgreSQL database")
print()

# Get instance
try:
    print("Initializing connection...")
    db = PostgresConnection.get_instance()
    
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
    
    # Write to file
    with open("test_postgres.py", "w") as f:
        f.write(test_script)
    print("Test script written to test_postgres.py")
    
    # Create a README file
    readme_content = """# CryoProtect v2 - Direct PostgreSQL Connection

This solution provides a direct connection to the Supabase PostgreSQL database,
bypassing DNS resolution issues.

## Files

- `direct_postgres_connect.py`: Script that tests and sets up the direct PostgreSQL connection
- `postgres_helper.py`: Helper class for connecting to the database
- `test_postgres.py`: Script to test the database connection

## Usage

To use the direct PostgreSQL connection in your code:

```python
from postgres_helper import PostgresConnection

# Get connection instance
db = PostgresConnection.get_instance()

# Execute query
result = db.execute_query("SELECT * FROM public.molecules LIMIT 10")
for molecule in result:
    print(f"Molecule: {molecule['name']}")

# Close all connections when done
db.close_all()
```

## Troubleshooting

If you have connection issues:

1. Make sure your Supabase project is active (not paused)
2. Verify your database password in the .env file
3. Add your IP address to the Supabase database allowed IP list
4. Check your network and firewall settings

## Connection Parameters

The connection uses the following environment variables:

- `SUPABASE_DB_HOST`: Database hostname or IP address
- `SUPABASE_DB_PORT`: Database port (default: 5432)
- `SUPABASE_DB_NAME`: Database name (default: postgres)
- `SUPABASE_DB_USER`: Database user (default: postgres)
- `SUPABASE_DB_PASSWORD`: Database password
"""
    
    # Write to file
    with open("POSTGRES_README.md", "w") as f:
        f.write(readme_content)
    print("README written to POSTGRES_README.md")

# Final summary
print("\n" + "="*80)
print("SUMMARY")
print("="*80)
if connection_ok:
    print("Direct PostgreSQL connection SUCCESSFUL!")
    print("\nCreated files:")
    print("- postgres_helper.py: Helper class for database connections")
    print("- test_postgres.py: Test script for the connection")
    print("- POSTGRES_README.md: Documentation")
    
    print("\nNext steps:")
    print("1. Run the test script: python test_postgres.py")
    print("2. Update your code to use the PostgresConnection class")
    print("3. Make sure to close connections with db.close_all() when done")
else:
    print("Direct PostgreSQL connection FAILED!")
    print("\nTroubleshooting tips:")
    print("1. Check if your Supabase project is active (not paused)")
    print("2. Verify your database password in the .env file")
    print("3. Add your IP address to the Supabase database allowed IP list")
    print("4. Try a different database IP address or hostname")
    print("5. Check your network and firewall settings")
    
print("\n" + "="*80)
