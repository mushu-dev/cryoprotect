#!/usr/bin/env python3
"""
Debug the entire connection flow from config loading to query execution.

This script provides end-to-end tracing of the database connection process,
helping identify where issues might be occurring in the connection chain.

Usage:
  python debug_connection_flow.py [--config-path PATH] [--verbose]
"""

import os
import sys
import json
import argparse
import traceback
from pathlib import Path
from datetime import datetime

try:
    import psycopg2
    from psycopg2.extras import DictCursor, RealDictCursor
except ImportError:
    print("Error: psycopg2 package not installed. Please install it with:")
    print("  pip install psycopg2-binary")
    sys.exit(1)

# Set up argument parsing
parser = argparse.ArgumentParser(description="Debug the database connection flow")
parser.add_argument('--config-path', help='Path to the configuration file')
parser.add_argument('--verbose', action='store_true', help='Enable verbose output')
args = parser.parse_args()

# Configure verbose output
VERBOSE = args.verbose

def log(message, level=0, always=False):
    """Log a message with appropriate indentation."""
    if always or VERBOSE:
        indent = "  " * level
        print(f"{indent}{message}")

def log_section(title):
    """Log a section header."""
    log("\n" + "=" * 50, always=True)
    log(f" {title} ", always=True)
    log("=" * 50, always=True)

def find_config_path():
    """Find the configuration path from various sources."""
    log_section("CONFIG PATH RESOLUTION")
    
    # Check command line argument
    if args.config_path:
        log(f"Using config path from command line: {args.config_path}", always=True)
        return args.config_path
    
    # Check environment variable
    env_var = os.environ.get('CRYOPROTECT_CONFIG_PATH')
    if env_var:
        log(f"Using config path from CRYOPROTECT_CONFIG_PATH: {env_var}", always=True)
        return env_var
    
    # Check default locations
    default_paths = [
        Path.home() / '.cryoprotect' / 'config.json',
        Path.cwd() / 'config.json',
        Path.cwd() / 'config' / 'config.json',
        Path(__file__).parent / 'config.json',
        Path(__file__).parent / 'config' / 'config.json',
    ]
    
    for path in default_paths:
        if path.exists():
            log(f"Using config from default location: {path}", always=True)
            return str(path)
    
    log("Could not find configuration file!", always=True)
    return None

def load_config(config_path):
    """Load the configuration file."""
    log_section("CONFIG LOADING")
    
    if not config_path:
        log("No configuration path available", always=True)
        return None
    
    try:
        log(f"Loading configuration from: {config_path}")
        with open(config_path, 'r') as f:
            config = json.load(f)
            
        # Log config structure (excluding sensitive info)
        log("Configuration structure:", level=1)
        safe_config = sanitize_config(config)
        log(json.dumps(safe_config, indent=2), level=2)
        
        return config
    except FileNotFoundError:
        log(f"Configuration file not found: {config_path}", always=True)
    except json.JSONDecodeError as e:
        log(f"Invalid JSON in configuration file: {e}", always=True)
    except Exception as e:
        log(f"Error loading configuration: {e}", always=True)
        if VERBOSE:
            traceback.print_exc()
    
    return None

def sanitize_config(config):
    """Create a copy of the config with passwords masked."""
    if not isinstance(config, dict):
        return config
    
    result = {}
    for key, value in config.items():
        if key.lower() in ('password', 'secret', 'key', 'token'):
            result[key] = '********'
        elif isinstance(value, dict):
            result[key] = sanitize_config(value)
        elif isinstance(value, list):
            result[key] = [sanitize_config(item) if isinstance(item, dict) else item 
                          for item in value]
        else:
            result[key] = value
    
    return result

def extract_connection_params(config):
    """Extract database connection parameters from the configuration."""
    log_section("CONNECTION PARAMETERS")
    
    if not config:
        log("No configuration available", always=True)
        return None
    
    try:
        # Try to extract database connection parameters
        log("Extracting connection parameters from config")
        
        # Check for database config
        if 'database' not in config:
            log("Error: 'database' section not found in config", always=True)
            return None
        
        db_config = config['database']
        
        # Check for connection config
        if 'connection' not in db_config:
            log("Error: 'connection' section not found in database config", always=True)
            return None
        
        conn_config = db_config['connection']
        
        # Check connection mode
        mode = conn_config.get('mode')
        log(f"Connection mode: {mode}", level=1)
        
        # Extract parameters based on mode
        params = {}
        
        if mode == 'supabase':
            log("Processing Supabase connection parameters", level=1)
            supabase_config = conn_config.get('supabase', {})
            
            # Required parameters
            required_params = ['host', 'port', 'database', 'user', 'password', 'project_id']
            for param in required_params:
                if param not in supabase_config:
                    log(f"Error: Missing required parameter '{param}' for Supabase connection", always=True)
                    return None
                
                params[param] = supabase_config[param]
                if param != 'password':
                    log(f"{param}: {params[param]}", level=2)
                else:
                    log(f"{param}: ********", level=2)
            
            # Default values
            params['application_name'] = supabase_config.get('application_name', 'CryoProtect')
            log(f"application_name: {params['application_name']}", level=2)
            
        elif mode == 'postgres':
            log("Processing PostgreSQL connection parameters", level=1)
            postgres_config = conn_config.get('postgres', {})
            
            # Required parameters
            required_params = ['host', 'port', 'database', 'user', 'password']
            for param in required_params:
                if param not in postgres_config:
                    log(f"Error: Missing required parameter '{param}' for PostgreSQL connection", always=True)
                    return None
                
                params[param] = postgres_config[param]
                if param != 'password':
                    log(f"{param}: {params[param]}", level=2)
                else:
                    log(f"{param}: ********", level=2)
            
            # Default values
            params['application_name'] = postgres_config.get('application_name', 'CryoProtect')
            log(f"application_name: {params['application_name']}", level=2)
            
        else:
            log(f"Error: Unsupported connection mode: {mode}", always=True)
            return None
        
        # Add SSL parameters if present
        ssl_config = conn_config.get('ssl', {})
        if ssl_config:
            log("SSL configuration found", level=1)
            params['sslmode'] = ssl_config.get('mode', 'prefer')
            log(f"sslmode: {params['sslmode']}", level=2)
            
            # Add SSL certificate paths if provided
            for ssl_param in ['sslrootcert', 'sslcert', 'sslkey']:
                if ssl_param in ssl_config:
                    params[ssl_param] = ssl_config[ssl_param]
                    log(f"{ssl_param}: {params[ssl_param]}", level=2)
        
        return params
    
    except Exception as e:
        log(f"Error extracting connection parameters: {e}", always=True)
        if VERBOSE:
            traceback.print_exc()
    
    return None

def test_connection(params):
    """Test direct connection to the database."""
    log_section("DIRECT CONNECTION TEST")
    
    if not params:
        log("No connection parameters available", always=True)
        return False
    
    try:
        log(f"Connecting to {params['host']}:{params['port']} as {params['user']}...")
        
        # Start timing
        start_time = datetime.now()
        
        # Connect directly using psycopg2
        connection = psycopg2.connect(
            host=params['host'],
            port=params['port'],
            database=params['database'],
            user=params['user'],
            password=params['password'],
            application_name=params['application_name']
        )
        
        # Calculate connection time
        connection_time = (datetime.now() - start_time).total_seconds()
        log(f"✓ Connection successful! ({connection_time:.2f} seconds)", always=True)
        
        # Get connection details
        log("Connection details:", level=1)
        log(f"Server version: {connection.server_version}", level=2)
        log(f"Protocol version: {connection.protocol_version}", level=2)
        
        # Test a simple query
        log("Testing simple query...", level=1)
        with connection.cursor(cursor_factory=RealDictCursor) as cursor:
            cursor.execute("SELECT current_database(), current_user, version();")
            db, user, version = cursor.fetchone()
            log(f"Database: {db}", level=2)
            log(f"User: {user}", level=2)
            log(f"Version: {version}", level=2)
        
        # Test transaction support
        log("Testing transaction support...", level=1)
        try:
            # Set autocommit mode first before any transactions
            connection.autocommit = True
            
            # Now test a transaction
            connection.autocommit = False
            with connection.cursor() as cursor:
                cursor.execute("SELECT 1;")
            connection.commit()
            log("✓ Transaction support verified", level=2)
        except Exception as e:
            log(f"✗ Transaction error: {e}", level=2, always=True)
        
        # Close the connection
        connection.close()
        log("Connection closed", level=1)
        
        return True
    
    except psycopg2.OperationalError as e:
        log(f"✗ Connection failed: {e}", always=True)
        log("This usually indicates network issues, incorrect credentials, or server unavailability", level=1, always=True)
        
        error_str = str(e)
        if "password authentication failed" in error_str:
            log("The password appears to be incorrect", level=2, always=True)
        elif "could not connect to server" in error_str:
            log("Could not reach the database server - check network and firewall settings", level=2, always=True)
        elif "database" in error_str and "does not exist" in error_str:
            log("The specified database doesn't exist", level=2, always=True)
        elif "role" in error_str and "does not exist" in error_str:
            log("The specified user doesn't exist", level=2, always=True)
    
    except Exception as e:
        log(f"✗ Error connecting to database: {e}", always=True)
        if VERBOSE:
            traceback.print_exc()
    
    return False

def test_adapter_connection(config, params):
    """Test connection using the adapter system."""
    log_section("ADAPTER CONNECTION TEST")
    
    if not config or not params:
        log("Missing configuration or parameters", always=True)
        return False
    
    try:
        # Try to import the adapter system
        log("Importing adapter system")
        
        try:
            from database.adapter_factory import ConnectionFactory
            log("✓ Adapter factory module imported successfully", level=1)
        except ImportError as e:
            log(f"✗ Could not import adapter factory: {e}", level=1, always=True)
            log("The adapter factory module may not be in the Python path", level=2)
            
            # Try fallback to old adapter module
            try:
                from database.adapter import get_adapter
                log("✓ Falling back to legacy adapter module", level=1)
                
                # Get an adapter
                adapter = get_adapter(config)
                
                # Test a simple query
                try:
                    cursor = adapter.execute_query("SELECT current_database(), current_user;")
                    log(f"✓ Legacy adapter connection successful", level=1, always=True)
                    log(f"Connection test successful", level=2)
                    return True
                except Exception as e:
                    log(f"✗ Failed to execute query through legacy adapter: {e}", level=1, always=True)
                    return False
            except ImportError:
                log("✗ Legacy adapter module not found either", level=1, always=True)
            except Exception as e:
                log(f"✗ Error with legacy adapter: {e}", level=1, always=True)
                if VERBOSE:
                    traceback.print_exc()
            
            return False
        
        # Create a connection factory
        log("Creating connection factory")
        try:
            factory = ConnectionFactory(config)
            log("✓ Connection factory created successfully", level=1)
        except Exception as e:
            log(f"✗ Failed to create connection factory: {e}", level=1, always=True)
            if VERBOSE:
                traceback.print_exc()
            return False
        
        # Get a connection
        log("Getting connection from factory")
        try:
            conn = factory.get_connection()
            log("✓ Got connection from factory", level=1)
        except Exception as e:
            log(f"✗ Failed to get connection from factory: {e}", level=1, always=True)
            if VERBOSE:
                traceback.print_exc()
            return False
        
        # Test the connection
        log("Testing adapter connection")
        try:
            with conn.cursor() as cursor:
                cursor.execute("SELECT current_database(), current_user;")
                db, user = cursor.fetchone()
                log(f"✓ Adapter connection successful", level=1, always=True)
                log(f"Database: {db}", level=2)
                log(f"User: {user}", level=2)
            
            # Return connection to pool instead of closing directly
            from database import db
            db.release_connection(conn)
            log("Connection returned to pool", level=1)
            
            return True
        except Exception as e:
            log(f"✗ Failed to execute query through adapter: {e}", level=1, always=True)
            if VERBOSE:
                traceback.print_exc()
    
    except Exception as e:
        log(f"✗ Error in adapter test: {e}", always=True)
        if VERBOSE:
            traceback.print_exc()
    
    return False

def test_database_modules(config):
    """Test the database modules."""
    log_section("DATABASE MODULES TEST")
    
    if not config:
        log("No configuration available", always=True)
        return False
    
    all_successful = True
    
    # Test the main database module
    log("Testing database.db module")
    try:
        from database import db
        
        # First make sure the connection pool is initialized
        if not hasattr(db, '_pool') or db._pool is None:
            log("Initializing db connection pool", level=1)
            pool_success = db.init_connection_pool(config=config)
            if not pool_success:
                log("✗ Failed to initialize db connection pool", level=1, always=True)
                all_successful = False
        
        conn = db.get_connection()
        if conn is None:
            log("✗ db.get_connection() returned None", level=1, always=True)
            all_successful = False
        else:
            with conn.cursor() as cursor:
                cursor.execute("SELECT current_database(), current_user;")
                db_name, user = cursor.fetchone()
                log(f"✓ database.db connection successful", level=1, always=True)
                log(f"Database: {db_name}", level=2)
                log(f"User: {user}", level=2)
            # Release connection back to pool
            db.release_connection(conn)
    except ImportError as e:
        log(f"✗ Could not import database.db: {e}", level=1, always=True)
        all_successful = False
    except Exception as e:
        log(f"✗ Error using database.db module: {e}", level=1, always=True)
        if VERBOSE:
            traceback.print_exc()
        all_successful = False
    
    # Test the service role module if it exists
    log("Testing database.db_service_role module")
    try:
        from database import db_service_role
        
        # Extract Supabase-specific connection parameters for service role
        service_role_config = None
        if 'database' in config and 'connection' in config['database']:
            connection_config = config['database']['connection']
            if 'supabase' in connection_config:
                service_role_config = connection_config['supabase'].copy()
                # Add pooling settings if available
                if 'pooling' in config['database']:
                    pooling = config['database']['pooling']
                    if 'min_connections' in pooling:
                        service_role_config['min_connections'] = pooling['min_connections']
                    if 'max_connections' in pooling:
                        service_role_config['max_connections'] = pooling['max_connections']
                # Make sure options for service role are set
                service_role_config['options'] = "-c role=service_role -c statement_timeout=60000"
        
        # First make sure the connection pool is initialized
        if not hasattr(db_service_role, '_pool') or db_service_role._pool is None:
            log("Initializing db_service_role connection pool", level=1)
            pool_success = db_service_role.init_connection_pool(config=service_role_config)
            if not pool_success:
                log("✗ Failed to initialize db_service_role connection pool", level=1, always=True)
                all_successful = False
        
        conn = db_service_role.get_connection()
        if conn is None:
            log("✗ db_service_role.get_connection() returned None", level=1, always=True)
            all_successful = False
        else:
            with conn.cursor() as cursor:
                cursor.execute("SELECT current_database(), current_user;")
                db_name, user = cursor.fetchone()
                log(f"✓ database.db_service_role connection successful", level=1, always=True)
                log(f"Database: {db_name}", level=2)
                log(f"User: {user}", level=2)
            # Release connection back to pool
            db_service_role.release_connection(conn)
    except ImportError as e:
        log(f"✗ Could not import database.db_service_role: {e}", level=1, always=True)
        log("This module may not exist in your codebase", level=2)
    except Exception as e:
        log(f"✗ Error using database.db_service_role module: {e}", level=1, always=True)
        if VERBOSE:
            traceback.print_exc()
        all_successful = False
    
    # Test the public module if it exists
    log("Testing database.db_public module")
    try:
        from database import db_public
        
        # Extract Supabase-specific connection parameters for public access
        public_config = None
        if 'database' in config and 'connection' in config['database']:
            connection_config = config['database']['connection']
            if 'supabase' in connection_config:
                public_config = connection_config['supabase'].copy()
                # Add pooling settings if available
                if 'pooling' in config['database']:
                    pooling = config['database']['pooling']
                    if 'min_connections' in pooling:
                        public_config['min_connections'] = pooling['min_connections']
                    if 'max_connections' in pooling:
                        public_config['max_connections'] = pooling['max_connections']
        
        # First make sure the connection pool is initialized
        if not hasattr(db_public, '_pool') or db_public._pool is None:
            log("Initializing db_public connection pool", level=1)
            pool_success = db_public.init_connection_pool(config=public_config)
            if not pool_success:
                log("✗ Failed to initialize db_public connection pool", level=1, always=True)
                all_successful = False
        
        conn = db_public.get_connection()
        if conn is None:
            log("✗ db_public.get_connection() returned None", level=1, always=True)
            all_successful = False
        else:
            with conn.cursor() as cursor:
                cursor.execute("SELECT current_database(), current_user;")
                db_name, user = cursor.fetchone()
                log(f"✓ database.db_public connection successful", level=1, always=True)
                log(f"Database: {db_name}", level=2)
                log(f"User: {user}", level=2)
            # Release connection back to pool
            db_public.release_connection(conn)
    except ImportError as e:
        log(f"✗ Could not import database.db_public: {e}", level=1, always=True)
        log("This module may not exist in your codebase", level=2)
    except Exception as e:
        log(f"✗ Error using database.db_public module: {e}", level=1, always=True)
        if VERBOSE:
            traceback.print_exc()
        all_successful = False
    
    return all_successful

def test_application_modules(config):
    """Test application modules that use the database."""
    log_section("APPLICATION MODULES TEST")
    
    if not config:
        log("No configuration available", always=True)
        return False
    
    all_successful = True
    
    # Test database queries directly
    log("Testing direct database queries")
    try:
        from database import db
        
        # Extract Supabase-specific connection parameters
        db_config = None
        if 'database' in config and 'connection' in config['database']:
            connection_config = config['database']['connection']
            if 'supabase' in connection_config:
                db_config = connection_config['supabase'].copy()
                # Add pooling settings if available
                if 'pooling' in config['database']:
                    pooling = config['database']['pooling']
                    if 'min_connections' in pooling:
                        db_config['min_connections'] = pooling['min_connections']
                    if 'max_connections' in pooling:
                        db_config['max_connections'] = pooling['max_connections']
        
        # Initialize the connection pool if needed
        if not hasattr(db, '_pool') or db._pool is None:
            log("Initializing db connection pool", level=1)
            db.init_connection_pool(config=db_config)
        
        conn = db.get_connection()
        if conn:
            # Test basic molecules query
            try:
                with conn.cursor(cursor_factory=RealDictCursor) as cursor:
                    cursor.execute("SELECT count(*) FROM molecules;")
                    result = cursor.fetchone()
                    if result:
                        count = result['count']
                        log(f"✓ Molecule count query successful: {count} molecules", level=1, always=True)
                    
                    # Try to get a specific molecule
                    cursor.execute("SELECT id, name FROM molecules LIMIT 1;")
                    result = cursor.fetchone()
                    if result:
                        molecule_id = result['id']
                        molecule_name = result['name']
                        log(f"✓ Retrieved molecule from database: {molecule_name} (ID: {molecule_id})", level=1, always=True)
                    else:
                        log("No molecules found in database", level=1)
            except Exception as e:
                log(f"✗ Error querying molecules: {e}", level=1, always=True)
                if VERBOSE:
                    traceback.print_exc()
                all_successful = False
            
            # Test querying other tables
            try:
                with conn.cursor(cursor_factory=RealDictCursor) as cursor:
                    # List table names in public schema
                    cursor.execute("""
                        SELECT table_name 
                        FROM information_schema.tables 
                        WHERE table_schema = 'public' 
                        ORDER BY table_name
                    """)
                    tables = [row['table_name'] for row in cursor.fetchall()]
                    log(f"✓ Database tables found: {', '.join(tables[:5])}{'...' if len(tables) > 5 else ''}", 
                       level=1, always=True)
            except Exception as e:
                log(f"✗ Error querying database tables: {e}", level=1, always=True)
                if VERBOSE:
                    traceback.print_exc()
                all_successful = False
            
            # Release the connection
            db.release_connection(conn)
        else:
            log("✗ Failed to get database connection", level=1, always=True)
            all_successful = False
    
    except Exception as e:
        log(f"✗ Error in direct database query test: {e}", level=1, always=True)
        if VERBOSE:
            traceback.print_exc()
        all_successful = False
    
    # Optionally try to import API modules safely
    try:
        log("Testing safe import of API modules", level=1)
        try:
            import api
            log(f"✓ Successfully imported api module", level=1)
            
            # Try to import models safely without executing __init__ code
            import importlib.util
            spec = importlib.util.find_spec('api.models')
            if spec:
                log(f"✓ api.models module exists", level=1)
            else:
                log(f"✗ api.models module not found", level=1)
        except ImportError as e:
            log(f"✗ Error importing api module: {e}", level=1)
    except Exception as e:
        log(f"✗ Error testing API modules: {e}", level=1)
        # Don't affect success status for this optional test
    
    return all_successful

def main():
    """Main entry point for the script."""
    log_section("DATABASE CONNECTION FLOW DEBUGGER")
    log(f"Started at: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}", always=True)
    
    # Find configuration path
    config_path = find_config_path()
    
    # Load configuration
    config = load_config(config_path)
    
    # Extract connection parameters
    params = extract_connection_params(config)
    
    # Test direct connection
    try:
        direct_success = test_connection(params)
    except Exception as e:
        log(f"✗ Unexpected error in direct connection test: {e}", always=True)
        if VERBOSE:
            traceback.print_exc()
        direct_success = False
    
    # Initialize other test results
    adapter_success = False
    module_success = False
    app_success = False
    
    # Only continue if direct connection was successful
    if direct_success:
        # Test adapter connection
        try:
            adapter_success = test_adapter_connection(config, params)
        except Exception as e:
            log(f"✗ Unexpected error in adapter connection test: {e}", always=True)
            if VERBOSE:
                traceback.print_exc()
            adapter_success = False
        
        # Test database modules
        try:
            module_success = test_database_modules(config)
        except Exception as e:
            log(f"✗ Unexpected error in database modules test: {e}", always=True)
            if VERBOSE:
                traceback.print_exc()
            module_success = False
        
        # Test application modules
        try:
            app_success = test_application_modules(config)
        except Exception as e:
            log(f"✗ Unexpected error in application modules test: {e}", always=True)
            if VERBOSE:
                traceback.print_exc()
            app_success = False
    
    # Summary
    log_section("CONNECTION FLOW SUMMARY")
    log(f"Direct connection: {'✓ SUCCESS' if direct_success else '✗ FAILED'}", always=True)
    
    if direct_success:
        log(f"Adapter connection: {'✓ SUCCESS' if adapter_success else '✗ FAILED'}", always=True)
        log(f"Database modules: {'✓ SUCCESS' if module_success else '✗ FAILED'}", always=True)
        log(f"Application modules: {'✓ SUCCESS' if app_success else '✗ FAILED'}", always=True)
        
        if direct_success and adapter_success and module_success and app_success:
            log("\n✓ All connection tests PASSED!", always=True)
            return 0
        else:
            log("\n⚠ Some connection tests FAILED", always=True)
            return 1
    else:
        log("✗ Direct connection failed, cannot continue with further tests", always=True)
        return 1

if __name__ == "__main__":
    sys.exit(main())