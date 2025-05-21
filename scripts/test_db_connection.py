#!/usr/bin/env python
"""
Test direct connection to Supabase PostgreSQL database from Fedora
"""
import os
import sys
import time
import argparse
from urllib.parse import urlparse
from dotenv import load_dotenv
import psycopg2

# ANSI colors
RED = '\033[0;31m'
GREEN = '\033[0;32m'
YELLOW = '\033[0;33m'
BLUE = '\033[0;34m'
RESET = '\033[0m'

def print_colored(color, message):
    """Print a colored message"""
    print(f"{color}{message}{RESET}")

def print_section(title):
    """Print a section title"""
    print("\n" + "=" * 60)
    print_colored(BLUE, f"  {title}")
    print("=" * 60)

def load_environment():
    """Load environment variables from .env file"""
    env_path = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), '.env')
    
    if not os.path.exists(env_path):
        print_colored(YELLOW, f"Warning: .env file not found at {env_path}")
        return {}
    
    load_dotenv(env_path)
    return {
        'supabase_url': os.getenv('SUPABASE_URL'),
        'supabase_key': os.getenv('SUPABASE_KEY'),
        'supabase_service_key': os.getenv('SUPABASE_SERVICE_KEY'),
        'supabase_db_password': os.getenv('SUPABASE_DB_PASSWORD'),
        'database_url': os.getenv('DATABASE_URL')
    }

def extract_db_connection_info(supabase_url, database_url=None):
    """Extract database connection info from Supabase URL"""
    # Parse Supabase URL to get hostname
    if not supabase_url:
        return None
        
    parsed_url = urlparse(supabase_url)
    hostname = parsed_url.netloc
    
    # If database_url is provided, parse it
    if database_url:
        try:
            parsed_db_url = urlparse(database_url)
            username = parsed_db_url.username
            password = parsed_db_url.password
            db_name = parsed_db_url.path.lstrip('/')
            port = parsed_db_url.port or 5432
            return {
                'hostname': parsed_db_url.hostname,
                'port': port,
                'username': username,
                'password': password,
                'db_name': db_name
            }
        except Exception as e:
            print_colored(YELLOW, f"Warning: Could not parse DATABASE_URL: {e}")
    
    # Default Supabase configuration
    return {
        'hostname': hostname,
        'port': 5432,
        'username': 'postgres',
        'password': None,  # Will be prompted for if not provided
        'db_name': 'postgres'
    }

def test_connection(conn_info, password):
    """Test connection to PostgreSQL database"""
    if not conn_info:
        print_colored(RED, "Error: Missing connection information")
        return False
        
    if not password and not conn_info.get('password'):
        print_colored(RED, "Error: No password provided")
        return False
    
    conn_password = conn_info.get('password') or password
    
    # Print connection info
    print_colored(BLUE, "Connection Information:")
    print(f"  Hostname: {conn_info['hostname']}")
    print(f"  Port: {conn_info['port']}")
    print(f"  Username: {conn_info['username']}")
    print(f"  Database: {conn_info['db_name']}")
    print(f"  Password: {'*' * 8 if conn_password else 'Not provided'}")
    
    try:
        print("\nAttempting to connect...")
        start_time = time.time()
        
        conn = psycopg2.connect(
            host=conn_info['hostname'],
            port=conn_info['port'],
            user=conn_info['username'],
            password=conn_password,
            dbname=conn_info['db_name'],
            connect_timeout=10
        )
        
        elapsed_time = time.time() - start_time
        print_colored(GREEN, f"✓ Connection successful! (Connected in {elapsed_time:.2f}s)")
        
        # Get PostgreSQL version
        cursor = conn.cursor()
        cursor.execute("SELECT version();")
        version = cursor.fetchone()[0]
        print_colored(BLUE, "PostgreSQL Server Information:")
        print(f"  {version}")
        
        # Test a simple query
        print("\nExecuting simple query...")
        cursor.execute("SELECT current_database(), current_user, inet_server_addr(), inet_server_port();")
        db_info = cursor.fetchone()
        
        print_colored(BLUE, "Connection Details:")
        print(f"  Current Database: {db_info[0]}")
        print(f"  Current User: {db_info[1]}")
        print(f"  Server Address: {db_info[2]}")
        print(f"  Server Port: {db_info[3]}")
        
        # Close connection
        cursor.close()
        conn.close()
        print_colored(GREEN, "Connection closed successfully")
        return True
        
    except psycopg2.OperationalError as e:
        print_colored(RED, f"✗ Connection failed: {e}")
        
        # Suggest possible fixes
        if "timeout" in str(e).lower():
            print_colored(YELLOW, "\nPossible reasons for timeout:")
            print("  1. Network connectivity issues")
            print("  2. Firewall blocking the connection")
            print("  3. Database server is not running or unreachable")
            print("  4. Wrong hostname or port")
            
        elif "password authentication failed" in str(e).lower():
            print_colored(YELLOW, "\nPossible reasons for authentication failure:")
            print("  1. Incorrect password")
            print("  2. User doesn't have access to this database")
            print("  3. User requires SSL connection")
            
        else:
            print_colored(YELLOW, "\nPossible fixes:")
            print("  1. Check your network connectivity")
            print("  2. Verify hostname and port")
            print("  3. Check your firewall settings")
            print("  4. Try connecting with SSL (sslmode=require)")
            
        return False
        
    except Exception as e:
        print_colored(RED, f"✗ Unexpected error: {e}")
        return False

def main():
    """Main function"""
    parser = argparse.ArgumentParser(description='Test Supabase PostgreSQL connection')
    parser.add_argument('-p', '--password', help='Database password')
    parser.add_argument('-u', '--url', help='Supabase URL')
    args = parser.parse_args()
    
    print_section("Supabase PostgreSQL Connection Test")
    print("Testing direct connection to Supabase PostgreSQL database from Fedora")
    
    # Load environment variables
    env_vars = load_environment()
    
    supabase_url = args.url or env_vars.get('supabase_url')
    
    if not supabase_url:
        print_colored(RED, "Error: Supabase URL not provided")
        print("Please provide a Supabase URL using --url or add it to your .env file")
        return 1
        
    # Extract connection info from Supabase URL
    conn_info = extract_db_connection_info(supabase_url, env_vars.get('database_url'))
    
    # Get password
    password = args.password or env_vars.get('supabase_db_password')
    
    if not password:
        print_colored(YELLOW, "No password provided in arguments or .env file")
        password = input("Enter database password: ")
    
    # Test connection
    if test_connection(conn_info, password):
        print_section("Connection Test: SUCCESS")
        print_colored(GREEN, "Your Fedora environment can successfully connect to the Supabase PostgreSQL database!")
        print("\nYou're ready to run the CryoProtect application with Podman:")
        print("  cd /home/mushu/Projects/CryoProtect && ./quickstart_podman.sh")
        return 0
    else:
        print_section("Connection Test: FAILED")
        print_colored(RED, "Unable to connect to the Supabase PostgreSQL database")
        print("\nTry the following:")
        print("1. Check your Supabase database connection settings")
        print("2. Verify network connectivity with: ./scripts/test_supabase_connectivity.sh")
        print("3. Ensure your firewall allows outbound connections to port 5432")
        print("4. Check if your Supabase database is running and accessible")
        return 1

if __name__ == '__main__':
    sys.exit(main())