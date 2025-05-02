"""
Supabase Direct Connection Diagnostic Tool

This script attempts to connect to a Supabase PostgreSQL database using different methods
to diagnose connection issues and provides a solution.
"""

import os
import sys
import socket
import logging
import psycopg2
from dotenv import load_dotenv
from typing import Dict, Any, Optional, Tuple

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

def load_env_variables() -> Dict[str, Any]:
    """Load environment variables and return them as a dictionary."""
    load_dotenv()
    
    return {
        'db_host': os.getenv('SUPABASE_DB_HOST'),
        'db_port': os.getenv('SUPABASE_DB_PORT', '5432'),
        'db_name': os.getenv('SUPABASE_DB_NAME', 'postgres'),
        'db_user': os.getenv('SUPABASE_DB_USER', 'postgres'),
        'db_password': os.getenv('SUPABASE_DB_PASSWORD'),
        'db_ip': os.getenv('SUPABASE_DB_IP_ADDRESS'),
        'project_id': os.getenv('SUPABASE_PROJECT_ID')
    }

def resolve_hostname(hostname: str) -> Tuple[Optional[str], str]:
    """
    Attempt to resolve a hostname to an IP address and return details.
    
    Returns:
        Tuple containing (resolved_ip, message)
    """
    try:
        ip_address = socket.gethostbyname(hostname)
        return ip_address, f"Successfully resolved {hostname} to {ip_address}"
    except socket.gaierror as e:
        return None, f"Failed to resolve {hostname}: {str(e)}"

def test_connection_with_hostname(env_vars: Dict[str, Any]) -> Tuple[bool, str]:
    """
    Test database connection using hostname.
    
    Returns:
        Tuple containing (success, message)
    """
    try:
        conn_string = (
            f"host={env_vars['db_host']} "
            f"port={env_vars['db_port']} "
            f"dbname={env_vars['db_name']} "
            f"user={env_vars['db_user']} "
            f"password={env_vars['db_password']} "
            f"connect_timeout=10"
        )
        
        logger.info(f"Attempting to connect to {env_vars['db_host']}:{env_vars['db_port']} "
                  f"as {env_vars['db_user']} using hostname")
        
        conn = psycopg2.connect(conn_string)
        with conn.cursor() as cursor:
            cursor.execute("SELECT 1")
            result = cursor.fetchone()
        
        conn.close()
        
        if result and result[0] == 1:
            return True, "Successfully connected using hostname"
        else:
            return False, "Connected but query execution failed"
    
    except Exception as e:
        return False, f"Failed to connect using hostname: {str(e)}"

def test_connection_with_ip(env_vars: Dict[str, Any], ip_address: str) -> Tuple[bool, str]:
    """
    Test database connection using IP address.
    
    Returns:
        Tuple containing (success, message)
    """
    try:
        conn_string = (
            f"host={ip_address} "
            f"port={env_vars['db_port']} "
            f"dbname={env_vars['db_name']} "
            f"user={env_vars['db_user']} "
            f"password={env_vars['db_password']} "
            f"connect_timeout=10"
        )
        
        logger.info(f"Attempting to connect to {ip_address}:{env_vars['db_port']} "
                  f"as {env_vars['db_user']} using IP address")
        
        conn = psycopg2.connect(conn_string)
        with conn.cursor() as cursor:
            cursor.execute("SELECT 1")
            result = cursor.fetchone()
        
        conn.close()
        
        if result and result[0] == 1:
            return True, f"Successfully connected using IP address {ip_address}"
        else:
            return False, "Connected but query execution failed"
    
    except Exception as e:
        return False, f"Failed to connect using IP address {ip_address}: {str(e)}"

def test_connection_with_supabase_cli_settings() -> Tuple[bool, str]:
    """
    Test database connection using Supabase CLI connection settings.
    
    Returns:
        Tuple containing (success, message)
    """
    try:
        # This uses the default Supabase local development settings
        conn_string = (
            "host=127.0.0.1 "
            "port=54322 "
            "dbname=postgres "
            "user=postgres "
            "password=postgres "
            "connect_timeout=10"
        )
        
        logger.info("Attempting to connect to local Supabase instance at 127.0.0.1:54322")
        
        conn = psycopg2.connect(conn_string)
        with conn.cursor() as cursor:
            cursor.execute("SELECT 1")
            result = cursor.fetchone()
        
        conn.close()
        
        if result and result[0] == 1:
            return True, "Successfully connected to local Supabase instance"
        else:
            return False, "Connected to local Supabase but query execution failed"
    
    except Exception as e:
        return False, f"Failed to connect to local Supabase instance: {str(e)}"

def generate_connection_url(host: str, port: str, dbname: str, user: str, password: str) -> str:
    """Generate a PostgreSQL connection URL."""
    # URL encode the password to handle special characters
    import urllib.parse
    encoded_password = urllib.parse.quote_plus(password)
    
    return f"postgresql://{user}:{encoded_password}@{host}:{port}/{dbname}"

def get_alternative_ip_addresses(project_id: str) -> Dict[str, str]:
    """
    Get alternative IP addresses that might work for Supabase connection.
    
    Args:
        project_id: The Supabase project ID
        
    Returns:
        Dictionary of connection type to IP address
    """
    # These are common Cloudflare IP ranges used by Supabase
    return {
        "Cloudflare IPv4 (Primary)": "172.64.149.246",
        "Cloudflare IPv4 (Alternate 1)": "172.64.128.246",
        "Cloudflare IPv4 (Alternate 2)": "172.64.148.246",
        "IPv4 Direct (Last Resort)": "52.71.23.87"
    }

def main():
    """Main function to diagnose Supabase connection issues."""
    print("\n=====================================================")
    print("Supabase Direct Connection Diagnostic Tool")
    print("=====================================================\n")
    
    # Load environment variables
    env_vars = load_env_variables()
    
    # Check required variables
    missing_vars = []
    for key in ['db_host', 'db_user', 'db_password']:
        if not env_vars.get(key):
            missing_vars.append(key.upper())
    
    if missing_vars:
        print(f"Error: Missing required environment variables: {', '.join(missing_vars)}")
        print("Please check your .env file and ensure these variables are set.")
        return 1
    
    print("Environment variables loaded successfully.")
    
    # Hostname Resolution Test
    print("\n----- Step 1: Testing Hostname Resolution -----")
    ip_address, resolution_message = resolve_hostname(env_vars['db_host'])
    print(resolution_message)
    
    # Standard hostname-based connection test
    print("\n----- Step 2: Testing Connection with Hostname -----")
    hostname_success, hostname_message = test_connection_with_hostname(env_vars)
    print(hostname_message)
    
    # IP-based connection test
    if ip_address:
        print("\n----- Step 3: Testing Connection with Resolved IP -----")
        ip_success, ip_message = test_connection_with_ip(env_vars, ip_address)
        print(ip_message)
    else:
        ip_success = False
        print("\nSkipping IP-based connection test due to hostname resolution failure.")
    
    # Try configured IP address if available
    if not ip_success and env_vars.get('db_ip'):
        print(f"\n----- Step 4: Testing Connection with Configured IP ({env_vars['db_ip']}) -----")
        configured_ip_success, configured_ip_message = test_connection_with_ip(env_vars, env_vars['db_ip'])
        print(configured_ip_message)
        ip_success = configured_ip_success
    
    # Try alternative IP addresses
    if not hostname_success and not ip_success and env_vars.get('project_id'):
        print("\n----- Step 5: Testing Connection with Alternative IPs -----")
        alternative_ips = get_alternative_ip_addresses(env_vars['project_id'])
        
        for ip_type, alt_ip in alternative_ips.items():
            print(f"\nTrying {ip_type}: {alt_ip}")
            alt_success, alt_message = test_connection_with_ip(env_vars, alt_ip)
            print(alt_message)
            
            if alt_success:
                print(f"\nSUCCESS! Found working IP address: {alt_ip}")
                print("To fix your connection, add the following to your .env file:")
                print(f"SUPABASE_DB_IP_ADDRESS={alt_ip}")
                
                # Update environment file
                with open(".env", "r") as f:
                    env_content = f.read()
                
                if "SUPABASE_DB_IP_ADDRESS=" not in env_content:
                    with open(".env", "a") as f:
                        f.write(f"\n# Added by Connection Diagnostic Tool\nSUPABASE_DB_IP_ADDRESS={alt_ip}\n")
                    print("\nThe .env file has been updated with the working IP address.")
                else:
                    print("\nYour .env file already contains a SUPABASE_DB_IP_ADDRESS entry. Please update it manually.")
                
                break
    
    # Try Supabase CLI local settings
    print("\n----- Step 6: Testing Connection with Supabase CLI Local Settings -----")
    cli_success, cli_message = test_connection_with_supabase_cli_settings()
    print(cli_message)
    
    # Summary and recommendations
    print("\n=====================================================")
    print("Connection Diagnostic Summary")
    print("=====================================================")
    
    if hostname_success:
        print("\n✅ Hostname-based connection successful!")
        print("You can use the standard hostname-based connection method.")
        connection_url = generate_connection_url(
            env_vars['db_host'], 
            env_vars['db_port'], 
            env_vars['db_name'], 
            env_vars['db_user'], 
            env_vars['db_password']
        )
        print(f"\nConnection URL: {connection_url}")
    elif ip_success:
        print("\n✅ IP-based connection successful!")
        print("You should use the IP-based connection method.")
        
        if ip_address:
            working_ip = ip_address
        elif env_vars.get('db_ip'):
            working_ip = env_vars['db_ip']
        else:
            working_ip = "your_working_ip"  # Placeholder, should not reach here
        
        connection_url = generate_connection_url(
            working_ip, 
            env_vars['db_port'], 
            env_vars['db_name'], 
            env_vars['db_user'], 
            env_vars['db_password']
        )
        print(f"\nConnection URL with IP: {connection_url}")
    elif cli_success:
        print("\n✅ Local Supabase CLI connection successful!")
        print("You can use the local Supabase instance for development.")
        connection_url = "postgresql://postgres:postgres@127.0.0.1:54322/postgres"
        print(f"\nLocal Connection URL: {connection_url}")
    else:
        print("\n❌ All connection attempts failed.")
        print("\nRecommendations:")
        print("1. Check if your Supabase project is active and not paused")
        print("2. Verify your IP is allowed in the Supabase project settings")
        print("3. Try connecting with the Supabase CLI using 'npx supabase start'")
        print("4. Contact Supabase support for additional help")
    
    print("\nDiagnostic complete!")
    return 0

if __name__ == "__main__":
    sys.exit(main())