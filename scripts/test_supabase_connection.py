#!/usr/bin/env python3
"""
Comprehensive Supabase Connection Test Script for CryoProtect

This script tests:
1. IPv4 and IPv6 connectivity to Supabase
2. Database connection using both protocols
3. Connection pooling with different configurations
4. API access using both protocols

Usage:
  python test_supabase_connection.py --url SUPABASE_URL --key SUPABASE_KEY [--anon-key ANON_KEY] [--prefer-ipv4]

Requirements:
  pip install psycopg2-binary requests supabase
"""

import argparse
import json
import os
import socket
import sys
import time
from urllib.parse import urlparse

try:
    import psycopg2
    from psycopg2 import pool
except ImportError:
    print("Error: psycopg2 module not found. Install with: pip install psycopg2-binary")
    sys.exit(1)

try:
    import requests
except ImportError:
    print("Error: requests module not found. Install with: pip install requests")
    sys.exit(1)

try:
    from supabase import create_client, Client
except ImportError:
    print("Warning: supabase-py module not found. Install with: pip install supabase")
    has_supabase_client = False
else:
    has_supabase_client = True


# ANSI color codes
class Colors:
    RED = '\033[0;31m'
    GREEN = '\033[0;32m'
    YELLOW = '\033[0;33m'
    BLUE = '\033[0;34m'
    MAGENTA = '\033[0;35m'
    CYAN = '\033[0;36m'
    BOLD = '\033[1m'
    NC = '\033[0m'  # No Color


def header(text):
    """Print a section header"""
    print(f"\n{Colors.BOLD}{Colors.MAGENTA}● {text}{Colors.NC}")
    print(f"{Colors.MAGENTA}{'-' * 50}{Colors.NC}")


def success(text):
    """Print a success message"""
    print(f"{Colors.GREEN}✓ {text}{Colors.NC}")


def warning(text):
    """Print a warning message"""
    print(f"{Colors.YELLOW}⚠ {text}{Colors.NC}")


def error(text):
    """Print an error message"""
    print(f"{Colors.RED}✗ {text}{Colors.NC}")


def info(text):
    """Print an info message"""
    print(f"{Colors.CYAN}{text}{Colors.NC}")


def resolve_hostname(hostname, family=socket.AF_UNSPEC):
    """Resolve hostname to IP address(es)"""
    try:
        addrinfo = socket.getaddrinfo(hostname, None, family, socket.SOCK_STREAM)
        return [addr[4][0] for addr in addrinfo]
    except socket.gaierror as e:
        return []


def test_dns_resolution(hostname):
    """Test DNS resolution for both IPv4 and IPv6"""
    header("DNS Resolution Test")
    
    ipv4_addresses = resolve_hostname(hostname, socket.AF_INET)
    ipv6_addresses = resolve_hostname(hostname, socket.AF_INET6)
    
    if ipv4_addresses:
        success(f"IPv4 resolution: {', '.join(ipv4_addresses)}")
    else:
        warning("No IPv4 addresses found")
    
    if ipv6_addresses:
        success(f"IPv6 resolution: {', '.join(ipv6_addresses)}")
    else:
        warning("No IPv6 addresses found")
    
    return ipv4_addresses, ipv6_addresses


def test_basic_connectivity(hostname, port=5432):
    """Test basic connectivity to hostname:port via both IPv4 and IPv6"""
    header("Basic Connectivity Test")
    
    ipv4_ok = False
    ipv6_ok = False
    
    # Test IPv4 connectivity
    ipv4_addresses = resolve_hostname(hostname, socket.AF_INET)
    if ipv4_addresses:
        try:
            s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
            s.settimeout(5)
            start_time = time.time()
            s.connect((ipv4_addresses[0], port))
            elapsed = time.time() - start_time
            s.close()
            success(f"IPv4 connection to {ipv4_addresses[0]}:{port} succeeded in {elapsed:.3f}s")
            ipv4_ok = True
        except Exception as e:
            error(f"IPv4 connection to {ipv4_addresses[0]}:{port} failed: {e}")
    
    # Test IPv6 connectivity
    ipv6_addresses = resolve_hostname(hostname, socket.AF_INET6)
    if ipv6_addresses:
        try:
            s = socket.socket(socket.AF_INET6, socket.SOCK_STREAM)
            s.settimeout(5)
            start_time = time.time()
            s.connect((ipv6_addresses[0], port))
            elapsed = time.time() - start_time
            s.close()
            success(f"IPv6 connection to {ipv6_addresses[0]}:{port} succeeded in {elapsed:.3f}s")
            ipv6_ok = True
        except Exception as e:
            error(f"IPv6 connection to {ipv6_addresses[0]}:{port} failed: {e}")
    
    return ipv4_ok, ipv6_ok


def test_postgres_connection(hostname, port, dbname, password, prefer_ipv4=True):
    """Test PostgreSQL connections using different strategies"""
    header("PostgreSQL Connection Test")
    
    results = {
        "ipv4_direct": None,
        "ipv6_direct": None,
        "hostname": None,
        "pooled": None
    }
    
    # Get IP addresses
    ipv4_addresses = resolve_hostname(hostname, socket.AF_INET)
    ipv6_addresses = resolve_hostname(hostname, socket.AF_INET6)
    
    # Test direct IPv4 connection
    if ipv4_addresses:
        try:
            conn_str = f"host={ipv4_addresses[0]} port={port} dbname={dbname} user=postgres password={password} sslmode=require"
            start_time = time.time()
            conn = psycopg2.connect(conn_str)
            elapsed = time.time() - start_time
            
            with conn.cursor() as cursor:
                cursor.execute("SELECT version()")
                version = cursor.fetchone()[0]
            
            conn.close()
            success(f"Direct IPv4 connection succeeded in {elapsed:.3f}s")
            info(f"  PostgreSQL version: {version.split(',')[0]}")
            results["ipv4_direct"] = True
        except Exception as e:
            error(f"Direct IPv4 connection failed: {e}")
            results["ipv4_direct"] = False
    
    # Test direct IPv6 connection
    if ipv6_addresses:
        try:
            # Note the square brackets around IPv6 address
            conn_str = f"host={ipv6_addresses[0]} port={port} dbname={dbname} user=postgres password={password} sslmode=require"
            start_time = time.time()
            conn = psycopg2.connect(conn_str)
            elapsed = time.time() - start_time()
            
            with conn.cursor() as cursor:
                cursor.execute("SELECT version()")
                version = cursor.fetchone()[0]
            
            conn.close()
            success(f"Direct IPv6 connection succeeded in {elapsed:.3f}s")
            info(f"  PostgreSQL version: {version.split(',')[0]}")
            results["ipv6_direct"] = True
        except Exception as e:
            error(f"Direct IPv6 connection failed: {e}")
            results["ipv6_direct"] = False
    
    # Test connection with hostname (DNS resolution)
    try:
        conn_str = f"host={hostname} port={port} dbname={dbname} user=postgres password={password} sslmode=require"
        start_time = time.time()
        conn = psycopg2.connect(conn_str)
        elapsed = time.time() - start_time
        
        with conn.cursor() as cursor:
            cursor.execute("SELECT version()")
            version = cursor.fetchone()[0]
            
            # Check which IP was actually used
            cursor.execute("SELECT inet_server_addr()")
            server_addr = cursor.fetchone()[0]
        
        conn.close()
        success(f"Hostname connection succeeded in {elapsed:.3f}s")
        info(f"  Connected to: {server_addr}")
        info(f"  PostgreSQL version: {version.split(',')[0]}")
        results["hostname"] = True
    except Exception as e:
        error(f"Hostname connection failed: {e}")
        results["hostname"] = False
    
    # Test connection pooling
    try:
        # Create a connection pool
        info("Testing connection pool...")
        
        # Choose the appropriate address family based on preference and availability
        if prefer_ipv4 and ipv4_addresses:
            host = ipv4_addresses[0]
            info(f"  Creating pool with IPv4 address: {host}")
        elif ipv6_addresses:
            host = ipv6_addresses[0]
            info(f"  Creating pool with IPv6 address: {host}")
        else:
            host = hostname
            info(f"  Creating pool with hostname: {host}")
        
        # Create the connection pool
        connpool = pool.ThreadedConnectionPool(
            minconn=2,
            maxconn=5,
            host=host,
            port=port,
            dbname=dbname,
            user="postgres",
            password=password,
            sslmode="require"
        )
        
        # Test getting connections from the pool
        connections = []
        for i in range(3):
            start_time = time.time()
            conn = connpool.getconn()
            elapsed = time.time() - start_time
            connections.append(conn)
            success(f"  Pool connection {i+1} acquired in {elapsed:.3f}s")
        
        # Return connections to the pool
        for conn in connections:
            connpool.putconn(conn)
        
        # Close the pool
        connpool.closeall()
        success("Connection pool test succeeded")
        results["pooled"] = True
    except Exception as e:
        error(f"Connection pool test failed: {e}")
        results["pooled"] = False
    
    return results


def test_api_access(url, anon_key):
    """Test API access using both IPv4 and IPv6"""
    header("API Access Test")
    
    # Extract hostname from URL
    parsed_url = urlparse(url)
    hostname = parsed_url.netloc
    
    api_url = f"{url}/rest/v1/"
    headers = {
        "apikey": anon_key,
        "Content-Type": "application/json"
    }
    
    # Test IPv4 API access
    try:
        start_time = time.time()
        response = requests.get(api_url, headers=headers, family=socket.AF_INET)
        elapsed = time.time() - start_time
        
        if response.status_code in (200, 401, 404):  # Any of these could be valid based on permissions
            success(f"IPv4 API access succeeded in {elapsed:.3f}s (Status: {response.status_code})")
            ipv4_api_ok = True
        else:
            error(f"IPv4 API access returned unexpected status: {response.status_code}")
            ipv4_api_ok = False
    except Exception as e:
        error(f"IPv4 API access failed: {e}")
        ipv4_api_ok = False
    
    # Test IPv6 API access
    try:
        start_time = time.time()
        response = requests.get(api_url, headers=headers, family=socket.AF_INET6)
        elapsed = time.time() - start_time
        
        if response.status_code in (200, 401, 404):  # Any of these could be valid based on permissions
            success(f"IPv6 API access succeeded in {elapsed:.3f}s (Status: {response.status_code})")
            ipv6_api_ok = True
        else:
            error(f"IPv6 API access returned unexpected status: {response.status_code}")
            ipv6_api_ok = False
    except Exception as e:
        error(f"IPv6 API access failed: {e}")
        ipv6_api_ok = False
    
    # Test Supabase client if available
    if has_supabase_client:
        info("Testing Supabase client...")
        try:
            start_time = time.time()
            supabase: Client = create_client(url, anon_key)
            elapsed = time.time() - start_time
            
            # Try a simple query to verify connection
            response = supabase.table("molecules").select("count").execute()
            
            success(f"Supabase client connection succeeded in {elapsed:.3f}s")
            client_ok = True
        except Exception as e:
            error(f"Supabase client connection failed: {e}")
            client_ok = False
    else:
        warning("Supabase client not available for testing")
        client_ok = None
    
    return {
        "ipv4_api": ipv4_api_ok,
        "ipv6_api": ipv6_api_ok,
        "client": client_ok
    }


def summarize_results(url, dns_results, connectivity_results, postgres_results, api_results):
    """Summarize all test results and provide recommendations"""
    header("Summary and Recommendations")
    
    # Extract hostname from URL
    parsed_url = urlparse(url)
    hostname = parsed_url.netloc
    
    # Count successful tests for each IP version
    ipv4_count = sum([
        bool(dns_results[0]),               # Has IPv4 DNS entries
        connectivity_results[0],            # Basic IPv4 connectivity
        postgres_results.get("ipv4_direct", False),  # Direct IPv4 PostgreSQL
        api_results.get("ipv4_api", False)  # IPv4 API access
    ])
    
    ipv6_count = sum([
        bool(dns_results[1]),               # Has IPv6 DNS entries
        connectivity_results[1],            # Basic IPv6 connectivity
        postgres_results.get("ipv6_direct", False),  # Direct IPv6 PostgreSQL
        api_results.get("ipv6_api", False)  # IPv6 API access
    ])
    
    # Overall assessment
    print(f"{Colors.BOLD}Overall Network Compatibility:{Colors.NC}")
    print(f"  IPv4 Compatibility: {ipv4_count}/4 tests passed")
    print(f"  IPv6 Compatibility: {ipv6_count}/4 tests passed")
    
    # Connection recommendations
    print(f"\n{Colors.BOLD}Connection Strategy Recommendation:{Colors.NC}")
    
    if ipv4_count >= 3 and ipv6_count >= 3:
        # Both IPv4 and IPv6 work well
        success("Dual-stack connectivity works well")
        print(f"Recommendation: Use flexible connection handler that can use both IPv4 and IPv6")
        print("""
Implementation example:
```python
def get_connection(hostname, prefer_ipv4=True):
    # Try preferred protocol first, then fall back
    families = [socket.AF_INET, socket.AF_INET6]
    if not prefer_ipv4:
        families.reverse()
    
    for family in families:
        try:
            addrinfo = socket.getaddrinfo(hostname, 5432, family, socket.SOCK_STREAM)
            if addrinfo:
                ip_addr = addrinfo[0][4][0]
                conn = psycopg2.connect(
                    host=ip_addr,
                    port=5432,
                    dbname="postgres",
                    user="postgres",
                    password="your_password",
                    sslmode="require"
                )
                return conn
        except Exception:
            continue
    
    # Fall back to hostname
    return psycopg2.connect(
        host=hostname,
        port=5432,
        dbname="postgres",
        user="postgres",
        password="your_password",
        sslmode="require"
    )
```
""")
    elif ipv4_count >= 3:
        # IPv4 works better
        warning("IPv4 connectivity works best")
        print(f"Recommendation: Configure for IPv4-preferred connections")
        print("""
Implementation example:
```python
def get_connection(hostname):
    # Force IPv4 connection
    try:
        addrinfo = socket.getaddrinfo(hostname, 5432, socket.AF_INET, socket.SOCK_STREAM)
        if addrinfo:
            ip_addr = addrinfo[0][4][0]
            return psycopg2.connect(
                host=ip_addr,
                port=5432,
                dbname="postgres",
                user="postgres",
                password="your_password",
                sslmode="require"
            )
    except Exception:
        # Fall back to hostname which will still attempt IPv4 first
        return psycopg2.connect(
            host=hostname,
            port=5432,
            dbname="postgres",
            user="postgres",
            password="your_password",
            sslmode="require"
        )
```
""")
    elif ipv6_count >= 3:
        # IPv6 works better
        warning("IPv6 connectivity works best")
        print(f"Recommendation: Configure for IPv6-preferred connections")
        print("""
Implementation example:
```python
def get_connection(hostname):
    # Force IPv6 connection
    try:
        addrinfo = socket.getaddrinfo(hostname, 5432, socket.AF_INET6, socket.SOCK_STREAM)
        if addrinfo:
            ip_addr = addrinfo[0][4][0]
            return psycopg2.connect(
                host=ip_addr,
                port=5432,
                dbname="postgres",
                user="postgres",
                password="your_password",
                sslmode="require"
            )
    except Exception:
        # Fall back to hostname
        return psycopg2.connect(
            host=hostname,
            port=5432,
            dbname="postgres",
            user="postgres",
            password="your_password",
            sslmode="require"
        )
```
""")
    else:
        # Network connectivity issues
        error("Limited connectivity detected")
        print(f"Recommendation: Check network configuration and firewall settings")
        print("""
Troubleshooting steps:
1. Verify network connectivity to Supabase services
2. Check firewall rules for ports 5432 and 443
3. Ensure DNS resolution is working correctly
4. Try explicitly setting IPv4 or IPv6 in connection strings
""")
    
    # Configuration file recommendation
    print(f"\n{Colors.BOLD}Configuration Example:{Colors.NC}")
    print("""
Add this to your config.py file:

```python
# Database configuration
DB_CONFIG = {
    # Connection details
    'host': 'yourproject.supabase.co',  # Will be resolved based on preference
    'port': 5432,
    'dbname': 'postgres',
    'user': 'postgres',
    'password': 'your_db_password',
    'sslmode': 'require',
    
    # Protocol preferences
    'prefer_ipv4': True,  # Set to False to prefer IPv6
    'fallback_to_alternate': True,  # Try alternate protocol if preferred fails
    
    # Connection pooling settings
    'min_connections': 5,
    'max_connections': 20,
    'connection_timeout': 30,
    'idle_in_transaction_timeout': 60,
    'validation_query': 'SELECT 1',
    
    # Error handling
    'max_retries': 3,
    'retry_delay': 1,  # seconds
}
```
""")


def main():
    """Main function to run the tests"""
    parser = argparse.ArgumentParser(description='Test Supabase connectivity')
    parser.add_argument('--url', required=True, help='Supabase URL (e.g., https://yourproject.supabase.co)')
    parser.add_argument('--key', required=True, help='Supabase database password or service role key')
    parser.add_argument('--anon-key', help='Supabase anonymous key for API testing')
    parser.add_argument('--prefer-ipv4', action='store_true', help='Prefer IPv4 over IPv6')
    parser.add_argument('--dbname', default='postgres', help='Database name (default: postgres)')
    
    args = parser.parse_args()
    
    # Extract hostname from URL
    parsed_url = urlparse(args.url)
    hostname = parsed_url.netloc
    
    print(f"{Colors.BOLD}{Colors.BLUE}============================================={Colors.NC}")
    print(f"{Colors.BOLD}{Colors.BLUE}  CryoProtect Supabase Connection Test       {Colors.NC}")
    print(f"{Colors.BOLD}{Colors.BLUE}============================================={Colors.NC}")
    print(f"Testing connectivity to: {args.url}")
    print(f"Database: {args.dbname}")
    print(f"Protocol preference: {'IPv4' if args.prefer_ipv4 else 'IPv6'}")
    print()
    
    # Run tests
    dns_results = test_dns_resolution(hostname)
    connectivity_results = test_basic_connectivity(hostname)
    postgres_results = test_postgres_connection(hostname, 5432, args.dbname, args.key, args.prefer_ipv4)
    
    # API tests are optional and require the anon key
    api_results = {}
    if args.anon_key:
        api_results = test_api_access(args.url, args.anon_key)
    else:
        warning("Skipping API tests (--anon-key not provided)")
    
    # Summarize results
    summarize_results(args.url, dns_results, connectivity_results, postgres_results, api_results)
    
    print(f"\n{Colors.BOLD}{Colors.BLUE}============================================={Colors.NC}")
    print(f"{Colors.BOLD}{Colors.BLUE}       Connection Test Completed               {Colors.NC}")
    print(f"{Colors.BOLD}{Colors.BLUE}============================================={Colors.NC}")


if __name__ == "__main__":
    main()