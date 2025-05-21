#!/usr/bin/env python3
"""
CryoProtect v2 - Supabase Connection Diagnostic Tool

This script performs a comprehensive diagnosis of Supabase connection issues,
checking DNS resolution, network connectivity, authentication, and database access.
It provides specific error details and suggests fixes for common issues.
"""

import os
import sys
import json
import time
import socket
import subprocess
import requests
import urllib.parse
from dotenv import load_dotenv
from urllib.parse import urlparse

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

print(f"{Colors.HEADER}{Colors.BOLD}CryoProtect v2 - Supabase Connection Diagnostic Tool{Colors.ENDC}")
print(f"{Colors.CYAN}This tool will diagnose connection issues with your Supabase project.{Colors.ENDC}")
print()

# Extract Supabase configuration
SUPABASE_URL = os.getenv("SUPABASE_URL")
SUPABASE_KEY = os.getenv("SUPABASE_KEY")
SUPABASE_PROJECT_ID = os.getenv("SUPABASE_PROJECT_ID", "")

if not SUPABASE_URL or not SUPABASE_KEY:
    print(f"{Colors.FAIL}ERROR: SUPABASE_URL and SUPABASE_KEY must be set in .env file{Colors.ENDC}")
    sys.exit(1)

print(f"{Colors.BOLD}Supabase Configuration:{Colors.ENDC}")
print(f"  URL: {SUPABASE_URL}")
print(f"  Project ID: {SUPABASE_PROJECT_ID or 'Not specified'}")
print(f"  Key Type: {'anon' if 'anon' in SUPABASE_KEY else 'service_role'}")
print()

# Parse URL for host details
parsed_url = urlparse(SUPABASE_URL)
hostname = parsed_url.netloc
if not hostname:
    print(f"{Colors.FAIL}ERROR: Invalid Supabase URL format. Should be https://your-project-id.supabase.co{Colors.ENDC}")
    sys.exit(1)

# Extract database connection details
DB_HOST = os.getenv("SUPABASE_DB_HOST") or f"db.{hostname}"
DB_PORT = os.getenv("SUPABASE_DB_PORT", "5432")
DB_NAME = os.getenv("SUPABASE_DB_NAME", "postgres")
DB_USER = os.getenv("SUPABASE_DB_USER", "postgres")
DB_PASSWORD = os.getenv("SUPABASE_DB_PASSWORD", "postgres")

# ===== DNS Test =====
def check_dns(hostname):
    print(f"{Colors.BOLD}[1] Checking DNS resolution...{Colors.ENDC}")
    
    # Extract domain from hostname (e.g., "project-id.supabase.co" from "project-id.supabase.co")
    domain = hostname
    if hostname.startswith("db."):
        domain = hostname[3:]  # Remove "db." prefix
    
    # Try to resolve hostname
    try:
        print(f"  Resolving API hostname ({hostname})...")
        ip_address = socket.gethostbyname(hostname)
        print(f"  {Colors.GREEN}Hostname resolved: {hostname} -> {ip_address}{Colors.ENDC}")
    except socket.gaierror as e:
        print(f"  {Colors.FAIL}✗ Failed to resolve hostname: {hostname}{Colors.ENDC}")
        print(f"  Error: {str(e)}")
        
        # Try to resolve the base supabase.co domain
        try:
            print(f"\n  Testing if supabase.co is resolvable...")
            base_ip = socket.gethostbyname("supabase.co")
            print(f"  {Colors.GREEN}✓ Base domain 'supabase.co' resolves to {base_ip}{Colors.ENDC}")
            print(f"  {Colors.WARNING}This suggests a specific project hostname issue, not a general DNS problem.{Colors.ENDC}")
        except socket.gaierror:
            print(f"  {Colors.FAIL}✗ Failed to resolve even 'supabase.co'{Colors.ENDC}")
            print(f"  {Colors.WARNING}This suggests a general DNS resolution problem.{Colors.ENDC}")
        
        # Provide specific recommendations
        print(f"\n  {Colors.BOLD}Recommended fixes:{Colors.ENDC}")
        print(f"  1. Verify your project ID is correct in SUPABASE_URL")
        print(f"  2. Check your internet connection and DNS settings")
        print(f"  3. Try adding the following entry to your hosts file:")
        print(f"     XX.XX.XX.XX {hostname} # Replace XX.XX.XX.XX with actual IP if known")
        
        return False
    
    # Try to resolve DB hostname
    db_hostname = f"db.{domain}"
    try:
        print(f"  Resolving database hostname ({db_hostname})...")
        db_ip_address = socket.gethostbyname(db_hostname)
        print(f"  {Colors.GREEN}✓ Database hostname resolved: {db_hostname} -> {db_ip_address}{Colors.ENDC}")
    except socket.gaierror as e:
        print(f"  {Colors.FAIL}✗ Failed to resolve database hostname: {db_hostname}{Colors.ENDC}")
        print(f"  Error: {str(e)}")
        print(f"  {Colors.WARNING}This may cause database connection issues.{Colors.ENDC}")
        
        return False
    
    return True

# ===== Network Test =====
def check_network(hostname):
    print(f"\n{Colors.BOLD}[2] Checking network connectivity...{Colors.ENDC}")
    
    # Check HTTP connectivity to API
    api_url = f"https://{hostname}/rest/v1/"
    try:
        print(f"  Testing HTTP connection to API ({api_url})...")
        start_time = time.time()
        response = requests.get(api_url, timeout=10)
        elapsed = time.time() - start_time
        
        if response.status_code in [200, 204, 401, 403]:  # These status codes indicate network is fine
            print(f"  {Colors.GREEN}✓ API is reachable (Status: {response.status_code}, Latency: {elapsed:.2f}s){Colors.ENDC}")
        else:
            print(f"  {Colors.WARNING}⚠ API returned unexpected status code: {response.status_code}{Colors.ENDC}")
            print(f"  Response: {response.text[:100]}")
    except requests.exceptions.RequestException as e:
        print(f"  {Colors.FAIL}✗ Failed to connect to API: {str(e)}{Colors.ENDC}")
        
        # Check if it's a certificate error
        if "SSLError" in str(e):
            print(f"  {Colors.WARNING}This appears to be an SSL certificate issue.{Colors.ENDC}")
            print(f"  {Colors.BOLD}Recommended fixes:{Colors.ENDC}")
            print(f"  1. Update your SSL certificates: pip install --upgrade certifi")
            print(f"  2. Check your system time (incorrect time can cause certificate validation failures)")
        elif "ConnectTimeout" in str(e):
            print(f"  {Colors.WARNING}This appears to be a network connectivity or firewall issue.{Colors.ENDC}")
            print(f"  {Colors.BOLD}Recommended fixes:{Colors.ENDC}")
            print(f"  1. Check if your network blocks outbound HTTPS connections")
            print(f"  2. Try connecting from a different network")
        
        return False
    
    # Check HTTP connectivity to database (this should fail with a protocol error, which is normal)
    # We're just checking if the port is reachable
    db_hostname = DB_HOST
    
    try:
        print(f"  Testing TCP connection to database port ({db_hostname}:{DB_PORT})...")
        
        # Try to connect to the database port
        sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        sock.settimeout(5)
        result = sock.connect_ex((db_hostname, int(DB_PORT)))
        sock.close()
        
        if result == 0:
            print(f"  {Colors.GREEN}✓ Database port is reachable{Colors.ENDC}")
        else:
            print(f"  {Colors.FAIL}✗ Database port is not reachable (Error code: {result}){Colors.ENDC}")
            print(f"  {Colors.WARNING}This suggests a network connectivity issue to the database.{Colors.ENDC}")
            print(f"  {Colors.BOLD}Recommended fixes:{Colors.ENDC}")
            print(f"  1. Verify your database hostname and port are correct")
            print(f"  2. Check if your network blocks outbound PostgreSQL connections")
            print(f"  3. For Supabase, direct database connections might require IP allow-listing")
            return False
    except socket.error as e:
        print(f"  {Colors.FAIL}✗ Failed to test database connection: {str(e)}{Colors.ENDC}")
        return False
    
    return True

# ===== Auth Test =====
def check_auth():
    print(f"\n{Colors.BOLD}[3] Checking Supabase authentication...{Colors.ENDC}")
    
    # Test with Supabase key
    headers = {
        "apikey": SUPABASE_KEY,
        "Authorization": f"Bearer {SUPABASE_KEY}",
        "Content-Type": "application/json"
    }
    
    api_url = f"{SUPABASE_URL}/rest/v1/molecules?limit=1"
    
    try:
        print(f"  Testing authentication with service_role key...")
        response = requests.get(api_url, headers=headers)
        
        if response.status_code in [200, 204]:
            print(f"  {Colors.GREEN}✓ Authentication successful{Colors.ENDC}")
            return True
        elif response.status_code == 401:
            print(f"  {Colors.FAIL}✗ Authentication failed: Unauthorized (401){Colors.ENDC}")
            print(f"  {Colors.WARNING}This suggests your API key is invalid or expired.{Colors.ENDC}")
            print(f"  {Colors.BOLD}Recommended fixes:{Colors.ENDC}")
            print(f"  1. Check if your SUPABASE_KEY is correct and hasn't expired")
            print(f"  2. Ensure you're using the service_role key, not the anon key")
            print(f"  3. Regenerate your API key in the Supabase dashboard if necessary")
            
            return False
        elif response.status_code == 403:
            print(f"  {Colors.FAIL}✗ Authentication failed: Forbidden (403){Colors.ENDC}")
            print(f"  {Colors.WARNING}This suggests you don't have permission to access this resource.{Colors.ENDC}")
            print(f"  {Colors.BOLD}Recommended fixes:{Colors.ENDC}")
            print(f"  1. Verify your RLS (Row Level Security) policies allow access")
            print(f"  2. Ensure you're using the service_role key, not the anon key")
            
            return False
        else:
            print(f"  {Colors.FAIL}✗ Authentication failed: Unexpected status code {response.status_code}{Colors.ENDC}")
            print(f"  Response: {response.text[:100]}")
            return False
    except requests.exceptions.RequestException as e:
        print(f"  {Colors.FAIL}✗ Failed to test authentication: {str(e)}{Colors.ENDC}")
        return False

# ===== Database Test =====
def check_database():
    print(f"\n{Colors.BOLD}[4] Checking database access via MCP...{Colors.ENDC}")
    
    try:
        # Use Supabase MCP to execute a simple query
        from supabase_mcp_tools import execute_sql_on_supabase
        
        print(f"  Testing SQL execution via MCP...")
        try:
            result = execute_sql_on_supabase(SUPABASE_PROJECT_ID, "SELECT COUNT(*) FROM public.molecules")
            if result:
                print(f"  {Colors.GREEN}✓ MCP SQL execution successful{Colors.ENDC}")
                print(f"  Result: {json.dumps(result, indent=2)[:100]}...")
                return True
            else:
                print(f"  {Colors.FAIL}✗ MCP SQL execution failed: No result returned{Colors.ENDC}")
                return False
        except Exception as e:
            print(f"  {Colors.FAIL}✗ MCP SQL execution failed: {str(e)}{Colors.ENDC}")
            
            if "command not found" in str(e) or "No such file or directory" in str(e):
                print(f"  {Colors.WARNING}This suggests an issue with the Node.js/npm environment or the MCP package.{Colors.ENDC}")
                print(f"  {Colors.BOLD}Recommended fixes:{Colors.ENDC}")
                print(f"  1. Ensure Node.js and npm are installed: node --version && npm --version")
                print(f"  2. Install the Supabase MCP package: npm install -g @supabase/mcp-server-supabase")
            elif "access token" in str(e).lower() or "authorization" in str(e).lower():
                print(f"  {Colors.WARNING}This suggests an issue with the MCP access token.{Colors.ENDC}")
                print(f"  {Colors.BOLD}Recommended fixes:{Colors.ENDC}")
                print(f"  1. Verify the token in supabase_mcp_tools.py is valid")
                print(f"  2. Re-authenticate with the Supabase CLI: npx supabase login")
            
            return False
    except ImportError:
        print(f"  {Colors.WARNING}⚠ Could not import supabase_mcp_tools module. Skipping MCP test.{Colors.ENDC}")
        return None

# Run all tests
dns_ok = check_dns(hostname)
network_ok = check_network(hostname) if dns_ok else False
auth_ok = check_auth() if network_ok else False
db_ok = check_database() if auth_ok else False

# Summary
print(f"\n{Colors.HEADER}{Colors.BOLD}Diagnostic Summary:{Colors.ENDC}")
print(f"  DNS Resolution:   {'✓' if dns_ok else '✗'}")
print(f"  Network Connectivity: {'✓' if network_ok else '✗'}")
print(f"  Authentication:   {'✓' if auth_ok else '✗'}")
print(f"  Database Access:  {'✓' if db_ok else '?' if db_ok is None else '✗'}")

print(f"\n{Colors.BOLD}Overall Status:{Colors.ENDC}")
if dns_ok and network_ok and auth_ok and (db_ok or db_ok is None):
    print(f"{Colors.GREEN}✓ All checks passed! Your Supabase connection appears to be working correctly.{Colors.ENDC}")
    print(f"\n{Colors.BOLD}Next Steps:{Colors.ENDC}")
    print(f"  - If you're still experiencing issues, check specific API endpoint implementations")
    print(f"  - Verify your application code for handling Supabase responses correctly")
    print(f"  - Consider reviewing your service role configuration or RLS policies")
else:
    print(f"{Colors.FAIL}✗ Some checks failed. See above for specific issues and recommended fixes.{Colors.ENDC}")
    
    print(f"\n{Colors.BOLD}Common Issues and Solutions:{Colors.ENDC}")
    
    if not dns_ok:
        print(f"\n{Colors.BOLD}DNS Issues:{Colors.ENDC}")
        print(f"  - If you're using a VPN, try disabling it temporarily")
        print(f"  - Try using an alternative DNS server like 8.8.8.8 (Google) or 1.1.1.1 (Cloudflare)")
        print(f"  - If on a corporate network, ask your IT team about DNS restrictions")
    
    if not network_ok and dns_ok:
        print(f"\n{Colors.BOLD}Network Issues:{Colors.ENDC}")
        print(f"  - Check firewall settings that might block connections to *.supabase.co")
        print(f"  - If you're on a restricted network, request outbound access to Supabase servers")
        print(f"  - Test your connection with a simple curl or wget command")
    
    if not auth_ok and network_ok:
        print(f"\n{Colors.BOLD}Authentication Issues:{Colors.ENDC}")
        print(f"  - Verify your API key in your .env file matches the key in your Supabase dashboard")
        print(f"  - Check if your project is in a paused state (free tier projects pause after inactivity)")
        print(f"  - Make sure you're using the service_role key for administrative operations")
    
    if not db_ok and auth_ok:
        print(f"\n{Colors.BOLD}Database Issues:{Colors.ENDC}")
        print(f"  - Check if your tables exist in the database")
        print(f"  - Verify your RLS policies allow the operation you're trying to perform")
        print(f"  - If using MCP, ensure Node.js, npm, and the Supabase MCP package are installed")

print(f"\n{Colors.CYAN}For further assistance, check the Supabase documentation at https://supabase.com/docs{Colors.ENDC}")

if __name__ == "__main__":
    pass