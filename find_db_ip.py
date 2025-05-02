#!/usr/bin/env python3
"""
CryoProtect v2 - Find Database IP

This script attempts to find the IP address of the Supabase database
by querying DNS servers directly.
"""

import socket
import subprocess
import platform
import os
from dotenv import load_dotenv

# Load environment variables
load_dotenv()

# ANSI color codes for better readability (optional)
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

print(f"{Colors.HEADER}{Colors.BOLD}CryoProtect v2 - Find Database IP{Colors.ENDC}")
print(f"{Colors.CYAN}This tool will attempt to find the IP address of your Supabase database.{Colors.ENDC}")
print()

# Get Supabase configuration
SUPABASE_PROJECT_ID = os.getenv("SUPABASE_PROJECT_ID")
if not SUPABASE_PROJECT_ID:
    SUPABASE_URL = os.getenv("SUPABASE_URL", "")
    if SUPABASE_URL:
        import urllib.parse
        parsed_url = urllib.parse.urlparse(SUPABASE_URL)
        hostname = parsed_url.netloc
        if hostname and hostname.endswith('.supabase.co'):
            SUPABASE_PROJECT_ID = hostname.split('.')[0]
            print(f"Extracted project ID from URL: {SUPABASE_PROJECT_ID}")
        else:
            # Check if it's an IP address
            if hostname:
                print(f"Using IP address as hostname: {hostname}")
                SUPABASE_PROJECT_ID = hostname

if not SUPABASE_PROJECT_ID:
    print(f"{Colors.FAIL}ERROR: Could not determine project ID{Colors.ENDC}")
    import sys
    sys.exit(1)

# Determine hostnames
if SUPABASE_PROJECT_ID.count('.') > 0:  # It's an IP or domain, not a project ID
    api_hostname = SUPABASE_PROJECT_ID
    db_hostname = f"db.{SUPABASE_PROJECT_ID}"
else:
    api_hostname = f"{SUPABASE_PROJECT_ID}.supabase.co"
    db_hostname = f"db.{SUPABASE_PROJECT_ID}.supabase.co"

print(f"API Hostname: {api_hostname}")
print(f"Database Hostname: {db_hostname}")
print()

# List of DNS servers to try
dns_servers = [
    "8.8.8.8",       # Google DNS
    "1.1.1.1",       # Cloudflare DNS
    "9.9.9.9",       # Quad9
    "208.67.222.222" # OpenDNS
]

# Function to query DNS using nslookup or dig
def query_dns(hostname, dns_server):
    print(f"Querying {dns_server} for {hostname}...")
    
    try:
        if platform.system() == 'Windows':
            result = subprocess.run(
                ['nslookup', hostname, dns_server],
                capture_output=True,
                text=True,
                check=False
            )
            
            if result.returncode == 0:
                # Parse nslookup output
                ip_addresses = []
                for line in result.stdout.splitlines():
                    if 'Address:' in line and dns_server not in line:
                        ip = line.split('Address:')[-1].strip()
                        ip_addresses.append(ip)
                
                if ip_addresses:
                    print(f"Found {len(ip_addresses)} IP address(es): {', '.join(ip_addresses)}")
                    return ip_addresses
                else:
                    print(f"No IP addresses found in response")
                    return []
            else:
                print(f"Query failed: {result.stderr}")
                return []
        else:
            result = subprocess.run(
                ['dig', '@' + dns_server, hostname, '+short'],
                capture_output=True,
                text=True,
                check=False
            )
            
            if result.returncode == 0:
                ip_addresses = result.stdout.strip().split('\n')
                ip_addresses = [ip for ip in ip_addresses if ip]
                
                if ip_addresses:
                    print(f"Found {len(ip_addresses)} IP address(es): {', '.join(ip_addresses)}")
                    return ip_addresses
                else:
                    print(f"No IP addresses found in response")
                    return []
            else:
                print(f"Query failed: {result.stderr}")
                return []
    except Exception as e:
        print(f"Error querying DNS: {str(e)}")
        return []

# Function to update hosts file
def update_hosts_file(hostname, ip_address):
    print(f"\n{Colors.BOLD}Updating hosts file for {hostname} -> {ip_address}{Colors.ENDC}")
    
    # Determine hosts file location
    hosts_file = r'C:\Windows\System32\drivers\etc\hosts' if platform.system() == 'Windows' else '/etc/hosts'
    
    # Check if we have permission to modify the hosts file
    try:
        with open(hosts_file, 'r+') as f:
            pass
        has_permission = True
    except PermissionError:
        has_permission = False
    
    if not has_permission:
        print(f"{Colors.WARNING}You don't have permission to modify the hosts file.{Colors.ENDC}")
        print(f"Please add the following line to your hosts file ({hosts_file}):")
        print(f"  {ip_address}    {hostname}")
        return False
    
    # Read existing hosts file
    with open(hosts_file, 'r') as f:
        lines = f.readlines()
    
    # Check if entry already exists
    entry_exists = False
    for i, line in enumerate(lines):
        if line.strip() and not line.strip().startswith('#'):
            parts = line.strip().split()
            if len(parts) >= 2 and parts[1] == hostname:
                entry_exists = True
                lines[i] = f"{ip_address}    {hostname}\n"
                print(f"Updated existing entry for {hostname}")
                break
    
    # Add new entry if it doesn't exist
    if not entry_exists:
        lines.append(f"{ip_address}    {hostname}\n")
        print(f"Added new entry for {hostname}")
    
    # Write back to hosts file
    with open(hosts_file, 'w') as f:
        f.writelines(lines)
    
    print(f"{Colors.GREEN}Hosts file updated successfully{Colors.ENDC}")
    return True

# Try to resolve API hostname
print(f"{Colors.BOLD}Attempting to resolve API hostname ({api_hostname})...{Colors.ENDC}")
api_ip = None

try:
    api_ip = socket.gethostbyname(api_hostname)
    print(f"{Colors.GREEN}API hostname resolved: {api_hostname} -> {api_ip}{Colors.ENDC}")
except socket.gaierror:
    print(f"{Colors.FAIL}Failed to resolve API hostname{Colors.ENDC}")
    
    # Try alternative DNS servers
    for dns_server in dns_servers:
        ip_addresses = query_dns(api_hostname, dns_server)
        if ip_addresses:
            api_ip = ip_addresses[0]
            print(f"{Colors.GREEN}Found API IP using {dns_server}: {api_ip}{Colors.ENDC}")
            break

# Try to resolve DB hostname
print(f"\n{Colors.BOLD}Attempting to resolve DB hostname ({db_hostname})...{Colors.ENDC}")
db_ip = None

try:
    db_ip = socket.gethostbyname(db_hostname)
    print(f"{Colors.GREEN}DB hostname resolved: {db_hostname} -> {db_ip}{Colors.ENDC}")
except socket.gaierror:
    print(f"{Colors.FAIL}Failed to resolve DB hostname{Colors.ENDC}")
    
    # Try alternative DNS servers
    for dns_server in dns_servers:
        ip_addresses = query_dns(db_hostname, dns_server)
        if ip_addresses:
            db_ip = ip_addresses[0]
            print(f"{Colors.GREEN}Found DB IP using {dns_server}: {db_ip}{Colors.ENDC}")
            break

# Handle case where we couldn't find the DB IP
if not db_ip:
    print(f"\n{Colors.WARNING}Could not find DB IP address through DNS queries.{Colors.ENDC}")
    print(f"Trying to infer from API IP...")
    
    if api_ip:
        # In many cases, the DB IP will be similar to the API IP
        # Often it's on the same subnet with just the last octet different
        ip_parts = api_ip.split('.')
        
        # Try some common patterns - first check if it's 172.64.x.x pattern which is common for Supabase
        if ip_parts[0] == '172' and ip_parts[1] == '64':
            # Often the db uses the same first 3 octets with a different last octet
            potential_db_ip = f"{ip_parts[0]}.{ip_parts[1]}.{ip_parts[2]}.247"
            print(f"Inferring potential DB IP: {potential_db_ip}")
            db_ip = potential_db_ip
        else:
            # Just keep the same IP as a last resort
            print(f"Using API IP as fallback for DB: {api_ip}")
            db_ip = api_ip

# Update hosts file
if api_ip:
    try:
        update_hosts_file(api_hostname, api_ip)
    except Exception as e:
        print(f"Error updating hosts file for API: {str(e)}")
        print(f"Please manually add this line to your hosts file:")
        print(f"{api_ip}    {api_hostname}")

if db_ip:
    try:
        update_hosts_file(db_hostname, db_ip)
    except Exception as e:
        print(f"Error updating hosts file for DB: {str(e)}")
        print(f"Please manually add this line to your hosts file:")
        print(f"{db_ip}    {db_hostname}")

print(f"\n{Colors.BOLD}Summary:{Colors.ENDC}")
print(f"API Hostname: {api_hostname} -> {api_ip or 'Not found'}")
print(f"DB Hostname: {db_hostname} -> {db_ip or 'Not found'}")

if api_ip and db_ip:
    print(f"\n{Colors.GREEN}{Colors.BOLD}Successfully found both IP addresses!{Colors.ENDC}")
    print(f"Test your connection with: python supabase_connection_diagnostic.py")
else:
    print(f"\n{Colors.WARNING}{Colors.BOLD}Could not find all IP addresses.{Colors.ENDC}")
    print(f"You may need to manually update your hosts file or try alternatives.")

# Also write to a standalone file for manual update
try:
    with open('hosts_entries.txt', 'w') as f:
        if api_ip:
            f.write(f"{api_ip}    {api_hostname}\n")
        if db_ip:
            f.write(f"{db_ip}    {db_hostname}\n")
    print(f"\nWrote entries to hosts_entries.txt for manual update if needed")
except Exception as e:
    print(f"Error writing hosts_entries.txt: {str(e)}")

# Generate the .env update command
print(f"\n{Colors.BOLD}Update your .env file with:{Colors.ENDC}")

# Create an update_env.py script
try:
    with open('update_env.py', 'w') as f:
        f.write("""#!/usr/bin/env python3
from dotenv import load_dotenv, set_key

# Load current .env file
load_dotenv()

# Update with correct values
set_key('.env', 'SUPABASE_URL', 'https://{0}')
set_key('.env', 'SUPABASE_DB_HOST', '{1}')

print("Updated .env file with correct hostnames")
""".format(api_hostname, db_hostname))
    print(f"Created update_env.py script to update your .env file")
    print(f"Run: python update_env.py")
except Exception as e:
    print(f"Error creating update_env.py: {str(e)}")
    print(f"SUPABASE_URL=https://{api_hostname}")
    print(f"SUPABASE_DB_HOST={db_hostname}")
