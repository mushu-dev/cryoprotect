#!/usr/bin/env python3
"""
CryoProtect v2 - Fix Supabase DNS and Connection Issues

This script attempts to fix common Supabase DNS and connection issues by:
1. Resolving hostnames through fallback DNS servers
2. Creating a local hosts file entry (Windows/WSL compatible)
3. Testing connection with the new DNS entries
4. Updating the environment with the direct IP if needed
"""

import os
import sys
import socket
import platform
import subprocess
import requests
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

print(f"{Colors.HEADER}{Colors.BOLD}CryoProtect v2 - Fix Supabase DNS and Connection Issues{Colors.ENDC}")
print(f"{Colors.CYAN}This tool will attempt to fix DNS and connection issues with your Supabase project.{Colors.ENDC}")
print()

# Extract Supabase URL
SUPABASE_URL = os.getenv("SUPABASE_URL")
if not SUPABASE_URL:
    print(f"{Colors.FAIL}ERROR: SUPABASE_URL must be set in .env file{Colors.ENDC}")
    sys.exit(1)

# Parse URL for hostname
import urllib.parse
parsed_url = urllib.parse.urlparse(SUPABASE_URL)
hostname = parsed_url.netloc
if not hostname:
    print(f"{Colors.FAIL}ERROR: Invalid Supabase URL format. Should be https://your-project-id.supabase.co{Colors.ENDC}")
    sys.exit(1)

# Database hostname
db_hostname = f"db.{hostname}"

# Print current settings
print(f"{Colors.BOLD}Current Supabase Settings:{Colors.ENDC}")
print(f"  API Hostname: {hostname}")
print(f"  Database Hostname: {db_hostname}")
print()

# Check if running in WSL
def is_wsl():
    try:
        with open('/proc/version', 'r') as f:
            return 'microsoft' in f.read().lower()
    except:
        return False

# Function to check if we have admin/root privileges
def has_admin_privileges():
    if platform.system() == 'Windows':
        import ctypes
        return ctypes.windll.shell32.IsUserAnAdmin() != 0
    else:
        return os.geteuid() == 0

# Function to resolve hostname with fallback DNS servers
def resolve_hostname_with_fallbacks(hostname):
    print(f"{Colors.BOLD}Attempting to resolve {hostname} using fallback DNS servers...{Colors.ENDC}")
    
    # List of fallback DNS servers to try
    dns_servers = [
        "8.8.8.8",       # Google Public DNS
        "1.1.1.1",       # Cloudflare DNS
        "9.9.9.9",       # Quad9
        "208.67.222.222" # OpenDNS
    ]
    
    # Try standard resolution first
    try:
        ip = socket.gethostbyname(hostname)
        print(f"  {Colors.GREEN}✓ Hostname resolved using system DNS: {hostname} -> {ip}{Colors.ENDC}")
        return ip
    except socket.gaierror:
        print(f"  {Colors.WARNING}⚠ Failed to resolve using system DNS{Colors.ENDC}")
    
    # If that fails, try each fallback DNS server
    for dns_server in dns_servers:
        print(f"  Trying DNS server {dns_server}...")
        
        try:
            # Use nslookup or dig to resolve hostname using specific DNS server
            if platform.system() == 'Windows':
                result = subprocess.run(
                    ['nslookup', hostname, dns_server],
                    capture_output=True,
                    text=True,
                    check=False
                )
                
                # Parse nslookup output
                for line in result.stdout.splitlines():
                    if 'Address:' in line and not dns_server in line:
                        ip = line.split('Address:')[-1].strip()
                        print(f"  {Colors.GREEN}✓ Hostname resolved using {dns_server}: {hostname} -> {ip}{Colors.ENDC}")
                        return ip
            else:
                result = subprocess.run(
                    ['dig', '@' + dns_server, hostname, '+short'],
                    capture_output=True,
                    text=True,
                    check=False
                )
                
                if result.stdout.strip():
                    ip = result.stdout.strip().split('\n')[0]
                    print(f"  {Colors.GREEN}✓ Hostname resolved using {dns_server}: {hostname} -> {ip}{Colors.ENDC}")
                    return ip
                
        except Exception as e:
            print(f"  Failed with {dns_server}: {str(e)}")
    
    # If we got here, all DNS servers failed
    print(f"  {Colors.FAIL}✗ Failed to resolve hostname with all DNS servers{Colors.ENDC}")
    return None

# Function to add entries to hosts file
def add_to_hosts_file(hostnames_with_ips):
    """Add hostname->IP mappings to the hosts file"""
    
    if not hostnames_with_ips:
        print(f"{Colors.WARNING}No hostnames to add to hosts file{Colors.ENDC}")
        return False
    
    # Determine hosts file location
    hosts_file = r'C:\Windows\System32\drivers\etc\hosts' if platform.system() == 'Windows' else '/etc/hosts'
    
    if is_wsl():
        print(f"{Colors.BOLD}Running in WSL, will modify both WSL and Windows hosts files{Colors.ENDC}")
        wsl_hosts_file = '/etc/hosts'
        win_hosts_file = r'C:\Windows\System32\drivers\etc\hosts'
    else:
        print(f"{Colors.BOLD}Modifying hosts file: {hosts_file}{Colors.ENDC}")
    
    if not has_admin_privileges():
        print(f"{Colors.WARNING}⚠ You don't have administrator privileges. Manual hosts file modification required.{Colors.ENDC}")
        print(f"\n{Colors.BOLD}Please add the following entries to your hosts file ({hosts_file}):{Colors.ENDC}")
        for hostname, ip in hostnames_with_ips.items():
            print(f"  {ip}    {hostname}")
        
        if platform.system() == 'Windows':
            print(f"\nOn Windows, you can do this by:")
            print(f"1. Run Notepad as administrator")
            print(f"2. Open file {hosts_file}")
            print(f"3. Add the lines above to the end of the file")
            print(f"4. Save and close Notepad")
        else:
            print(f"\nOn Linux/MacOS, you can do this by:")
            print(f"1. Run: sudo nano {hosts_file}")
            print(f"2. Add the lines above to the end of the file")
            print(f"3. Press Ctrl+O to save, then Ctrl+X to exit")
        
        user_input = input(f"\n{Colors.BOLD}Press Enter once you've manually updated your hosts file, or type 'skip' to continue without updating: {Colors.ENDC}")
        return user_input.lower() != 'skip'
    
    # We have admin privileges, so we can modify the hosts file
    try:
        # Read existing hosts file
        with open(hosts_file, 'r') as f:
            content = f.readlines()
        
        # Check if entries already exist
        existing_entries = {}
        for line in content:
            if line.strip() and not line.strip().startswith('#'):
                parts = line.strip().split()
                if len(parts) >= 2:
                    existing_entries[parts[1]] = parts[0]
        
        # Add new entries
        new_content = content.copy()
        added_entries = []
        
        for hostname, ip in hostnames_with_ips.items():
            if hostname in existing_entries:
                if existing_entries[hostname] == ip:
                    print(f"  Entry already exists: {ip}    {hostname}")
                else:
                    print(f"  Updating entry: {existing_entries[hostname]} -> {ip}    {hostname}")
                    # Remove existing entry
                    new_content = [line for line in new_content if not (not line.strip().startswith('#') and hostname in line.split())]
                    # Add new entry
                    new_content.append(f"{ip}    {hostname}\n")
                    added_entries.append(hostname)
            else:
                print(f"  Adding new entry: {ip}    {hostname}")
                new_content.append(f"{ip}    {hostname}\n")
                added_entries.append(hostname)
        
        if not added_entries:
            print(f"  {Colors.GREEN}✓ No changes needed to hosts file{Colors.ENDC}")
            return True
        
        # Write back to hosts file
        with open(hosts_file, 'w') as f:
            f.writelines(new_content)
        
        print(f"  {Colors.GREEN}✓ Successfully updated hosts file{Colors.ENDC}")
        
        # If running in WSL, also update Windows hosts file
        if is_wsl():
            print(f"  {Colors.BOLD}Updating Windows hosts file...{Colors.ENDC}")
            
            # Create a PowerShell script to modify the Windows hosts file
            ps_script = f"""
$hostFile = '{win_hosts_file}'
$content = Get-Content $hostFile
$newContent = $content

"""
            for hostname, ip in hostnames_with_ips.items():
                ps_script += f"""
# Check for and remove any existing entry for {hostname}
$newContent = $newContent | Where-Object {{ $_ -notmatch '{hostname}\\s*$' }}
# Add new entry
$newContent += "{ip}    {hostname}"
"""
            
            ps_script += """
# Write back to hosts file
$newContent | Set-Content $hostFile
"""
            
            # Write PowerShell script to a temporary file
            with open('/tmp/update_hosts.ps1', 'w') as f:
                f.write(ps_script)
            
            # Execute the PowerShell script as administrator
            try:
                subprocess.run(
                    ['powershell.exe', '-ExecutionPolicy', 'Bypass', '-Command', 'Start-Process powershell -ArgumentList "-ExecutionPolicy Bypass -File /tmp/update_hosts.ps1" -Verb RunAs'],
                    check=True
                )
                print(f"  {Colors.GREEN}✓ Successfully updated Windows hosts file{Colors.ENDC}")
            except Exception as e:
                print(f"  {Colors.FAIL}✗ Failed to update Windows hosts file: {str(e)}{Colors.ENDC}")
                print(f"  {Colors.WARNING}You may need to manually update your Windows hosts file.{Colors.ENDC}")
        
        return True
        
    except Exception as e:
        print(f"  {Colors.FAIL}✗ Failed to modify hosts file: {str(e)}{Colors.ENDC}")
        return False

# Function to test connection after DNS changes
def test_connection(hostname):
    print(f"{Colors.BOLD}Testing connection to {hostname}...{Colors.ENDC}")
    
    try:
        # Try to resolve hostname again to verify hosts file changes
        ip = socket.gethostbyname(hostname)
        print(f"  {Colors.GREEN}✓ Hostname now resolves: {hostname} -> {ip}{Colors.ENDC}")
        
        # For API hostname, test HTTP connection
        if not hostname.startswith('db.'):
            url = f"https://{hostname}/rest/v1/"
            try:
                response = requests.get(url, timeout=5)
                print(f"  {Colors.GREEN}✓ HTTP connection successful (Status: {response.status_code}){Colors.ENDC}")
                return True
            except requests.exceptions.RequestException as e:
                print(f"  {Colors.FAIL}✗ HTTP connection failed: {str(e)}{Colors.ENDC}")
                return False
        # For DB hostname, test socket connection
        else:
            sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
            sock.settimeout(5)
            result = sock.connect_ex((hostname, 5432))
            sock.close()
            
            if result == 0:
                print(f"  {Colors.GREEN}✓ Socket connection successful{Colors.ENDC}")
                return True
            else:
                print(f"  {Colors.FAIL}✗ Socket connection failed (Error code: {result}){Colors.ENDC}")
                return False
    except socket.gaierror:
        print(f"  {Colors.FAIL}✗ Hostname still doesn't resolve{Colors.ENDC}")
        return False

# Function to update .env file with direct IP
def update_env_with_direct_ip(hostname, ip):
    print(f"{Colors.BOLD}Updating .env file with direct IP for {hostname}...{Colors.ENDC}")
    
    try:
        # Read .env file
        with open('.env', 'r') as f:
            content = f.readlines()
        
        # Update or add entry
        updated = False
        for i, line in enumerate(content):
            if line.strip().startswith('SUPABASE_URL='):
                protocol = 'https://'
                path = ''
                
                # Parse current URL
                current_url = line.strip().split('=', 1)[1].strip()
                if current_url.startswith(protocol):
                    current_url = current_url[len(protocol):]
                
                # Check if there's a path component
                if '/' in current_url:
                    current_hostname, path = current_url.split('/', 1)
                    path = '/' + path
                else:
                    current_hostname = current_url
                
                # Replace hostname with IP but keep protocol and path
                if current_hostname == hostname:
                    new_url = f"{protocol}{ip}{path}"
                    content[i] = f"SUPABASE_URL={new_url}\n"
                    updated = True
                    print(f"  {Colors.GREEN}✓ Updated SUPABASE_URL to use direct IP: {new_url}{Colors.ENDC}")
                    break
        
        if not updated:
            print(f"  {Colors.WARNING}⚠ Did not find SUPABASE_URL entry in .env file{Colors.ENDC}")
            return False
        
        # Add comment explaining the change
        content.insert(content.index(f"SUPABASE_URL={new_url}\n"), 
                      f"# Direct IP workaround for DNS issues with {hostname}\n")
        
        # Write back to .env file
        with open('.env', 'w') as f:
            f.writelines(content)
        
        print(f"  {Colors.GREEN}✓ Successfully updated .env file{Colors.ENDC}")
        return True
        
    except Exception as e:
        print(f"  {Colors.FAIL}✗ Failed to update .env file: {str(e)}{Colors.ENDC}")
        return False

# Main function
def main():
    # Step 1: Resolve hostnames with fallback DNS
    api_ip = resolve_hostname_with_fallbacks(hostname)
    db_ip = resolve_hostname_with_fallbacks(db_hostname)
    
    if not api_ip or not db_ip:
        print(f"{Colors.FAIL}Failed to resolve one or both hostnames. Cannot proceed with fix.{Colors.ENDC}")
        return
    
    # Step 2: Add entries to hosts file
    hosts_entries = {
        hostname: api_ip,
        db_hostname: db_ip
    }
    
    hosts_updated = add_to_hosts_file(hosts_entries)
    
    # Step 3: Test connection
    if hosts_updated:
        print(f"\n{Colors.BOLD}Testing connections after hosts file update...{Colors.ENDC}")
        api_connection = test_connection(hostname)
        db_connection = test_connection(db_hostname)
        
        if api_connection and db_connection:
            print(f"\n{Colors.GREEN}{Colors.BOLD}✓ DNS issue fixed successfully!{Colors.ENDC}")
            print(f"{Colors.GREEN}Both API and database connections are now working.{Colors.ENDC}")
            return
        else:
            print(f"\n{Colors.WARNING}⚠ Hosts file updated, but connection test failed.{Colors.ENDC}")
    
    # Step 4: If hosts file update didn't work or wasn't done, try direct IP in .env
    print(f"\n{Colors.BOLD}Attempting direct IP workaround...{Colors.ENDC}")
    env_updated = update_env_with_direct_ip(hostname, api_ip)
    
    if env_updated:
        print(f"\n{Colors.GREEN}{Colors.BOLD}✓ Direct IP workaround applied!{Colors.ENDC}")
        print(f"{Colors.GREEN}Your application will now use the direct IP instead of the hostname.{Colors.ENDC}")
        print(f"{Colors.WARNING}Note: This is a temporary workaround. Consider fixing your DNS configuration.{Colors.ENDC}")
    else:
        print(f"\n{Colors.FAIL}✗ Failed to apply direct IP workaround.{Colors.ENDC}")
    
    # Final instructions
    print(f"\n{Colors.BOLD}Restart your application to apply changes.{Colors.ENDC}")
    print(f"If issues persist, consider:")
    print(f"1. Using a different network or VPN")
    print(f"2. Contacting your network administrator")
    print(f"3. Using Supabase CLI with npx -y @supabase/mcp-server-supabase@latest")

if __name__ == "__main__":
    main()