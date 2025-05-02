#!/usr/bin/env python3
"""
IP Resolution Module for CryoProtect v2

This module provides functions to resolve hostnames to IP addresses using multiple methods:
1. Standard DNS resolution
2. Alternative DNS servers (8.8.8.8, 1.1.1.1)
3. MCP resolution (placeholder)
4. Heuristic method (placeholder)

The module handles OS-specific differences for DNS resolution commands.
"""

import socket
import subprocess
import platform
import logging
import re
import os
from typing import Optional, List, Dict
from dotenv import load_dotenv

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Alternative DNS servers to try
ALTERNATIVE_DNS_SERVERS = ["8.8.8.8", "1.1.1.1"]

def resolve_ip_address(hostname: str) -> Optional[str]:
    """
    Attempt to resolve IP address using multiple methods.
    
    Args:
        hostname: The hostname to resolve
        
    Returns:
        The resolved IP address as a string, or None if resolution fails
    """
    # Method 1: Standard resolution
    ip_address = _standard_dns_resolution(hostname)
    if ip_address:
        return ip_address
    
    # Method 2: Alternative DNS servers
    ip_address = _alternative_dns_resolution(hostname)
    if ip_address:
        return ip_address
    
    # Method 3: MCP Resolution (placeholder)
    ip_address = _mcp_resolution(hostname)
    if ip_address:
        return ip_address
    
    # Method 4: Heuristic fallback (placeholder)
    ip_address = _heuristic_resolution(hostname)
    if ip_address:
        return ip_address
    
    # All methods failed
    logger.error(f"All IP resolution methods failed for hostname: {hostname}")
    return None

def _standard_dns_resolution(hostname: str) -> Optional[str]:
    """
    Attempt to resolve hostname using standard socket.gethostbyname().
    
    Args:
        hostname: The hostname to resolve
        
    Returns:
        The resolved IP address as a string, or None if resolution fails
    """
    try:
        ip_address = socket.gethostbyname(hostname)
        logger.info(f"Standard DNS resolution succeeded: {hostname} -> {ip_address}")
        return ip_address
    except socket.gaierror as e:
        logger.warning(f"Standard DNS resolution failed for {hostname}: {str(e)}")
        return None
    except Exception as e:
        logger.warning(f"Unexpected error during standard DNS resolution for {hostname}: {str(e)}")
        return None

def _alternative_dns_resolution(hostname: str) -> Optional[str]:
    """
    Attempt to resolve hostname using alternative DNS servers.
    Uses OS-specific commands (nslookup for Windows, dig for Unix).
    
    Args:
        hostname: The hostname to resolve
        
    Returns:
        The resolved IP address as a string, or None if resolution fails
    """
    for dns_server in ALTERNATIVE_DNS_SERVERS:
        try:
            # Windows vs Unix implementation differences
            if platform.system() == "Windows":
                ip_address = _resolve_with_nslookup(hostname, dns_server)
            else:
                ip_address = _resolve_with_dig(hostname, dns_server)
                
            if ip_address:
                logger.info(f"Alternative DNS resolution succeeded using {dns_server}: {hostname} -> {ip_address}")
                return ip_address
        except Exception as e:
            logger.warning(f"Alternative DNS resolution failed using {dns_server} for {hostname}: {str(e)}")
            continue
    
    logger.warning(f"All alternative DNS resolution attempts failed for {hostname}")
    return None

def _resolve_with_nslookup(hostname: str, dns_server: str) -> Optional[str]:
    """
    Resolve hostname using nslookup command (Windows).
    
    Args:
        hostname: The hostname to resolve
        dns_server: The DNS server to use
        
    Returns:
        The resolved IP address as a string, or None if resolution fails
    """
    try:
        result = subprocess.run(
            ["nslookup", hostname, dns_server],
            capture_output=True, text=True, timeout=5
        )
        
        # Check if the command was successful
        if result.returncode != 0:
            logger.warning(f"nslookup command failed with return code {result.returncode}")
            return None
        
        # Parse the output to extract IP address
        # nslookup output format example:
        # Server:  dns.google
        # Address:  8.8.8.8
        #
        # Name:    example.com
        # Address:  93.184.216.34
        
        output_lines = result.stdout.strip().split('\n')
        ip_pattern = re.compile(r'Address:\s+(\d+\.\d+\.\d+\.\d+)')
        
        # Skip the first IP address (which is the DNS server's address)
        found_first_address = False
        
        for line in output_lines:
            match = ip_pattern.search(line)
            if match:
                if found_first_address:
                    return match.group(1)
                found_first_address = True
        
        logger.warning(f"No IP address found in nslookup output for {hostname}")
        return None
    except subprocess.TimeoutExpired:
        logger.warning(f"nslookup command timed out for {hostname} using {dns_server}")
        return None
    except Exception as e:
        logger.warning(f"Error executing nslookup for {hostname} using {dns_server}: {str(e)}")
        return None

def _resolve_with_dig(hostname: str, dns_server: str) -> Optional[str]:
    """
    Resolve hostname using dig command (Unix).
    
    Args:
        hostname: The hostname to resolve
        dns_server: The DNS server to use
        
    Returns:
        The resolved IP address as a string, or None if resolution fails
    """
    try:
        result = subprocess.run(
            ["dig", "+short", hostname, f"@{dns_server}"],
            capture_output=True, text=True, timeout=5
        )
        
        # Check if the command was successful
        if result.returncode != 0:
            logger.warning(f"dig command failed with return code {result.returncode}")
            return None
        
        # Parse the output to extract IP address
        # dig +short output format is typically just the IP address on a line
        output = result.stdout.strip()
        if output:
            # Take the first line if there are multiple results
            ip_address = output.split('\n')[0]
            # Verify it's a valid IP address
            if re.match(r'^\d+\.\d+\.\d+\.\d+$', ip_address):
                return ip_address
        
        logger.warning(f"No valid IP address found in dig output for {hostname}")
        return None
    except subprocess.TimeoutExpired:
        logger.warning(f"dig command timed out for {hostname} using {dns_server}")
        return None
    except Exception as e:
        logger.warning(f"Error executing dig for {hostname} using {dns_server}: {str(e)}")
        return None

def _mcp_resolution(hostname: str) -> Optional[str]:
    """
    Placeholder for MCP-based resolution.
    This would use the Supabase MCP to get project info including IP.
    
    Args:
        hostname: The hostname to resolve
        
    Returns:
        The resolved IP address as a string, or None if resolution fails
    """
    # This is a placeholder implementation
    # In a real implementation, this would use the Supabase MCP to get project info
    logger.info(f"MCP resolution not fully implemented for {hostname}")
    
    # Extract project ID from hostname
    # Example: db.tsdlmynydfuypiugmkev.supabase.co -> tsdlmynydfuypiugmkev
    project_id_match = re.search(r'db\.([a-z0-9]+)\.supabase\.co', hostname)
    if not project_id_match:
        logger.warning(f"Could not extract project ID from hostname: {hostname}")
        return None
    
    project_id = project_id_match.group(1)
    logger.info(f"Extracted project ID from hostname: {project_id}")
    
    # In a real implementation, we would use the Supabase MCP to get the project info
    # For now, we'll just return None to indicate this method failed
    return None

def _heuristic_resolution(hostname: str) -> Optional[str]:
    """
    Placeholder for heuristic-based resolution.
    This would use patterns in IP assignment as a fallback.
    
    Args:
        hostname: The hostname to resolve
        
    Returns:
        The resolved IP address as a string, or None if resolution fails
    """
    # This is a placeholder implementation
    # In a real implementation, this might use known IP patterns for Supabase projects
    logger.info(f"Heuristic resolution not fully implemented for {hostname}")
    
    # Extract project ID from hostname
    # Example: db.tsdlmynydfuypiugmkev.supabase.co -> tsdlmynydfuypiugmkev
    project_id_match = re.search(r'db\.([a-z0-9]+)\.supabase\.co', hostname)
    if not project_id_match:
        logger.warning(f"Could not extract project ID from hostname: {hostname}")
        return None
    
    # In a real implementation, we might have some heuristics based on project ID
    # For now, we'll just return None to indicate this method failed
    return None

def update_env_with_ip(hostname: str, ip_address: str) -> bool:
    """
    Update the .env file with the resolved IP address.
    
    Args:
        hostname: The hostname that was resolved
        ip_address: The resolved IP address
        
    Returns:
        True if the .env file was updated successfully, False otherwise
    """
    try:
        # Load current .env file
        load_dotenv()
        
        # Check if .env file exists
        env_path = '.env'
        if not os.path.exists(env_path):
            logger.warning(f".env file not found at {env_path}")
            return False
        
        # Read current .env file
        with open(env_path, 'r') as f:
            env_lines = f.readlines()
        
        # Check if SUPABASE_DB_IP already exists
        ip_var_exists = False
        for i, line in enumerate(env_lines):
            if line.startswith('SUPABASE_DB_IP='):
                # Update existing variable
                env_lines[i] = f'SUPABASE_DB_IP={ip_address}\n'
                ip_var_exists = True
                break
        
        # Add new variable if it doesn't exist
        if not ip_var_exists:
            env_lines.append(f'\n# Added by IP Resolution Module\n')
            env_lines.append(f'SUPABASE_DB_IP={ip_address}\n')
        
        # Write updated .env file
        with open(env_path, 'w') as f:
            f.writelines(env_lines)
        
        logger.info(f"Updated .env file with SUPABASE_DB_IP={ip_address}")
        return True
    except Exception as e:
        logger.error(f"Failed to update .env file: {str(e)}")
        return False

def resolve_and_update_env(hostname: str) -> Dict[str, str]:
    """
    Resolve hostname to IP address and update .env file.
    
    Args:
        hostname: The hostname to resolve
        
    Returns:
        Dictionary with results: {'hostname': hostname, 'ip': ip_address, 'updated': True/False}
    """
    result = {
        'hostname': hostname,
        'ip': None,
        'updated': False
    }
    
    # Resolve IP address
    ip_address = resolve_ip_address(hostname)
    if not ip_address:
        logger.error(f"Failed to resolve IP address for {hostname}")
        return result
    
    result['ip'] = ip_address
    
    # Update .env file
    updated = update_env_with_ip(hostname, ip_address)
    result['updated'] = updated
    
    return result

if __name__ == "__main__":
    # Example usage
    import sys
    
    if len(sys.argv) > 1:
        hostname = sys.argv[1]
    else:
        # Default to Supabase hostname from project
        hostname = "db.tsdlmynydfuypiugmkev.supabase.co"
    
    print(f"Resolving IP address for {hostname}...")
    result = resolve_and_update_env(hostname)
    
    if result['ip']:
        print(f"✅ Successfully resolved {hostname} to {result['ip']}")
        if result['updated']:
            print(f"✅ Updated .env file with SUPABASE_DB_IP={result['ip']}")
        else:
            print(f"❌ Failed to update .env file")
    else:
        print(f"❌ Failed to resolve IP address for {hostname}")