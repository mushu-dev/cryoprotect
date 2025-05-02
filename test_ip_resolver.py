#!/usr/bin/env python3
"""
Test script for the IP Resolution module.
This script tests the various resolution methods implemented in ip_resolver.py.
"""

import logging
import sys
from ip_resolver import (
    resolve_ip_address,
    _standard_dns_resolution,
    _alternative_dns_resolution,
    resolve_and_update_env
)

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def test_standard_resolution(hostname):
    """Test standard DNS resolution."""
    print("\n=== Testing Standard DNS Resolution ===")
    ip_address = _standard_dns_resolution(hostname)
    if ip_address:
        print(f"✅ Standard resolution succeeded: {hostname} -> {ip_address}")
    else:
        print(f"❌ Standard resolution failed for {hostname}")
    return ip_address

def test_alternative_resolution(hostname):
    """Test alternative DNS resolution."""
    print("\n=== Testing Alternative DNS Resolution ===")
    ip_address = _alternative_dns_resolution(hostname)
    if ip_address:
        print(f"✅ Alternative resolution succeeded: {hostname} -> {ip_address}")
    else:
        print(f"❌ Alternative resolution failed for {hostname}")
    return ip_address

def test_all_methods(hostname):
    """Test all resolution methods."""
    print("\n=== Testing All Resolution Methods ===")
    ip_address = resolve_ip_address(hostname)
    if ip_address:
        print(f"✅ Resolution succeeded: {hostname} -> {ip_address}")
    else:
        print(f"❌ Resolution failed for {hostname}")
    return ip_address

def test_env_update(hostname):
    """Test resolving and updating .env file."""
    print("\n=== Testing .env Update ===")
    result = resolve_and_update_env(hostname)
    if result['ip']:
        print(f"✅ Resolution succeeded: {hostname} -> {result['ip']}")
        if result['updated']:
            print(f"✅ Updated .env file with SUPABASE_DB_IP={result['ip']}")
        else:
            print(f"❌ Failed to update .env file")
    else:
        print(f"❌ Resolution failed for {hostname}")
    return result

def main():
    """Run all tests."""
    print("IP Resolution Module Test")
    print("========================\n")
    
    # Get hostname from command line or use default
    if len(sys.argv) > 1:
        hostname = sys.argv[1]
    else:
        # Default to Supabase hostname from project
        hostname = "db.tsdlmynydfuypiugmkev.supabase.co"
    
    print(f"Testing with hostname: {hostname}\n")
    
    # Test each method individually
    std_ip = test_standard_resolution(hostname)
    alt_ip = test_alternative_resolution(hostname)
    all_ip = test_all_methods(hostname)
    env_result = test_env_update(hostname)
    
    # Print summary
    print("\nTest Summary:")
    print(f"  Standard Resolution:    {'✅' if std_ip else '❌'}")
    print(f"  Alternative Resolution: {'✅' if alt_ip else '❌'}")
    print(f"  All Methods:            {'✅' if all_ip else '❌'}")
    print(f"  .env Update:            {'✅' if env_result['updated'] else '❌'}")
    
    # Test with a known resolvable hostname as a control
    print("\nControl Test with google.com:")
    control_ip = resolve_ip_address("google.com")
    print(f"  google.com Resolution:  {'✅' if control_ip else '❌'}")
    if control_ip:
        print(f"  Resolved to: {control_ip}")

if __name__ == "__main__":
    main()