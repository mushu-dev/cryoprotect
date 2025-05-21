#!/usr/bin/env python
"""
Simple connectivity test for Supabase without psycopg2 dependency
"""
import os
import sys
import json
import socket
import time
import argparse
from urllib.parse import urlparse
from dotenv import load_dotenv
import requests

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
    }

def test_socket_connection(host, port=5432):
    """Test basic socket connection to a host:port"""
    try:
        start_time = time.time()
        sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        sock.settimeout(5)
        result = sock.connect_ex((host, port))
        elapsed_time = time.time() - start_time
        
        if result == 0:
            print_colored(GREEN, f"✓ Socket connection successful to {host}:{port} ({elapsed_time:.2f}s)")
            return True
        else:
            print_colored(RED, f"✗ Socket connection failed to {host}:{port} - Error code: {result}")
            return False
    except socket.gaierror:
        print_colored(RED, f"✗ Hostname could not be resolved: {host}")
        return False
    except socket.error as e:
        print_colored(RED, f"✗ Socket error: {e}")
        return False
    finally:
        sock.close()

def test_http_connection(url, api_key):
    """Test HTTP connection to Supabase REST API"""
    try:
        headers = {
            "apikey": api_key,
            "Authorization": f"Bearer {api_key}"
        }
        
        start_time = time.time()
        response = requests.get(f"{url}/rest/v1/", headers=headers, timeout=5)
        elapsed_time = time.time() - start_time
        
        print_colored(BLUE, f"HTTP Response Status: {response.status_code}")
        print(f"Response Time: {elapsed_time:.2f}s")
        
        if response.status_code < 300:
            print_colored(GREEN, "✓ Successfully connected to Supabase REST API")
            return True
        elif response.status_code == 401:
            print_colored(YELLOW, "⚠ Connected to API, but authentication failed (check your API key)")
            return False
        else:
            print_colored(RED, f"✗ API connection failed - Status Code: {response.status_code}")
            try:
                print(f"Response: {json.dumps(response.json(), indent=2)}")
            except:
                print(f"Response: {response.text[:200]}...")
            return False
    except requests.RequestException as e:
        print_colored(RED, f"✗ HTTP request failed: {e}")
        return False

def main():
    """Main function"""
    parser = argparse.ArgumentParser(description='Test Supabase connection without psycopg2')
    parser.add_argument('-u', '--url', help='Supabase URL')
    parser.add_argument('-k', '--key', help='Supabase API key')
    args = parser.parse_args()
    
    print_section("Supabase Connection Test (Simple)")
    
    # Load environment variables
    env_vars = load_environment()
    
    supabase_url = args.url or env_vars.get('supabase_url')
    supabase_key = args.key or env_vars.get('supabase_key')
    
    if not supabase_url:
        print_colored(RED, "Error: Supabase URL not provided")
        print("Please provide a Supabase URL using --url or add it to your .env file")
        return 1
        
    if not supabase_key:
        print_colored(RED, "Error: Supabase API key not provided")
        print("Please provide a Supabase API key using --key or add it to your .env file")
        return 1
    
    # Extract hostname from URL
    parsed_url = urlparse(supabase_url)
    hostname = parsed_url.netloc
    
    print_colored(BLUE, f"Testing connection to Supabase at {hostname}")
    
    # Test 1: Socket connection to PostgreSQL port
    print_section("PostgreSQL Socket Connection Test")
    db_connected = test_socket_connection(hostname, 5432)
    
    # Test 2: Socket connection to HTTP port
    print_section("HTTP Socket Connection Test")
    http_socket_connected = test_socket_connection(hostname, 443)
    
    # Test 3: HTTP connection to REST API
    print_section("HTTP API Connection Test")
    http_api_connected = test_http_connection(supabase_url, supabase_key)
    
    # Summary
    print_section("Connection Test Summary")
    
    if db_connected:
        print_colored(GREEN, "✓ PostgreSQL port is reachable")
    else:
        print_colored(RED, "✗ PostgreSQL port is not reachable")
    
    if http_socket_connected:
        print_colored(GREEN, "✓ HTTPS port is reachable")
    else:
        print_colored(RED, "✗ HTTPS port is not reachable")
    
    if http_api_connected:
        print_colored(GREEN, "✓ Supabase REST API is accessible")
    else:
        print_colored(RED, "✗ Supabase REST API is not accessible")
    
    # Overall result
    if db_connected and http_socket_connected and http_api_connected:
        print_colored(GREEN, "\n✓ All tests passed! Your Supabase instance is fully accessible.")
        return 0
    elif http_socket_connected or http_api_connected:
        print_colored(YELLOW, "\n⚠ Partial connectivity. Some features may work.")
        return 1
    else:
        print_colored(RED, "\n✗ Connection tests failed. Check your network settings and Supabase configuration.")
        return 2

if __name__ == '__main__':
    sys.exit(main())