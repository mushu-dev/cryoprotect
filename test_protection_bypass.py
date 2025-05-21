#!/usr/bin/env python
"""
Test script for the Vercel protection bypass functionality.
This script tests API access with and without the bypass token.
"""

import requests
import argparse
import sys

# Default token from our configuration
DEFAULT_TOKEN = "TAt23KbtFE8dkZobJU3hpgTP4L5ja07V"

def test_api_access(url, token=None, use_header=True, verbose=False):
    """
    Test API access with and without the protection bypass.
    
    Args:
        url: The base URL of the Vercel deployment
        token: The protection bypass token (optional)
        use_header: Whether to use the token in header (vs query param)
        verbose: Whether to print verbose output
    
    Returns:
        bool: True if tests pass, False otherwise
    """
    # Ensure URL doesn't end with a slash
    if url.endswith('/'):
        url = url[:-1]
    
    # First let's check if the base URL is accessible at all
    try:
        test_response = requests.get(url, timeout=10)
        if verbose:
            print(f"Base URL access test: {test_response.status_code}")
            print(f"Base URL response: {test_response.text[:200]}...")
        
        if test_response.status_code >= 400:
            print(f"❌ WARNING: Base URL returned status {test_response.status_code}")
    except Exception as e:
        print(f"❌ WARNING: Cannot access base URL: {str(e)}")
    
    # Build our request parameters
    headers = {'Accept': 'application/json'}
    params = {}
    
    if token:
        if use_header:
            headers['x-protection-bypass'] = token
            method = "header"
        else:
            params['bypass'] = token
            method = "query parameter"
        
        print(f"Testing with bypass token via {method}...")
    else:
        print("Testing without bypass token...")
    
    # Test regular API endpoint (should be protected)
    api_url = f"{url}/api/v1/health"
    try:
        if verbose:
            print(f"Requesting {api_url}")
            print(f"Headers: {headers}")
            print(f"Params: {params}")
        
        response = requests.get(api_url, headers=headers, params=params, timeout=10)
        
        if verbose:
            print(f"Response status: {response.status_code}")
            print(f"Response body: {response.text[:200]}...")
        
        if token:
            # With token, should succeed
            if response.status_code == 200:
                print(f"✅ SUCCESS: Protected API endpoint accessible with token")
                return True
            else:
                print(f"❌ FAILURE: Protected API endpoint returned {response.status_code}")
                return False
        else:
            # Without token, should fail
            if response.status_code == 403:
                print(f"✅ SUCCESS: Protected API endpoint correctly returned 403 without token")
                return True
            else:
                print(f"❌ FAILURE: Protected API endpoint did not block access (returned {response.status_code})")
                return False
    
    except Exception as e:
        print(f"❌ ERROR: {str(e)}")
        return False

def main():
    parser = argparse.ArgumentParser(description='Test Vercel protection bypass functionality')
    parser.add_argument('url', help='Base URL of the Vercel deployment to test')
    parser.add_argument('-t', '--token', default=DEFAULT_TOKEN, help='Protection bypass token')
    parser.add_argument('-q', '--query', action='store_true', help='Use query parameter instead of header')
    parser.add_argument('-v', '--verbose', action='store_true', help='Verbose output')
    args = parser.parse_args()
    
    print(f"Testing protection bypass for: {args.url}")
    
    # Test without token (should fail)
    no_token_test = test_api_access(args.url, token=None, verbose=args.verbose)
    
    # Test with token (should succeed)
    with_token_test = test_api_access(
        args.url, 
        token=args.token, 
        use_header=not args.query,
        verbose=args.verbose
    )
    
    # Overall result
    if no_token_test and with_token_test:
        print("\n✅ OVERALL: Protection bypass is working correctly!")
        return 0
    else:
        print("\n❌ OVERALL: Protection bypass tests failed")
        return 1

if __name__ == "__main__":
    sys.exit(main())