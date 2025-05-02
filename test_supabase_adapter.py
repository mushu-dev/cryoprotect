#!/usr/bin/env python3
"""
Test script for the SupabaseAdapter class.

This script demonstrates how to use the SupabaseAdapter to normalize
different types of Supabase responses.
"""

import json
from supabase_adapter import SupabaseAdapter


def main():
    """Run tests for the SupabaseAdapter class."""
    print("Testing SupabaseAdapter with various response formats...\n")
    
    # Create adapter instance
    adapter = SupabaseAdapter()
    
    # Test case 1: Standard successful response with data
    print("Test Case 1: Standard successful response with data")
    response1 = MockResponse(data=[{"id": 1, "name": "Test"}], error=None)
    result1 = adapter.normalize_response(response1)
    print_result(result1)
    
    # Test case 2: Error response with message
    print("\nTest Case 2: Error response with message")
    response2 = MockResponse(data=None, error={"message": "Not found"})
    result2 = adapter.normalize_response(response2)
    print_result(result2)
    
    # Test case 3: Dictionary with result/status format
    print("\nTest Case 3: Dictionary with result/status format")
    response3 = {"result": [{"id": 2, "name": "Another Test"}], "status": "ok"}
    result3 = adapter.normalize_response(response3)
    print_result(result3)
    
    # Test case 4: Dictionary with error message
    print("\nTest Case 4: Dictionary with error message")
    response4 = {"error": {"detail": "Permission denied", "code": 403}}
    result4 = adapter.normalize_response(response4)
    print_result(result4)
    
    # Test case 5: Simple list response
    print("\nTest Case 5: Simple list response")
    response5 = [{"id": 3, "name": "Item 1"}, {"id": 4, "name": "Item 2"}]
    result5 = adapter.normalize_response(response5)
    print_result(result5)
    
    # Test case 6: Empty response
    print("\nTest Case 6: Empty response")
    response6 = None
    result6 = adapter.normalize_response(response6)
    print_result(result6)


def print_result(result):
    """Print the normalized result in a readable format."""
    print(f"Success: {result['success']}")
    print(f"Data: {json.dumps(result['data'], indent=2) if result['data'] else None}")
    print(f"Error: {result['error']}")


class MockResponse:
    """Mock Supabase response object with data and error attributes."""
    
    def __init__(self, data=None, error=None):
        self.data = data
        self.error = error


if __name__ == "__main__":
    main()