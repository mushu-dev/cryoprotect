#!/usr/bin/env python
"""
CryoProtect Analyzer - Supabase API Integration Test Runner

This script runs the Supabase API integration tests and generates a detailed report.
"""

import os
import sys
import unittest
import json
import datetime
from unittest import TextTestRunner, TestSuite
from unittest.loader import TestLoader

# Add the parent directory to the path so we can import the api package
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

# Import the test case
from tests.test_supabase_api_integration import TestSupabaseAPIIntegration

def run_tests():
    """Run the Supabase API integration tests and generate a report."""
    # Create test suite
    loader = TestLoader()
    suite = TestSuite()
    
    # Add tests to the suite
    suite.addTest(loader.loadTestsFromTestCase(TestSupabaseAPIIntegration))
    
    # Create a test runner
    runner = TextTestRunner(verbosity=2)
    
    # Run the tests
    print("\n===== Running Supabase API Integration Tests =====\n")
    result = runner.run(suite)
    
    # Generate report
    generate_report(result)
    
    return result.wasSuccessful()

def generate_report(result):
    """Generate a detailed test report."""
    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    report_filename = f"SUPABASE_API_TEST_REPORT_{timestamp}.md"
    report_path = os.path.join("reports", report_filename)
    
    # Create reports directory if it doesn't exist
    os.makedirs("reports", exist_ok=True)
    
    # Collect test results
    test_results = {
        "total": result.testsRun,
        "passed": result.testsRun - len(result.failures) - len(result.errors),
        "failed": len(result.failures),
        "errors": len(result.errors),
        "failures": [
            {
                "test": str(test),
                "message": err
            } for test, err in result.failures
        ],
        "error_details": [
            {
                "test": str(test),
                "message": err
            } for test, err in result.errors
        ]
    }
    
    # Calculate success rate
    success_rate = (test_results["passed"] / test_results["total"]) * 100 if test_results["total"] > 0 else 0
    
    # Write report
    with open(report_path, "w") as f:
        f.write("# Supabase API Integration Test Report\n\n")
        f.write(f"**Date:** {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        
        f.write("## Summary\n\n")
        f.write(f"- **Total Tests:** {test_results['total']}\n")
        f.write(f"- **Passed:** {test_results['passed']}\n")
        f.write(f"- **Failed:** {test_results['failed']}\n")
        f.write(f"- **Errors:** {test_results['errors']}\n")
        f.write(f"- **Success Rate:** {success_rate:.2f}%\n\n")
        
        if test_results["failures"]:
            f.write("## Failures\n\n")
            for i, failure in enumerate(test_results["failures"], 1):
                f.write(f"### {i}. {failure['test']}\n\n")
                f.write("```\n")
                f.write(failure["message"])
                f.write("\n```\n\n")
        
        if test_results["error_details"]:
            f.write("## Errors\n\n")
            for i, error in enumerate(test_results["error_details"], 1):
                f.write(f"### {i}. {error['test']}\n\n")
                f.write("```\n")
                f.write(error["message"])
                f.write("\n```\n\n")
        
        f.write("## Authentication Test Results\n\n")
        f.write("- Unauthenticated access to public endpoints: " + 
                ("✅ PASSED" if "test_unauthenticated_get_molecules" not in [f["test"] for f in test_results["failures"]] else "❌ FAILED") + "\n")
        f.write("- Unauthenticated access restriction: " + 
                ("✅ PASSED" if "test_unauthenticated_cannot_create_mixture" not in [f["test"] for f in test_results["failures"]] else "❌ FAILED") + "\n")
        f.write("- Authenticated user access to own data: " + 
                ("✅ PASSED" if "test_authenticated_user_update_own_mixture" not in [f["test"] for f in test_results["failures"]] else "❌ FAILED") + "\n")
        f.write("- Authenticated user access restriction to others' data: " + 
                ("✅ PASSED" if "test_authenticated_user_cannot_update_others_mixture" not in [f["test"] for f in test_results["failures"]] else "❌ FAILED") + "\n")
        f.write("- Admin user access to any data: " + 
                ("✅ PASSED" if "test_admin_user_can_update_any_mixture" not in [f["test"] for f in test_results["failures"]] else "❌ FAILED") + "\n")
        f.write("- Service role access to any data: " + 
                ("✅ PASSED" if "test_service_role_can_access_any_data" not in [f["test"] for f in test_results["failures"]] else "❌ FAILED") + "\n\n")
        
        f.write("## RLS Policy Enforcement\n\n")
        f.write("Based on the test results, the Row Level Security (RLS) policies are:\n\n")
        
        if ("test_authenticated_user_cannot_update_others_mixture" not in [f["test"] for f in test_results["failures"]] and
            "test_admin_user_can_update_any_mixture" not in [f["test"] for f in test_results["failures"]]):
            f.write("✅ **Correctly Enforced**\n\n")
            f.write("The RLS policies are properly restricting access based on user roles:\n")
            f.write("- Regular users can only access their own data\n")
            f.write("- Admin users can access any data\n")
            f.write("- Service roles can bypass RLS policies\n")
        else:
            f.write("❌ **Incorrectly Enforced**\n\n")
            f.write("There are issues with the RLS policy enforcement:\n")
            if "test_authenticated_user_cannot_update_others_mixture" in [f["test"] for f in test_results["failures"]]:
                f.write("- Regular users can access other users' data\n")
            if "test_admin_user_can_update_any_mixture" in [f["test"] for f in test_results["failures"]]:
                f.write("- Admin users cannot access all data\n")
        
        f.write("\n## Recommendations\n\n")
        if test_results["failures"] or test_results["error_details"]:
            f.write("Based on the test results, the following issues should be addressed:\n\n")
            
            if any("test_unauthenticated" in f["test"] for f in test_results["failures"]):
                f.write("1. **Public Access Issues**: Review and fix public access to endpoints\n")
            
            if any("test_authenticated_user_cannot" in f["test"] for f in test_results["failures"]):
                f.write("2. **RLS Policy Issues**: Strengthen RLS policies to prevent unauthorized access\n")
            
            if any("test_admin_user" in f["test"] for f in test_results["failures"]):
                f.write("3. **Admin Access Issues**: Fix admin role permissions\n")
            
            if any("test_service_role" in f["test"] for f in test_results["failures"]):
                f.write("4. **Service Role Issues**: Review service role implementation\n")
            
            if any("test_validation_error" in f["test"] for f in test_results["failures"]):
                f.write("5. **Validation Issues**: Improve input validation\n")
        else:
            f.write("All tests passed successfully. The API endpoints with Supabase integration are working correctly.\n")
            f.write("The authentication system and RLS policies are properly enforced.\n")
    
    print(f"\nTest report generated: {report_path}")
    
    # Save JSON results
    json_report_path = os.path.join("reports", f"supabase_api_test_results_{timestamp}.json")
    with open(json_report_path, "w") as f:
        json.dump(test_results, f, indent=2)
    
    print(f"JSON results saved: {json_report_path}")

if __name__ == "__main__":
    success = run_tests()
    sys.exit(0 if success else 1)