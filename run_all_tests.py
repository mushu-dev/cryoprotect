#!/usr/bin/env python3
"""
CryoProtect v2 - Test Runner

This script runs all the test scripts and generates a comprehensive test report.
"""

import os
import sys
import json
import logging
import subprocess
from datetime import datetime
from typing import Dict, List, Any, Optional, Tuple

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler("test_runner.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

class TestRunner:
    """Test runner class."""

    def __init__(self):
        """Initialize the test runner."""
        self.test_scripts = [
            "test_database_schema.py",
            "test_rls_policies.py",
            "test_api_integration.py",
            "test_rdkit_integration.py",
            "test_end_to_end_workflows.py"
        ]
        self.test_results = {
            "status": "Not Started",
            "total_tests": 0,
            "passed_tests": 0,
            "failed_tests": 0,
            "skipped_tests": 0,
            "component_results": []
        }

    def run_tests(self) -> Dict[str, Any]:
        """Run all test scripts."""
        self.test_results["status"] = "Running"
        self.test_results["start_time"] = datetime.now().isoformat()

        # Run each test script
        for script in self.test_scripts:
            self.run_test_script(script)

        # Calculate overall results
        self.test_results["total_tests"] = sum(cr["total_tests"] for cr in self.test_results["component_results"])
        self.test_results["passed_tests"] = sum(cr["passed_tests"] for cr in self.test_results["component_results"])
        self.test_results["failed_tests"] = sum(cr["failed_tests"] for cr in self.test_results["component_results"])
        self.test_results["skipped_tests"] = sum(cr["skipped_tests"] for cr in self.test_results["component_results"])

        if self.test_results["failed_tests"] == 0:
            self.test_results["status"] = "Passed"
        else:
            self.test_results["status"] = "Failed"

        self.test_results["end_time"] = datetime.now().isoformat()

        return self.test_results

    def run_test_script(self, script: str) -> None:
        """Run a test script and collect the results."""
        logger.info(f"Running {script}...")

        try:
            # Run the test script
            result = subprocess.run(
                [sys.executable, script],
                capture_output=True,
                text=True
            )

            # Check if the script ran successfully
            if result.returncode == 0:
                logger.info(f"{script} passed")
                status = "Passed"
            else:
                logger.error(f"{script} failed with exit code {result.returncode}")
                logger.error(f"Error output: {result.stderr}")
                status = "Failed"

            # Try to load the JSON results file
            json_filename = script.replace(".py", "_results.json")
            if os.path.exists(json_filename):
                with open(json_filename, 'r') as f:
                    component_results = json.load(f)
            else:
                # Create a placeholder result if the JSON file doesn't exist
                component_results = {
                    "status": status,
                    "total_tests": 0,
                    "passed_tests": 0,
                    "failed_tests": 1 if status == "Failed" else 0,
                    "skipped_tests": 0,
                    "test_cases": []
                }

            # Add the component results to the overall results
            self.test_results["component_results"].append({
                "script": script,
                "status": status,
                "total_tests": component_results.get("total_tests", 0),
                "passed_tests": component_results.get("passed_tests", 0),
                "failed_tests": component_results.get("failed_tests", 0),
                "skipped_tests": component_results.get("skipped_tests", 0),
                "stdout": result.stdout,
                "stderr": result.stderr
            })

        except Exception as e:
            logger.error(f"Error running {script}: {str(e)}")
            self.test_results["component_results"].append({
                "script": script,
                "status": "Error",
                "total_tests": 0,
                "passed_tests": 0,
                "failed_tests": 1,
                "skipped_tests": 0,
                "error": str(e)
            })

    def save_results(self, filename: str) -> None:
        """Save the test results to a file."""
        try:
            with open(filename, 'w') as f:
                json.dump(self.test_results, f, indent=2)
            logger.info(f"Test results saved to {filename}")
        except Exception as e:
            logger.error(f"Failed to save test results: {str(e)}")

    def generate_report(self, filename: str) -> None:
        """Generate a Markdown report of the test results."""
        try:
            with open(filename, 'w') as f:
                f.write("# CryoProtect v2 Test Report\n\n")
                f.write(f"**Date:** {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")

                f.write("## Summary\n\n")
                f.write(f"- **Total Tests:** {self.test_results['total_tests']}\n")
                f.write(f"- **Passed Tests:** {self.test_results['passed_tests']}\n")
                f.write(f"- **Failed Tests:** {self.test_results['failed_tests']}\n")
                f.write(f"- **Skipped Tests:** {self.test_results['skipped_tests']}\n")
                f.write(f"- **Overall Status:** {self.test_results['status']}\n\n")

                f.write("## Component Results\n\n")
                f.write("| Component | Status | Total Tests | Passed | Failed | Skipped |\n")
                f.write("|-----------|--------|-------------|--------|--------|--------|\n")

                for component in self.test_results["component_results"]:
                    f.write(f"| {component['script']} | {component['status']} | {component['total_tests']} | {component['passed_tests']} | {component['failed_tests']} | {component['skipped_tests']} |\n")

                f.write("\n## Database Schema and Data Integrity\n\n")
                db_component = next((c for c in self.test_results["component_results"] if c["script"] == "test_database_schema.py"), None)
                if db_component:
                    f.write(f"**Status:** {db_component['status']}\n\n")
                    f.write(f"**Tests Passed:** {db_component['passed_tests']}/{db_component['total_tests']}\n\n")
                    
                    if db_component["status"] == "Failed":
                        f.write("### Issues Found\n\n")
                        f.write("The following issues were found in the database schema and data integrity:\n\n")
                        
                        # Try to load the detailed results
                        try:
                            with open("database_schema_test_results.json", 'r') as f:
                                detailed_results = json.load(f)
                                
                                for test_case in detailed_results.get("test_cases", []):
                                    if test_case["status"] == "Failed":
                                        f.write(f"- **{test_case['name']}**: {test_case['message']}\n")
                        except Exception:
                            f.write("- Detailed results not available\n")
                    else:
                        f.write("No issues found in the database schema and data integrity.\n")

                f.write("\n## RLS Policy Effectiveness\n\n")
                rls_component = next((c for c in self.test_results["component_results"] if c["script"] == "test_rls_policies.py"), None)
                if rls_component:
                    f.write(f"**Status:** {rls_component['status']}\n\n")
                    f.write(f"**Tests Passed:** {rls_component['passed_tests']}/{rls_component['total_tests']}\n\n")
                    
                    if rls_component["status"] == "Failed":
                        f.write("### Issues Found\n\n")
                        f.write("The following issues were found in the RLS policies:\n\n")
                        
                        # Try to load the detailed results
                        try:
                            with open("rls_policy_test_results.json", 'r') as f:
                                detailed_results = json.load(f)
                                
                                for test_case in detailed_results.get("test_cases", []):
                                    if test_case["status"] == "Failed":
                                        f.write(f"- **{test_case['name']}**: {test_case['message']}\n")
                        except Exception:
                            f.write("- Detailed results not available\n")
                    else:
                        f.write("No issues found in the RLS policies.\n")

                f.write("\n## API Endpoint Functionality\n\n")
                api_component = next((c for c in self.test_results["component_results"] if c["script"] == "test_api_integration.py"), None)
                if api_component:
                    f.write(f"**Status:** {api_component['status']}\n\n")
                    f.write(f"**Tests Passed:** {api_component['passed_tests']}/{api_component['total_tests']}\n\n")
                    
                    if api_component["status"] == "Failed":
                        f.write("### Issues Found\n\n")
                        f.write("The following issues were found in the API endpoints:\n\n")
                        
                        # Try to load the detailed results
                        try:
                            with open("api_integration_test_results.json", 'r') as f:
                                detailed_results = json.load(f)
                                
                                for test_case in detailed_results.get("test_cases", []):
                                    if test_case["status"] == "Failed":
                                        f.write(f"- **{test_case['name']}**: {test_case['message']}\n")
                        except Exception:
                            f.write("- Detailed results not available\n")
                    else:
                        f.write("No issues found in the API endpoints.\n")

                f.write("\n## RDKit Integration\n\n")
                rdkit_component = next((c for c in self.test_results["component_results"] if c["script"] == "test_rdkit_integration.py"), None)
                if rdkit_component:
                    f.write(f"**Status:** {rdkit_component['status']}\n\n")
                    f.write(f"**Tests Passed:** {rdkit_component['passed_tests']}/{rdkit_component['total_tests']}\n\n")
                    
                    if rdkit_component["status"] == "Failed":
                        f.write("### Issues Found\n\n")
                        f.write("The following issues were found in the RDKit integration:\n\n")
                        
                        # Try to load the detailed results
                        try:
                            with open("rdkit_integration_test_results.json", 'r') as f:
                                detailed_results = json.load(f)
                                
                                for test_case in detailed_results.get("test_cases", []):
                                    if test_case["status"] == "Failed":
                                        f.write(f"- **{test_case['name']}**: {test_case['message']}\n")
                        except Exception:
                            f.write("- Detailed results not available\n")
                    else:
                        f.write("No issues found in the RDKit integration.\n")

                f.write("\n## End-to-End Workflows\n\n")
                e2e_component = next((c for c in self.test_results["component_results"] if c["script"] == "test_end_to_end_workflows.py"), None)
                if e2e_component:
                    f.write(f"**Status:** {e2e_component['status']}\n\n")
                    f.write(f"**Tests Passed:** {e2e_component['passed_tests']}/{e2e_component['total_tests']}\n\n")
                    
                    if e2e_component["status"] == "Failed":
                        f.write("### Issues Found\n\n")
                        f.write("The following issues were found in the end-to-end workflows:\n\n")
                        
                        # Try to load the detailed results
                        try:
                            with open("end_to_end_test_results.json", 'r') as f:
                                detailed_results = json.load(f)
                                
                                for test_case in detailed_results.get("test_cases", []):
                                    if test_case["status"] == "Failed":
                                        f.write(f"- **{test_case['name']}**: {test_case['message']}\n")
                        except Exception:
                            f.write("- Detailed results not available\n")
                    else:
                        f.write("No issues found in the end-to-end workflows.\n")

                f.write("\n## Recommendations\n\n")
                
                if self.test_results["failed_tests"] > 0:
                    f.write("Based on the test results, the following recommendations are made:\n\n")
                    
                    # Database recommendations
                    db_component = next((c for c in self.test_results["component_results"] if c["script"] == "test_database_schema.py"), None)
                    if db_component and db_component["status"] == "Failed":
                        f.write("### Database Recommendations\n\n")
                        f.write("1. Fix the database schema issues identified in the test results\n")
                        f.write("2. Ensure all foreign key relationships are properly defined\n")
                        f.write("3. Check for and fix any data inconsistencies\n\n")
                    
                    # RLS recommendations
                    rls_component = next((c for c in self.test_results["component_results"] if c["script"] == "test_rls_policies.py"), None)
                    if rls_component and rls_component["status"] == "Failed":
                        f.write("### RLS Policy Recommendations\n\n")
                        f.write("1. Ensure RLS is enabled on all tables\n")
                        f.write("2. Fix any issues with user access control\n")
                        f.write("3. Verify that service role can bypass RLS as expected\n\n")
                    
                    # API recommendations
                    api_component = next((c for c in self.test_results["component_results"] if c["script"] == "test_api_integration.py"), None)
                    if api_component and api_component["status"] == "Failed":
                        f.write("### API Recommendations\n\n")
                        f.write("1. Fix any non-functional API endpoints\n")
                        f.write("2. Ensure proper error handling in all endpoints\n")
                        f.write("3. Verify that authentication and authorization are working correctly\n\n")
                    
                    # RDKit recommendations
                    rdkit_component = next((c for c in self.test_results["component_results"] if c["script"] == "test_rdkit_integration.py"), None)
                    if rdkit_component and rdkit_component["status"] == "Failed":
                        f.write("### RDKit Recommendations\n\n")
                        f.write("1. Ensure RDKit is properly installed and configured\n")
                        f.write("2. Fix any issues with molecular property calculations\n")
                        f.write("3. Verify that visualization and substructure searching work correctly\n\n")
                    
                    # End-to-end recommendations
                    e2e_component = next((c for c in self.test_results["component_results"] if c["script"] == "test_end_to_end_workflows.py"), None)
                    if e2e_component and e2e_component["status"] == "Failed":
                        f.write("### Workflow Recommendations\n\n")
                        f.write("1. Fix any issues with the molecule creation and analysis workflow\n")
                        f.write("2. Ensure the mixture creation and analysis workflow works correctly\n")
                        f.write("3. Verify that prediction and experiment comparison works as expected\n\n")
                else:
                    f.write("All tests passed! No immediate recommendations are necessary.\n\n")
                    f.write("For future improvements, consider:\n\n")
                    f.write("1. Expanding test coverage to include more edge cases\n")
                    f.write("2. Implementing automated testing as part of CI/CD pipeline\n")
                    f.write("3. Adding performance testing to ensure the system scales well\n")

                f.write("\n## Validation Status\n\n")
                
                if self.test_results["status"] == "Passed":
                    f.write("✅ **VALIDATED**: The system meets all requirements specified in the project documentation.\n")
                else:
                    f.write("❌ **NOT VALIDATED**: The system does not meet all requirements specified in the project documentation. Please address the issues identified in this report.\n")

            logger.info(f"Report generated at {filename}")
        except Exception as e:
            logger.error(f"Failed to generate report: {str(e)}")

def main():
    """Main function."""
    logger.info("Starting test runner")
    
    # Create and run the tests
    runner = TestRunner()
    results = runner.run_tests()
    
    # Save the results
    runner.save_results("test_runner_results.json")
    
    # Generate a report
    report_filename = f"TEST_RESULTS_REPORT_{datetime.now().strftime('%Y%m%d_%H%M%S')}.md"
    runner.generate_report(report_filename)
    
    # Print summary
    logger.info(f"Test Status: {results['status']}")
    logger.info(f"Total Tests: {results['total_tests']}")
    logger.info(f"Passed Tests: {results['passed_tests']}")
    logger.info(f"Failed Tests: {results['failed_tests']}")
    logger.info(f"Skipped Tests: {results['skipped_tests']}")
    
    # Exit with appropriate status code
    if results["status"] == "Passed":
        logger.info("All tests passed!")
        sys.exit(0)
    else:
        logger.error("Some tests failed. See log for details.")
        sys.exit(1)

if __name__ == "__main__":
    main()