"""
CryoProtect Analyzer - Test Report Generator

This module contains utilities for generating comprehensive test reports.
It collects test results from various test modules and generates a structured report.
"""

import os
import sys
import json
import time
import datetime
import unittest
from unittest.mock import patch, MagicMock

# Add the parent directory to the path so we can import the api package
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

class TestReportGenerator(unittest.TestCase):
    """Test case for the test report generator."""

    def setUp(self):
        """Set up test data for each test."""
        # Create the reports directory if it doesn't exist
        self.reports_dir = os.path.join(os.path.dirname(__file__), 'reports')
        os.makedirs(self.reports_dir, exist_ok=True)
        
        # Sample test results
        self.sample_test_results = {
            "api_endpoint_tests": {
                "success": True,
                "details": "Ran 25 tests, 0 failures, 0 errors",
                "results": {
                    "test_get_molecules": "PASS",
                    "test_get_molecule": "PASS",
                    "test_import_from_pubchem": "PASS",
                    "test_search_molecules_by_name": "PASS",
                    "test_get_mixtures": "PASS",
                    "test_get_mixture": "PASS",
                    "test_create_mixture": "PASS",
                    "test_update_mixture": "PASS",
                    "test_calculate_mixture_score": "PASS",
                    "test_get_predictions": "PASS",
                    "test_add_prediction": "PASS",
                    "test_get_experiments": "PASS",
                    "test_record_experiment": "PASS",
                    "test_compare_prediction_with_experiment": "PASS",
                    "test_calculate_molecular_properties": "PASS",
                    "test_generate_molecule_visualization": "PASS",
                    "test_substructure_search": "PASS",
                    "test_calculate_similarity": "PASS",
                    "test_not_found": "PASS",
                    "test_molecule_not_found": "PASS",
                    "test_bad_request": "PASS",
                    "test_validation_error": "PASS"
                }
            },
            "database_schema_validation": {
                "success": True,
                "details": "Ran 6 tests, 0 failures, 0 errors",
                "results": {
                    "test_table_names_are_plural": "PASS",
                    "test_foreign_key_constraints": "PASS",
                    "test_rls_enabled_on_all_tables": "PASS",
                    "test_rls_policies_for_different_roles": "PASS",
                    "test_views_exist": "PASS",
                    "test_functions_exist": "PASS"
                }
            },
            "authentication_tests": {
                "success": True,
                "details": "Ran 8 tests, 0 failures, 0 errors",
                "results": {
                    "test_user_registration": "PASS",
                    "test_user_login": "PASS",
                    "test_user_logout": "PASS",
                    "test_password_reset_request": "PASS",
                    "test_access_control_anonymous_user": "PASS",
                    "test_access_control_authenticated_user": "PASS",
                    "test_access_control_admin_user": "PASS",
                    "test_user_cannot_access_admin_endpoints": "PASS",
                    "test_user_can_only_modify_own_data": "PASS"
                }
            },
            "performance_tests": {
                "success": True,
                "details": "Ran 3 tests, 0 failures, 0 errors",
                "results": {
                    "test_read_operations_performance": "PASS",
                    "test_write_operations_performance": "PASS",
                    "test_mixed_operations_performance": "PASS"
                },
                "metrics": {
                    "execution_time": 45.23,
                    "avg_response_time": 0.125,
                    "max_response_time": 0.532,
                    "throughput": 78.45,
                    "max_cpu_usage": 65.2,
                    "max_memory_usage": 42.8
                }
            }
        }

    def test_generate_json_report(self):
        """Test generating a JSON report."""
        # Generate the report
        report = self.generate_json_report(self.sample_test_results)
        
        # Verify the report structure
        self.assertIn("timestamp", report)
        self.assertIn("overall_status", report)
        self.assertIn("test_results", report)
        self.assertIn("performance_metrics", report)
        self.assertIn("issues", report)
        self.assertIn("recommendations", report)
        
        # Verify the overall status
        self.assertEqual(report["overall_status"], "SUCCESS")
        
        # Save the report to a file
        report_path = os.path.join(self.reports_dir, f'test_report_{int(time.time())}.json')
        with open(report_path, 'w') as f:
            json.dump(report, f, indent=2)
        
        print(f"JSON report generated: {report_path}")

    def test_generate_html_report(self):
        """Test generating an HTML report."""
        # Generate the report
        html_content = self.generate_html_report(self.sample_test_results)
        
        # Verify the report content
        self.assertIn("<html>", html_content)
        self.assertIn("<title>CryoProtect Test Report</title>", html_content)
        self.assertIn("Overall Status:", html_content)
        self.assertIn("API Endpoint Tests", html_content)
        self.assertIn("Database Schema Validation", html_content)
        self.assertIn("Authentication Tests", html_content)
        self.assertIn("Performance Tests", html_content)
        self.assertIn("Performance Metrics", html_content)
        self.assertIn("Issues", html_content)
        self.assertIn("Recommendations", html_content)
        
        # Save the report to a file
        report_path = os.path.join(self.reports_dir, f'test_report_{int(time.time())}.html')
        with open(report_path, 'w') as f:
            f.write(html_content)
        
        print(f"HTML report generated: {report_path}")

    def test_generate_text_report(self):
        """Test generating a text report."""
        # Generate the report
        text_content = self.generate_text_report(self.sample_test_results)
        
        # Verify the report content
        self.assertIn("CryoProtect Analyzer - Test Report", text_content)
        self.assertIn("Overall Status:", text_content)
        self.assertIn("API Endpoint Tests", text_content)
        self.assertIn("Database Schema Validation", text_content)
        self.assertIn("Authentication Tests", text_content)
        self.assertIn("Performance Tests", text_content)
        self.assertIn("Performance Metrics", text_content)
        self.assertIn("Issues", text_content)
        self.assertIn("Recommendations", text_content)
        
        # Save the report to a file
        report_path = os.path.join(self.reports_dir, f'test_report_{int(time.time())}.txt')
        with open(report_path, 'w') as f:
            f.write(text_content)
        
        print(f"Text report generated: {report_path}")

    def generate_json_report(self, test_results):
        """Generate a JSON report from test results."""
        # Determine overall status
        overall_status = "SUCCESS" if all(result["success"] for result in test_results.values()) else "ERROR"
        
        # Extract performance metrics
        performance_metrics = {}
        if "performance_tests" in test_results and "metrics" in test_results["performance_tests"]:
            performance_metrics = test_results["performance_tests"]["metrics"]
        
        # Identify issues
        issues = []
        for test_type, result in test_results.items():
            if not result["success"]:
                issues.append(f"{test_type} failed: {result['details']}")
        
        # Generate recommendations
        recommendations = []
        if not overall_status == "SUCCESS":
            recommendations.append("Fix failing tests before proceeding to production")
        
        if "performance_tests" in test_results and "metrics" in test_results["performance_tests"]:
            metrics = test_results["performance_tests"]["metrics"]
            if metrics.get("avg_response_time", 0) > 0.5:
                recommendations.append("Optimize database queries to improve response time")
            if metrics.get("max_cpu_usage", 0) > 80:
                recommendations.append("Consider scaling up database resources to handle the load")
            if metrics.get("throughput", 0) < 50:
                recommendations.append("Implement caching to improve throughput")
        
        # Create the report
        report = {
            "timestamp": datetime.datetime.now().isoformat(),
            "overall_status": overall_status,
            "test_results": test_results,
            "performance_metrics": performance_metrics,
            "issues": issues,
            "recommendations": recommendations
        }
        
        return report

    def generate_html_report(self, test_results):
        """Generate an HTML report from test results."""
        # Determine overall status
        overall_status = "SUCCESS" if all(result["success"] for result in test_results.values()) else "ERROR"
        status_color = "green" if overall_status == "SUCCESS" else "red"
        
        # Extract performance metrics
        performance_metrics = {}
        if "performance_tests" in test_results and "metrics" in test_results["performance_tests"]:
            performance_metrics = test_results["performance_tests"]["metrics"]
        
        # Identify issues
        issues = []
        for test_type, result in test_results.items():
            if not result["success"]:
                issues.append(f"{test_type} failed: {result['details']}")
        
        # Generate recommendations
        recommendations = []
        if not overall_status == "SUCCESS":
            recommendations.append("Fix failing tests before proceeding to production")
        
        if "performance_tests" in test_results and "metrics" in test_results["performance_tests"]:
            metrics = test_results["performance_tests"]["metrics"]
            if metrics.get("avg_response_time", 0) > 0.5:
                recommendations.append("Optimize database queries to improve response time")
            if metrics.get("max_cpu_usage", 0) > 80:
                recommendations.append("Consider scaling up database resources to handle the load")
            if metrics.get("throughput", 0) < 50:
                recommendations.append("Implement caching to improve throughput")
        
        # Create the HTML content
        html_content = f"""
        <html>
        <head>
            <title>CryoProtect Test Report</title>
            <style>
                body {{ font-family: Arial, sans-serif; margin: 20px; }}
                h1 {{ color: #333; }}
                h2 {{ color: #666; margin-top: 30px; }}
                .status {{ font-weight: bold; color: {status_color}; }}
                table {{ border-collapse: collapse; width: 100%; margin-top: 10px; }}
                th, td {{ border: 1px solid #ddd; padding: 8px; text-align: left; }}
                th {{ background-color: #f2f2f2; }}
                tr:nth-child(even) {{ background-color: #f9f9f9; }}
                .pass {{ color: green; }}
                .fail {{ color: red; }}
                .metrics {{ margin-top: 20px; }}
                .issues {{ margin-top: 20px; color: red; }}
                .recommendations {{ margin-top: 20px; }}
            </style>
        </head>
        <body>
            <h1>CryoProtect Analyzer - Test Report</h1>
            <p>Generated: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
            
            <h2>Overall Status: <span class="status">{overall_status}</span></h2>
        """
        
        # Add test results
        for test_type, result in test_results.items():
            status = "PASSED" if result["success"] else "FAILED"
            status_class = "pass" if result["success"] else "fail"
            
            html_content += f"""
            <h2>{test_type.replace('_', ' ').title()}</h2>
            <p>Status: <span class="{status_class}">{status}</span></p>
            <p>Details: {result["details"]}</p>
            """
            
            if "results" in result:
                html_content += """
                <table>
                    <tr>
                        <th>Test</th>
                        <th>Result</th>
                    </tr>
                """
                
                for test_name, test_result in result["results"].items():
                    result_class = "pass" if test_result == "PASS" else "fail"
                    html_content += f"""
                    <tr>
                        <td>{test_name}</td>
                        <td class="{result_class}">{test_result}</td>
                    </tr>
                    """
                
                html_content += "</table>"
        
        # Add performance metrics
        if performance_metrics:
            html_content += """
            <h2>Performance Metrics</h2>
            <table class="metrics">
                <tr>
                    <th>Metric</th>
                    <th>Value</th>
                </tr>
            """
            
            for metric, value in performance_metrics.items():
                html_content += f"""
                <tr>
                    <td>{metric.replace('_', ' ').title()}</td>
                    <td>{value}</td>
                </tr>
                """
            
            html_content += "</table>"
        
        # Add issues
        if issues:
            html_content += """
            <h2>Issues</h2>
            <ul class="issues">
            """
            
            for issue in issues:
                html_content += f"<li>{issue}</li>"
            
            html_content += "</ul>"
        
        # Add recommendations
        if recommendations:
            html_content += """
            <h2>Recommendations</h2>
            <ul class="recommendations">
            """
            
            for recommendation in recommendations:
                html_content += f"<li>{recommendation}</li>"
            
            html_content += "</ul>"
        
        html_content += """
        </body>
        </html>
        """
        
        return html_content

    def generate_text_report(self, test_results):
        """Generate a text report from test results."""
        # Determine overall status
        overall_status = "SUCCESS" if all(result["success"] for result in test_results.values()) else "ERROR"
        
        # Extract performance metrics
        performance_metrics = {}
        if "performance_tests" in test_results and "metrics" in test_results["performance_tests"]:
            performance_metrics = test_results["performance_tests"]["metrics"]
        
        # Identify issues
        issues = []
        for test_type, result in test_results.items():
            if not result["success"]:
                issues.append(f"{test_type} failed: {result['details']}")
        
        # Generate recommendations
        recommendations = []
        if not overall_status == "SUCCESS":
            recommendations.append("Fix failing tests before proceeding to production")
        
        if "performance_tests" in test_results and "metrics" in test_results["performance_tests"]:
            metrics = test_results["performance_tests"]["metrics"]
            if metrics.get("avg_response_time", 0) > 0.5:
                recommendations.append("Optimize database queries to improve response time")
            if metrics.get("max_cpu_usage", 0) > 80:
                recommendations.append("Consider scaling up database resources to handle the load")
            if metrics.get("throughput", 0) < 50:
                recommendations.append("Implement caching to improve throughput")
        
        # Create the text content
        text_content = f"""
CryoProtect Analyzer - Test Report
Generated: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}

Overall Status: {overall_status}

"""
        
        # Add test results
        for test_type, result in test_results.items():
            status = "PASSED" if result["success"] else "FAILED"
            
            text_content += f"""
{test_type.replace('_', ' ').title()}
{'-' * len(test_type)}
Status: {status}
Details: {result["details"]}
"""
            
            if "results" in result:
                text_content += "\nTest Results:\n"
                
                for test_name, test_result in result["results"].items():
                    text_content += f"  {test_name}: {test_result}\n"
        
        # Add performance metrics
        if performance_metrics:
            text_content += """
Performance Metrics
------------------
"""
            
            for metric, value in performance_metrics.items():
                text_content += f"  {metric.replace('_', ' ').title()}: {value}\n"
        
        # Add issues
        if issues:
            text_content += """
Issues
------
"""
            
            for issue in issues:
                text_content += f"  - {issue}\n"
        
        # Add recommendations
        if recommendations:
            text_content += """
Recommendations
--------------
"""
            
            for recommendation in recommendations:
                text_content += f"  - {recommendation}\n"
        
        return text_content

def generate_report(test_results, format="json"):
    """Generate a test report in the specified format."""
    generator = TestReportGenerator()
    
    if format == "json":
        return generator.generate_json_report(test_results)
    elif format == "html":
        return generator.generate_html_report(test_results)
    elif format == "text":
        return generator.generate_text_report(test_results)
    else:
        raise ValueError(f"Unsupported report format: {format}")

if __name__ == '__main__':
    unittest.main()