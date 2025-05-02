#!/usr/bin/env python3
"""
CryoProtect v2 - Run All Database Performance Tests

This script runs all database performance tests and generates a comprehensive report.
"""

import os
import sys
import json
import time
import subprocess
from datetime import datetime

# Set up logging
import logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler("performance_testing.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

def run_test_script(script_name, description):
    """Run a test script and return its exit code."""
    logger.info(f"Running {description}...")
    
    start_time = time.time()
    result = subprocess.run([sys.executable, script_name], capture_output=True, text=True)
    end_time = time.time()
    
    if result.returncode == 0:
        logger.info(f"{description} completed successfully in {end_time - start_time:.2f} seconds")
    else:
        logger.error(f"{description} failed with exit code {result.returncode}")
        logger.error(f"STDOUT: {result.stdout}")
        logger.error(f"STDERR: {result.stderr}")
    
    return result.returncode == 0

def load_json_report(filename):
    """Load a JSON report file."""
    try:
        with open(filename, 'r') as f:
            return json.load(f)
    except Exception as e:
        logger.error(f"Error loading {filename}: {str(e)}")
        return None

def load_text_report(filename):
    """Load a text report file."""
    try:
        with open(filename, 'r') as f:
            return f.read()
    except Exception as e:
        logger.error(f"Error loading {filename}: {str(e)}")
        return None

def generate_comprehensive_report(performance_report, query_plan_report, concurrent_load_report):
    """Generate a comprehensive report from all test results."""
    report = {
        "timestamp": datetime.now().isoformat(),
        "performance_test": performance_report,
        "query_plan_analysis": query_plan_report,
        "concurrent_load_test": concurrent_load_report
    }
    
    # Save the comprehensive report
    with open("comprehensive_performance_report.json", "w") as f:
        json.dump(report, f, indent=2)
    
    # Generate a text summary
    summary = []
    
    summary.append("=" * 80)
    summary.append("CryoProtect v2 Comprehensive Database Performance Report")
    summary.append("=" * 80)
    summary.append(f"Timestamp: {report['timestamp']}")
    summary.append("")
    
    # Performance test summary
    if performance_report:
        summary.append("-" * 80)
        summary.append("Database Performance Test Summary")
        summary.append("-" * 80)
        
        # Resource usage
        if "metrics" in performance_report and "resource_usage" in performance_report["metrics"]:
            cpu = performance_report["metrics"]["resource_usage"]["cpu"]
            memory = performance_report["metrics"]["resource_usage"]["memory"]
            
            summary.append(f"CPU Usage: Min={cpu['min']:.2f}%, Max={cpu['max']:.2f}%, Avg={cpu['avg']:.2f}%")
            summary.append(f"Memory Usage: Min={memory['min']:.2f}%, Max={memory['max']:.2f}%, Avg={memory['avg']:.2f}%")
            summary.append("")
        
        # Operation performance
        if "metrics" in performance_report and "operations" in performance_report["metrics"]:
            summary.append("Top 5 Slowest Operations:")
            
            # Sort operations by p95 response time
            sorted_ops = sorted(
                performance_report["metrics"]["operations"].items(),
                key=lambda x: x[1]["p95_ms"],
                reverse=True
            )
            
            for i, (op_name, op_stats) in enumerate(sorted_ops[:5]):
                summary.append(f"{i+1}. {op_name}: {op_stats['p95_ms']:.2f} ms (P95)")
            
            summary.append("")
        
        # Bottlenecks
        if "analysis" in performance_report and "bottlenecks" in performance_report["analysis"]:
            bottlenecks = performance_report["analysis"]["bottlenecks"]
            if bottlenecks:
                summary.append(f"Identified {len(bottlenecks)} potential bottlenecks")
                for bottleneck in bottlenecks:
                    summary.append(f"- {bottleneck['operation']}: {bottleneck['type']}")
                summary.append("")
        
        # Recommendations
        if "recommendations" in performance_report:
            summary.append("Key Recommendations:")
            for i, rec in enumerate(performance_report["recommendations"][:5]):
                summary.append(f"{i+1}. {rec['recommendation']}")
            
            if len(performance_report["recommendations"]) > 5:
                summary.append(f"   ... and {len(performance_report['recommendations']) - 5} more recommendations")
            
            summary.append("")
    
    # Query plan analysis summary
    if query_plan_report:
        summary.append("-" * 80)
        summary.append("Query Plan Analysis Summary")
        summary.append("-" * 80)
        
        # Slow queries
        slow_queries = []
        if "queries" in query_plan_report:
            for query in query_plan_report["queries"]:
                if query["execution_time_ms"] and query["execution_time_ms"] > 100:
                    slow_queries.append({
                        "name": query["name"],
                        "time_ms": query["execution_time_ms"]
                    })
        
        if slow_queries:
            summary.append("Slow Queries (>100ms):")
            for i, query in enumerate(sorted(slow_queries, key=lambda x: x["time_ms"], reverse=True)):
                summary.append(f"{i+1}. {query['name']}: {query['time_ms']:.2f} ms")
            summary.append("")
        
        # Collect all unique suggestions
        all_suggestions = []
        if "queries" in query_plan_report:
            for query in query_plan_report["queries"]:
                if "suggestions" in query:
                    for suggestion in query["suggestions"]:
                        if suggestion not in all_suggestions:
                            all_suggestions.append(suggestion)
        
        if all_suggestions:
            summary.append(f"Found {len(all_suggestions)} optimization suggestions:")
            for i, suggestion in enumerate(all_suggestions[:5]):
                summary.append(f"{i+1}. {suggestion['description']}")
            
            if len(all_suggestions) > 5:
                summary.append(f"   ... and {len(all_suggestions) - 5} more suggestions")
            
            summary.append("")
    
    # Concurrent load test summary
    if concurrent_load_report:
        summary.append("-" * 80)
        summary.append("Concurrent Load Test Summary")
        summary.append("-" * 80)
        
        if "summary" in concurrent_load_report:
            test_summary = concurrent_load_report["summary"]
            
            summary.append(f"Test Duration: {test_summary.get('test_duration', 'N/A'):.2f} seconds")
            summary.append(f"Total Operations: {test_summary.get('total_operations', 'N/A')}")
            summary.append(f"Average Response Time: {test_summary.get('avg_response_time_ms', 'N/A'):.2f} ms")
            summary.append(f"Overall Throughput: {test_summary.get('throughput_per_second', 'N/A'):.2f} operations/second")
            summary.append("")
            
            # Operation performance
            if "operations" in test_summary:
                summary.append("Operation Performance:")
                
                # Sort operations by throughput
                sorted_ops = sorted(
                    test_summary["operations"].items(),
                    key=lambda x: x[1]["throughput_per_second"],
                    reverse=True
                )
                
                for op_name, op_stats in sorted_ops[:5]:
                    summary.append(f"- {op_name}:")
                    summary.append(f"  Count: {op_stats['count']} (Success: {op_stats['success_count']}, Error: {op_stats['error_count']})")
                    summary.append(f"  Avg: {op_stats['avg_ms']:.2f} ms, P95: {op_stats['p95_ms']:.2f} ms")
                    summary.append(f"  Throughput: {op_stats['throughput_per_second']:.2f} ops/sec")
                
                if len(test_summary["operations"]) > 5:
                    summary.append(f"   ... and {len(test_summary['operations']) - 5} more operations")
                
                summary.append("")
            
            # Errors
            if "errors" in test_summary and test_summary["errors"]:
                total_errors = sum(len(errors) for errors in test_summary["errors"].values())
                summary.append(f"Total Errors: {total_errors}")
                
                for op_name, errors in list(test_summary["errors"].items())[:3]:
                    summary.append(f"- {op_name}: {len(errors)} errors")
                
                if len(test_summary["errors"]) > 3:
                    summary.append(f"   ... and errors in {len(test_summary['errors']) - 3} more operations")
                
                summary.append("")
    
    # Overall assessment
    summary.append("=" * 80)
    summary.append("Overall Assessment")
    summary.append("=" * 80)
    
    # Determine if the database is ready for production
    is_ready = True
    critical_issues = []
    
    # Check for critical performance issues
    if performance_report and "analysis" in performance_report:
        # Check for high CPU or memory usage
        if "resource_issues" in performance_report["analysis"]:
            for issue in performance_report["analysis"]["resource_issues"]:
                if issue["type"] in ["high_cpu_usage", "high_memory_usage"] and issue["value"] > 90:
                    is_ready = False
                    critical_issues.append(f"Critical {issue['type']} ({issue['value']:.2f}%)")
        
        # Check for very slow operations
        if "slow_operations" in performance_report["analysis"]:
            for op in performance_report["analysis"]["slow_operations"]:
                if op["p95_ms"] > 2000:  # Operations taking more than 2 seconds
                    is_ready = False
                    critical_issues.append(f"Critical slow operation: {op['operation']} ({op['p95_ms']:.2f} ms)")
    
    # Check for critical issues in concurrent load test
    if concurrent_load_report and "summary" in concurrent_load_report:
        test_summary = concurrent_load_report["summary"]
        
        # Check for high error rates
        if "operations" in test_summary:
            for op_name, op_stats in test_summary["operations"].items():
                if op_stats["error_count"] > op_stats["success_count"]:
                    is_ready = False
                    critical_issues.append(f"High error rate in {op_name} ({op_stats['error_count']} errors, {op_stats['success_count']} successes)")
        
        # Check for very slow average response time
        if "avg_response_time_ms" in test_summary and test_summary["avg_response_time_ms"] > 1000:
            is_ready = False
            critical_issues.append(f"High average response time ({test_summary['avg_response_time_ms']:.2f} ms)")
    
    # Final assessment
    if is_ready:
        summary.append("ASSESSMENT: The database appears to be READY for production use.")
        summary.append("The performance tests indicate that the database can handle the expected production load.")
        summary.append("However, consider implementing the optimization recommendations to further improve performance.")
    else:
        summary.append("ASSESSMENT: The database is NOT READY for production use.")
        summary.append("The following critical issues need to be addressed:")
        for issue in critical_issues:
            summary.append(f"- {issue}")
        summary.append("Review the detailed reports and implement the recommended optimizations before deploying to production.")
    
    summary.append("")
    summary.append("=" * 80)
    
    # Save the summary report
    with open("comprehensive_performance_report.txt", "w") as f:
        f.write("\n".join(summary))
    
    return "\n".join(summary)

def main():
    """Main function."""
    print("\n" + "=" * 80)
    print("CryoProtect v2 - Running All Database Performance Tests")
    print("=" * 80 + "\n")
    
    # Track test results
    test_results = {
        "performance_test": False,
        "query_plan_analysis": False,
        "concurrent_load_test": False
    }
    
    # Run database performance test
    test_results["performance_test"] = run_test_script(
        "test_database_performance.py",
        "Database Performance Test"
    )
    
    # Run query plan analysis
    test_results["query_plan_analysis"] = run_test_script(
        "analyze_query_plans.py",
        "Query Plan Analysis"
    )
    
    # Run concurrent load test
    test_results["concurrent_load_test"] = run_test_script(
        "simulate_concurrent_load.py",
        "Concurrent Load Test"
    )
    
    # Load test reports
    performance_report = load_json_report("database_performance_report.json") if test_results["performance_test"] else None
    query_plan_report = load_json_report("query_plan_analysis_report.json") if test_results["query_plan_analysis"] else None
    concurrent_load_report = load_json_report("concurrent_load_test_report.json") if test_results["concurrent_load_test"] else None
    
    # Generate comprehensive report
    summary = generate_comprehensive_report(performance_report, query_plan_report, concurrent_load_report)
    
    print("\n" + "=" * 80)
    print("CryoProtect v2 - All Database Performance Tests Completed")
    print("=" * 80)
    print(f"Performance Test: {'SUCCESS' if test_results['performance_test'] else 'FAILED'}")
    print(f"Query Plan Analysis: {'SUCCESS' if test_results['query_plan_analysis'] else 'FAILED'}")
    print(f"Concurrent Load Test: {'SUCCESS' if test_results['concurrent_load_test'] else 'FAILED'}")
    print("=" * 80)
    print(f"Comprehensive report saved to comprehensive_performance_report.txt and comprehensive_performance_report.json")
    print("=" * 80 + "\n")
    
    # Print the summary report
    print(summary)
    
    return 0

if __name__ == "__main__":
    exit(main())