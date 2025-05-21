#!/usr/bin/env python3
"""
CryoProtect v2 - API Endpoint Performance Benchmark

This script benchmarks the performance of all API endpoints and generates a comprehensive report.
"""

import os
import sys
import json
import time
import statistics
import requests
import concurrent.futures
import argparse
from datetime import datetime
import psutil
import logging
from typing import Dict, List, Any, Tuple, Optional

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler("api_benchmark.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

# Default configuration
DEFAULT_CONFIG = {
    "base_url": "http://localhost:5000",
    "num_requests": 50,  # Number of requests per endpoint
    "concurrency": 10,   # Number of concurrent requests
    "timeout": 10,       # Request timeout in seconds
    "warmup_requests": 5,  # Number of warmup requests before benchmarking
    "auth_token": None,  # Authentication token (if required)
    "endpoints_file": "memory-bank/api_endpoints.json",  # Path to endpoints file
    "output_file": "API_Performance_Benchmark_Report.md",  # Output report file
    "performance_thresholds": {
        "excellent": 100,   # ms
        "good": 300,        # ms
        "acceptable": 1000  # ms
    }
}

class APIBenchmark:
    """API Benchmark class to test API endpoint performance."""
    
    def __init__(self, config: Dict[str, Any]):
        """Initialize the benchmark with the given configuration."""
        self.config = config
        self.base_url = config["base_url"]
        self.num_requests = config["num_requests"]
        self.concurrency = config["concurrency"]
        self.timeout = config["timeout"]
        self.warmup_requests = config["warmup_requests"]
        self.auth_token = config["auth_token"]
        self.endpoints_file = config["endpoints_file"]
        self.output_file = config["output_file"]
        self.performance_thresholds = config["performance_thresholds"]
        
        # Load endpoints
        self.endpoints = self._load_endpoints()
        
        # Results storage
        self.results = {}
        
    def _load_endpoints(self) -> Dict[str, Any]:
        """Load API endpoints from the endpoints file."""
        try:
            with open(self.endpoints_file, 'r') as f:
                data = json.load(f)
                return data.get("endpoints", {})
        except Exception as e:
            logger.error(f"Error loading endpoints file: {str(e)}")
            sys.exit(1)
    
    def _create_session(self) -> requests.Session:
        """Create a requests session with common headers."""
        session = requests.Session()
        headers = {
            "Content-Type": "application/json",
            "Accept": "application/json"
        }
        
        if self.auth_token:
            headers["Authorization"] = f"Bearer {self.auth_token}"
            
        session.headers.update(headers)
        return session
    
    def _get_endpoint_url(self, endpoint_name: str, endpoint_data: Dict[str, Any]) -> str:
        """Get the full URL for an endpoint."""
        path = endpoint_data.get("path", "")
        return f"{self.base_url}{path}"
    
    def _get_request_data(self, endpoint_name: str, endpoint_data: Dict[str, Any]) -> Dict[str, Any]:
        """Get the request data for an endpoint."""
        method = endpoint_data.get("method", "GET")
        url = self._get_endpoint_url(endpoint_name, endpoint_data)
        
        request_data = {
            "method": method,
            "url": url,
            "timeout": self.timeout
        }
        
        # Add request body for POST/PUT methods
        if method in ["POST", "PUT"] and "example_request" in endpoint_data:
            # Extract JSON body from example request if available
            example = endpoint_data.get("example_request", "")
            if "{" in example and "}" in example:
                json_start = example.find("{")
                json_end = example.rfind("}") + 1
                try:
                    json_body = json.loads(example[json_start:json_end])
                    request_data["json"] = json_body
                except:
                    # If parsing fails, check if there's a request_body field
                    if "request_body" in endpoint_data:
                        # Create a minimal valid request body
                        json_body = {}
                        for field, field_data in endpoint_data["request_body"].items():
                            if field_data.get("required", False):
                                # Use example value or a default based on type
                                field_type = field_data.get("type", "string")
                                if field_type == "string":
                                    json_body[field] = "test"
                                elif field_type == "integer":
                                    json_body[field] = 1
                                elif field_type == "number":
                                    json_body[field] = 1.0
                                elif field_type == "boolean":
                                    json_body[field] = True
                                elif field_type == "array":
                                    json_body[field] = []
                                elif field_type == "object":
                                    json_body[field] = {}
                        
                        if json_body:
                            request_data["json"] = json_body
        
        # Add path parameters if needed
        if "{" in url and "}" in url and "path_parameters" in endpoint_data:
            for param, param_data in endpoint_data["path_parameters"].items():
                placeholder = "{" + param + "}"
                if placeholder in url:
                    # Use a default value based on type
                    param_type = param_data.get("type", "string")
                    if "UUID" in param_type or "uuid" in param_type:
                        # Use a fixed UUID for testing
                        param_value = "00000000-0000-0000-0000-000000000000"
                    elif param_type == "integer":
                        param_value = "1"
                    else:
                        param_value = "test"
                    
                    url = url.replace(placeholder, param_value)
                    request_data["url"] = url
        
        return request_data
    
    def _make_request(self, session: requests.Session, request_data: Dict[str, Any]) -> Tuple[float, int, bool]:
        """Make a request and return response time, status code, and success flag."""
        method = request_data["method"]
        url = request_data["url"]
        
        start_time = time.time()
        success = False
        status_code = 0
        
        try:
            if method == "GET":
                response = session.get(url, timeout=request_data["timeout"])
            elif method == "POST":
                response = session.post(url, json=request_data.get("json"), timeout=request_data["timeout"])
            elif method == "PUT":
                response = session.put(url, json=request_data.get("json"), timeout=request_data["timeout"])
            elif method == "DELETE":
                response = session.delete(url, timeout=request_data["timeout"])
            else:
                logger.warning(f"Unsupported method: {method}")
                return 0, 0, False
            
            status_code = response.status_code
            success = 200 <= status_code < 300
        except requests.exceptions.Timeout:
            logger.warning(f"Request timed out: {url}")
            status_code = 408
        except requests.exceptions.ConnectionError:
            logger.warning(f"Connection error: {url}")
            status_code = 503
        except Exception as e:
            logger.warning(f"Request error: {url} - {str(e)}")
            status_code = 500
        
        end_time = time.time()
        response_time = (end_time - start_time) * 1000  # Convert to milliseconds
        
        return response_time, status_code, success
    
    def _benchmark_endpoint(self, endpoint_name: str, endpoint_data: Dict[str, Any]) -> Dict[str, Any]:
        """Benchmark a single endpoint."""
        logger.info(f"Benchmarking endpoint: {endpoint_name}")
        
        request_data = self._get_request_data(endpoint_name, endpoint_data)
        session = self._create_session()
        
        # Perform warmup requests
        for _ in range(self.warmup_requests):
            self._make_request(session, request_data)
        
        # Perform benchmark requests
        response_times = []
        status_codes = []
        successes = []
        
        start_time = time.time()
        
        with concurrent.futures.ThreadPoolExecutor(max_workers=self.concurrency) as executor:
            futures = [executor.submit(self._make_request, session, request_data) for _ in range(self.num_requests)]
            
            for future in concurrent.futures.as_completed(futures):
                response_time, status_code, success = future.result()
                response_times.append(response_time)
                status_codes.append(status_code)
                successes.append(success)
        
        end_time = time.time()
        total_time = end_time - start_time
        
        # Calculate metrics
        success_count = sum(successes)
        error_count = len(successes) - success_count
        error_rate = error_count / len(successes) if successes else 0
        
        if response_times:
            min_time = min(response_times)
            max_time = max(response_times)
            avg_time = statistics.mean(response_times)
            median_time = statistics.median(response_times)
            p95_time = sorted(response_times)[int(len(response_times) * 0.95)] if response_times else 0
            throughput = len(response_times) / total_time if total_time > 0 else 0
        else:
            min_time = max_time = avg_time = median_time = p95_time = throughput = 0
        
        # Determine performance rating
        if avg_time <= self.performance_thresholds["excellent"]:
            rating = "Excellent"
        elif avg_time <= self.performance_thresholds["good"]:
            rating = "Good"
        elif avg_time <= self.performance_thresholds["acceptable"]:
            rating = "Acceptable"
        else:
            rating = "Poor"
        
        # Collect status code distribution
        status_code_counts = {}
        for code in status_codes:
            status_code_counts[code] = status_code_counts.get(code, 0) + 1
        
        return {
            "endpoint": endpoint_name,
            "method": endpoint_data.get("method", "GET"),
            "path": endpoint_data.get("path", ""),
            "authentication_required": endpoint_data.get("authentication_required", False),
            "min_time": min_time,
            "max_time": max_time,
            "avg_time": avg_time,
            "median_time": median_time,
            "p95_time": p95_time,
            "throughput": throughput,
            "success_count": success_count,
            "error_count": error_count,
            "error_rate": error_rate,
            "status_codes": status_code_counts,
            "rating": rating,
            "total_requests": len(response_times),
            "total_time": total_time
        }
    
    def run_benchmarks(self) -> None:
        """Run benchmarks for all endpoints."""
        logger.info(f"Starting API benchmarks for {len(self.endpoints)} endpoints")
        
        # Track system resource usage
        cpu_usage = []
        memory_usage = []
        
        # Monitor resource usage in a separate thread
        def monitor_resources():
            while self.benchmarking:
                cpu_usage.append(psutil.cpu_percent(interval=1))
                memory_usage.append(psutil.virtual_memory().percent)
        
        import threading
        self.benchmarking = True
        monitor_thread = threading.Thread(target=monitor_resources)
        monitor_thread.daemon = True
        monitor_thread.start()
        
        start_time = time.time()
        
        # Benchmark each endpoint
        for endpoint_name, endpoint_data in self.endpoints.items():
            try:
                result = self._benchmark_endpoint(endpoint_name, endpoint_data)
                self.results[endpoint_name] = result
            except Exception as e:
                logger.error(f"Error benchmarking endpoint {endpoint_name}: {str(e)}")
                self.results[endpoint_name] = {
                    "endpoint": endpoint_name,
                    "error": str(e)
                }
        
        end_time = time.time()
        self.benchmarking = False
        
        # Calculate overall metrics
        self.overall_metrics = {
            "total_endpoints": len(self.endpoints),
            "successful_endpoints": sum(1 for r in self.results.values() if "error" not in r),
            "total_time": end_time - start_time,
            "cpu_usage": {
                "min": min(cpu_usage) if cpu_usage else 0,
                "max": max(cpu_usage) if cpu_usage else 0,
                "avg": statistics.mean(cpu_usage) if cpu_usage else 0
            },
            "memory_usage": {
                "min": min(memory_usage) if memory_usage else 0,
                "max": max(memory_usage) if memory_usage else 0,
                "avg": statistics.mean(memory_usage) if memory_usage else 0
            }
        }
        
        logger.info(f"Completed API benchmarks in {self.overall_metrics['total_time']:.2f} seconds")
    
    def generate_report(self) -> None:
        """Generate a Markdown report with benchmark results."""
        logger.info(f"Generating benchmark report: {self.output_file}")
        
        report_lines = []
        
        # Report header
        report_lines.append("# CryoProtect v2 API Performance Benchmark Report")
        report_lines.append("")
        report_lines.append(f"**Date:** {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        report_lines.append(f"**Base URL:** {self.base_url}")
        report_lines.append(f"**Requests per endpoint:** {self.num_requests}")
        report_lines.append(f"**Concurrency level:** {self.concurrency}")
        report_lines.append("")
        
        # Overall metrics
        report_lines.append("## Overall Metrics")
        report_lines.append("")
        report_lines.append(f"**Total endpoints tested:** {self.overall_metrics['total_endpoints']}")
        report_lines.append(f"**Successfully tested endpoints:** {self.overall_metrics['successful_endpoints']}")
        report_lines.append(f"**Total benchmark time:** {self.overall_metrics['total_time']:.2f} seconds")
        report_lines.append("")
        
        # System resource usage
        report_lines.append("### System Resource Usage")
        report_lines.append("")
        report_lines.append("| Resource | Minimum | Maximum | Average |")
        report_lines.append("|----------|---------|---------|---------|")
        report_lines.append(f"| CPU | {self.overall_metrics['cpu_usage']['min']:.2f}% | {self.overall_metrics['cpu_usage']['max']:.2f}% | {self.overall_metrics['cpu_usage']['avg']:.2f}% |")
        report_lines.append(f"| Memory | {self.overall_metrics['memory_usage']['min']:.2f}% | {self.overall_metrics['memory_usage']['max']:.2f}% | {self.overall_metrics['memory_usage']['avg']:.2f}% |")
        report_lines.append("")
        
        # Performance summary
        report_lines.append("## Performance Summary")
        report_lines.append("")
        
        # Group endpoints by rating
        ratings = {"Excellent": [], "Good": [], "Acceptable": [], "Poor": [], "Error": []}
        
        for endpoint_name, result in self.results.items():
            if "error" in result:
                ratings["Error"].append(endpoint_name)
            else:
                ratings[result["rating"]].append(endpoint_name)
        
        report_lines.append("| Rating | Count | Percentage |")
        report_lines.append("|--------|-------|------------|")
        
        for rating in ["Excellent", "Good", "Acceptable", "Poor", "Error"]:
            count = len(ratings[rating])
            percentage = (count / self.overall_metrics['total_endpoints']) * 100
            report_lines.append(f"| {rating} | {count} | {percentage:.2f}% |")
        
        report_lines.append("")
        
        # Top 5 fastest endpoints
        report_lines.append("### Top 5 Fastest Endpoints")
        report_lines.append("")
        report_lines.append("| Endpoint | Method | Path | Avg Response Time (ms) | Throughput (req/s) |")
        report_lines.append("|----------|--------|------|------------------------|---------------------|")
        
        fastest_endpoints = sorted(
            [r for r in self.results.values() if "error" not in r],
            key=lambda x: x["avg_time"]
        )[:5]
        
        for result in fastest_endpoints:
            report_lines.append(f"| {result['endpoint']} | {result['method']} | {result['path']} | {result['avg_time']:.2f} | {result['throughput']:.2f} |")
        
        report_lines.append("")
        
        # Top 5 slowest endpoints
        report_lines.append("### Top 5 Slowest Endpoints")
        report_lines.append("")
        report_lines.append("| Endpoint | Method | Path | Avg Response Time (ms) | Throughput (req/s) |")
        report_lines.append("|----------|--------|------|------------------------|---------------------|")
        
        slowest_endpoints = sorted(
            [r for r in self.results.values() if "error" not in r],
            key=lambda x: x["avg_time"],
            reverse=True
        )[:5]
        
        for result in slowest_endpoints:
            report_lines.append(f"| {result['endpoint']} | {result['method']} | {result['path']} | {result['avg_time']:.2f} | {result['throughput']:.2f} |")
        
        report_lines.append("")
        
        # Endpoints with highest error rates
        report_lines.append("### Endpoints with Highest Error Rates")
        report_lines.append("")
        report_lines.append("| Endpoint | Method | Path | Error Rate | Common Status Codes |")
        report_lines.append("|----------|--------|------|------------|---------------------|")
        
        error_endpoints = sorted(
            [r for r in self.results.values() if "error" not in r and r["error_rate"] > 0],
            key=lambda x: x["error_rate"],
            reverse=True
        )[:5]
        
        for result in error_endpoints:
            status_codes = ", ".join([f"{code}: {count}" for code, count in result["status_codes"].items() if code >= 400])
            report_lines.append(f"| {result['endpoint']} | {result['method']} | {result['path']} | {result['error_rate']:.2f} | {status_codes} |")
        
        report_lines.append("")
        
        # Detailed results for each endpoint
        report_lines.append("## Detailed Endpoint Results")
        report_lines.append("")
        
        for endpoint_name, result in sorted(self.results.items()):
            report_lines.append(f"### {endpoint_name}")
            report_lines.append("")
            
            if "error" in result:
                report_lines.append(f"**Error:** {result['error']}")
                report_lines.append("")
                continue
            
            report_lines.append(f"**Method:** {result['method']}")
            report_lines.append(f"**Path:** {result['path']}")
            report_lines.append(f"**Authentication Required:** {result['authentication_required']}")
            report_lines.append(f"**Performance Rating:** {result['rating']}")
            report_lines.append("")
            
            report_lines.append("#### Response Time Metrics")
            report_lines.append("")
            report_lines.append("| Metric | Value |")
            report_lines.append("|--------|-------|")
            report_lines.append(f"| Minimum | {result['min_time']:.2f} ms |")
            report_lines.append(f"| Maximum | {result['max_time']:.2f} ms |")
            report_lines.append(f"| Average | {result['avg_time']:.2f} ms |")
            report_lines.append(f"| Median | {result['median_time']:.2f} ms |")
            report_lines.append(f"| 95th Percentile | {result['p95_time']:.2f} ms |")
            report_lines.append("")
            
            report_lines.append("#### Throughput and Error Metrics")
            report_lines.append("")
            report_lines.append("| Metric | Value |")
            report_lines.append("|--------|-------|")
            report_lines.append(f"| Throughput | {result['throughput']:.2f} requests/second |")
            report_lines.append(f"| Success Count | {result['success_count']} |")
            report_lines.append(f"| Error Count | {result['error_count']} |")
            report_lines.append(f"| Error Rate | {result['error_rate']:.2f} |")
            report_lines.append("")
            
            report_lines.append("#### Status Code Distribution")
            report_lines.append("")
            report_lines.append("| Status Code | Count | Percentage |")
            report_lines.append("|-------------|-------|------------|")
            
            for code, count in sorted(result["status_codes"].items()):
                percentage = (count / result["total_requests"]) * 100
                report_lines.append(f"| {code} | {count} | {percentage:.2f}% |")
            
            report_lines.append("")
        
        # Recommendations
        report_lines.append("## Recommendations")
        report_lines.append("")
        
        # Find endpoints that need improvement
        poor_endpoints = [r for r in self.results.values() if "error" not in r and r["rating"] == "Poor"]
        high_error_endpoints = [r for r in self.results.values() if "error" not in r and r["error_rate"] > 0.1]
        
        if poor_endpoints:
            report_lines.append("### Endpoints Requiring Performance Optimization")
            report_lines.append("")
            report_lines.append("The following endpoints have poor performance (average response time > 1000ms) and should be optimized:")
            report_lines.append("")
            
            for result in poor_endpoints:
                report_lines.append(f"1. **{result['endpoint']}** ({result['method']} {result['path']})")
                report_lines.append(f"   - Average response time: {result['avg_time']:.2f} ms")
                report_lines.append(f"   - 95th percentile: {result['p95_time']:.2f} ms")
                report_lines.append(f"   - Possible optimizations:")
                report_lines.append(f"     - Review database queries for optimization opportunities")
                report_lines.append(f"     - Consider adding caching for frequently accessed data")
                report_lines.append(f"     - Check for N+1 query problems")
                report_lines.append("")
        
        if high_error_endpoints:
            report_lines.append("### Endpoints with High Error Rates")
            report_lines.append("")
            report_lines.append("The following endpoints have high error rates (>10%) and should be investigated:")
            report_lines.append("")
            
            for result in high_error_endpoints:
                report_lines.append(f"1. **{result['endpoint']}** ({result['method']} {result['path']})")
                report_lines.append(f"   - Error rate: {result['error_rate']:.2f}")
                report_lines.append(f"   - Common status codes: {', '.join([f'{code}: {count}' for code, count in result['status_codes'].items() if code >= 400])}")
                report_lines.append(f"   - Recommendations:")
                report_lines.append(f"     - Review error handling")
                report_lines.append(f"     - Check input validation")
                report_lines.append(f"     - Verify authentication and authorization logic")
                report_lines.append("")
        
        # General recommendations
        report_lines.append("### General Recommendations")
        report_lines.append("")
        report_lines.append("1. **Implement Caching:** Consider implementing caching for frequently accessed data to reduce database load and improve response times.")
        report_lines.append("2. **Optimize Database Queries:** Review and optimize database queries, especially for endpoints with poor performance.")
        report_lines.append("3. **Connection Pooling:** Ensure database connection pooling is properly configured to handle concurrent requests efficiently.")
        report_lines.append("4. **Rate Limiting:** Implement or adjust rate limiting to prevent abuse and ensure fair resource allocation.")
        report_lines.append("5. **Load Testing:** Conduct regular load testing to identify performance bottlenecks before they impact users.")
        report_lines.append("6. **Monitoring:** Set up continuous monitoring of API performance to detect and address issues proactively.")
        report_lines.append("")
        
        # Conclusion
        report_lines.append("## Conclusion")
        report_lines.append("")
        
        # Calculate overall rating
        excellent_percent = len(ratings["Excellent"]) / self.overall_metrics['total_endpoints'] * 100
        good_percent = len(ratings["Good"]) / self.overall_metrics['total_endpoints'] * 100
        acceptable_percent = len(ratings["Acceptable"]) / self.overall_metrics['total_endpoints'] * 100
        poor_percent = len(ratings["Poor"]) / self.overall_metrics['total_endpoints'] * 100
        error_percent = len(ratings["Error"]) / self.overall_metrics['total_endpoints'] * 100
        
        if excellent_percent + good_percent > 80:
            overall_assessment = "The API performance is generally good, with most endpoints responding quickly and reliably."
        elif excellent_percent + good_percent + acceptable_percent > 80:
            overall_assessment = "The API performance is acceptable, but there is room for improvement in several endpoints."
        else:
            overall_assessment = "The API performance needs significant improvement, with many endpoints showing poor response times or high error rates."
        
        report_lines.append(overall_assessment)
        report_lines.append("")
        
        if poor_endpoints or high_error_endpoints:
            report_lines.append("Priority should be given to addressing the specific endpoints highlighted in the recommendations section.")
        else:
            report_lines.append("While no critical performance issues were identified, continued monitoring and optimization is recommended to maintain and improve performance as the application evolves.")
        
        # Write report to file
        with open(self.output_file, 'w') as f:
            f.write("\n".join(report_lines))
        
        logger.info(f"Benchmark report generated: {self.output_file}")

def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="API Endpoint Performance Benchmark")
    
    parser.add_argument("--base-url", default=DEFAULT_CONFIG["base_url"],
                        help=f"Base URL of the API (default: {DEFAULT_CONFIG['base_url']})")
    
    parser.add_argument("--num-requests", type=int, default=DEFAULT_CONFIG["num_requests"],
                        help=f"Number of requests per endpoint (default: {DEFAULT_CONFIG['num_requests']})")
    
    parser.add_argument("--concurrency", type=int, default=DEFAULT_CONFIG["concurrency"],
                        help=f"Number of concurrent requests (default: {DEFAULT_CONFIG['concurrency']})")
    
    parser.add_argument("--timeout", type=int, default=DEFAULT_CONFIG["timeout"],
                        help=f"Request timeout in seconds (default: {DEFAULT_CONFIG['timeout']})")
    
    parser.add_argument("--warmup-requests", type=int, default=DEFAULT_CONFIG["warmup_requests"],
                        help=f"Number of warmup requests (default: {DEFAULT_CONFIG['warmup_requests']})")
    
    parser.add_argument("--auth-token", default=DEFAULT_CONFIG["auth_token"],
                        help="Authentication token (if required)")
    
    parser.add_argument("--endpoints-file", default=DEFAULT_CONFIG["endpoints_file"],
                        help=f"Path to endpoints file (default: {DEFAULT_CONFIG['endpoints_file']})")
    
    parser.add_argument("--output-file", default=DEFAULT_CONFIG["output_file"],
                        help=f"Output report file (default: {DEFAULT_CONFIG['output_file']})")
    
    return parser.parse_args()

def main():
    """Main function."""
    args = parse_args()
    
    # Create configuration from arguments
    config = {
        "base_url": args.base_url,
        "num_requests": args.num_requests,
        "concurrency": args.concurrency,
        "timeout": args.timeout,
        "warmup_requests": args.warmup_requests,
        "auth_token": args.auth_token,
        "endpoints_file": args.endpoints_file,
        "output_file": args.output_file,
        "performance_thresholds": DEFAULT_CONFIG["performance_thresholds"]
    }
    
    # Create and run benchmark
    benchmark = APIBenchmark(config)
    benchmark.run_benchmarks()
    benchmark.generate_report()
    
    return 0

if __name__ == "__main__":
    sys.exit(main())