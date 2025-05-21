#!/usr/bin/env python3
"""
Stress testing tool for the Convex connection pool.

This script performs stress testing on the Convex connection pool, simulating
various load patterns and concurrency levels to verify performance and resilience.
"""

import os
import sys
import time
import random
import logging
import argparse
import threading
import json
import concurrent.futures
from typing import Dict, Any, List, Tuple, Optional

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)

logger = logging.getLogger(__name__)

# Add the project directory to the Python path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from database.convex_pool import get_convex_pool, execute_convex_query, convex_connection

# Test configuration defaults
DEFAULT_CONCURRENCY = 10
DEFAULT_REQUESTS = 1000
DEFAULT_DURATION = 60  # seconds
DEFAULT_PATTERN = "steady"  # steady, ramp, spike, oscillate
DEFAULT_ERROR_RATE = 0  # percentage of simulated errors

class StressTestResults:
    """Class to track and summarize stress test results."""
    
    def __init__(self):
        """Initialize results tracking."""
        self.total_requests = 0
        self.successful_requests = 0
        self.failed_requests = 0
        self.total_time = 0
        self.min_time = float('inf')
        self.max_time = 0
        self.start_time = time.time()
        self.end_time = None
        self.status_codes = {}
        self.errors = {}
        self.lock = threading.RLock()
    
    def record_request(self, duration: float, success: bool, status_code: Optional[int] = None, error: Optional[str] = None):
        """
        Record a request result.
        
        Args:
            duration: Request duration in seconds
            success: Whether the request was successful
            status_code: HTTP status code (for successful requests)
            error: Error message (for failed requests)
        """
        with self.lock:
            self.total_requests += 1
            self.total_time += duration
            
            if success:
                self.successful_requests += 1
                if status_code is not None:
                    self.status_codes[status_code] = self.status_codes.get(status_code, 0) + 1
            else:
                self.failed_requests += 1
                if error:
                    error_type = error.split(":")[0] if ":" in error else error
                    self.errors[error_type] = self.errors.get(error_type, 0) + 1
            
            # Update min/max times
            self.min_time = min(self.min_time, duration)
            self.max_time = max(self.max_time, duration)
    
    def complete(self):
        """Mark the test as complete and record the end time."""
        self.end_time = time.time()
    
    def generate_report(self) -> Dict[str, Any]:
        """
        Generate a report of the test results.
        
        Returns:
            Dict with test results and statistics
        """
        with self.lock:
            # Calculate statistics
            total_time = self.end_time - self.start_time if self.end_time else time.time() - self.start_time
            avg_time = self.total_time / max(1, self.total_requests)
            requests_per_second = self.total_requests / max(1, total_time)
            success_rate = (self.successful_requests / max(1, self.total_requests)) * 100
            
            # Generate report
            return {
                "summary": {
                    "total_requests": self.total_requests,
                    "successful_requests": self.successful_requests,
                    "failed_requests": self.failed_requests,
                    "duration_seconds": round(total_time, 2),
                    "requests_per_second": round(requests_per_second, 2),
                    "success_rate_percent": round(success_rate, 2)
                },
                "timing": {
                    "average_request_time": round(avg_time, 3),
                    "min_request_time": round(self.min_time, 3) if self.min_time != float('inf') else 0,
                    "max_request_time": round(self.max_time, 3)
                },
                "status_codes": self.status_codes,
                "errors": self.errors
            }


def run_worker(
    worker_id: int,
    target_url: str,
    num_requests: int,
    results: StressTestResults,
    simulate_errors: bool = False,
    error_rate: float = 0.0,
    delay_range: Tuple[float, float] = (0.0, 0.1)
):
    """
    Worker function to generate load on the connection pool.
    
    Args:
        worker_id: Worker identifier
        target_url: Convex URL
        num_requests: Number of requests to make
        results: Results tracking object
        simulate_errors: Whether to simulate errors
        error_rate: Percentage of requests that should fail (0-100)
        delay_range: Range of delays between requests (min, max)
    """
    logger.info(f"Worker {worker_id} starting with {num_requests} requests")
    
    # Get pool
    pool = get_convex_pool()
    
    for i in range(num_requests):
        start_time = time.time()
        
        try:
            # Simulate artificial errors if configured
            if simulate_errors and random.random() < error_rate / 100:
                error_type = random.choice(["timeout", "connection", "server"])
                
                if error_type == "timeout":
                    # Simulate timeout
                    raise TimeoutError("Simulated timeout error")
                elif error_type == "connection":
                    # Simulate connection error
                    raise ConnectionError("Simulated connection error")
                else:
                    # Simulate server error
                    raise Exception("Simulated server error")
            
            # Make an actual request to Convex
            with pool.connection() as session:
                # Use a simple health check path or similar
                response = session.get(f"{target_url}/api", timeout=5)
                
                # Record result
                duration = time.time() - start_time
                results.record_request(
                    duration=duration,
                    success=True,
                    status_code=response.status_code
                )
            
            # Add a small delay between requests
            if delay_range[1] > 0:
                time.sleep(random.uniform(delay_range[0], delay_range[1]))
        
        except Exception as e:
            # Record failure
            duration = time.time() - start_time
            results.record_request(
                duration=duration,
                success=False,
                error=str(e)
            )
            
            logger.debug(f"Worker {worker_id} request {i+1}/{num_requests} failed: {str(e)}")
    
    logger.info(f"Worker {worker_id} completed {num_requests} requests")


def steady_load_test(
    url: str,
    concurrency: int,
    num_requests: int,
    duration: int,
    error_rate: float = 0.0
) -> StressTestResults:
    """
    Run a steady load test with consistent concurrency.
    
    Args:
        url: Convex URL
        concurrency: Number of concurrent workers
        num_requests: Total number of requests to make
        duration: Maximum test duration in seconds
        error_rate: Percentage of simulated errors
        
    Returns:
        StressTestResults: Test results
    """
    logger.info(f"Starting steady load test with {concurrency} workers and {num_requests} requests")
    
    # Initialize results
    results = StressTestResults()
    
    # Calculate requests per worker
    requests_per_worker = num_requests // concurrency
    remainder = num_requests % concurrency
    
    # Create threadpool
    with concurrent.futures.ThreadPoolExecutor(max_workers=concurrency) as executor:
        futures = []
        
        # Submit worker tasks
        for i in range(concurrency):
            # Add an extra request to some workers to distribute the remainder
            worker_requests = requests_per_worker + (1 if i < remainder else 0)
            
            future = executor.submit(
                run_worker,
                i + 1,
                url,
                worker_requests,
                results,
                True,
                error_rate,
                (0.01, 0.05)  # Small delay between requests
            )
            futures.append(future)
        
        # Wait for completion or timeout
        deadline = time.time() + duration
        for future in concurrent.futures.as_completed(futures):
            if time.time() > deadline:
                logger.warning("Test duration exceeded, stopping remaining workers")
                # Cancel any pending futures
                for f in futures:
                    if not f.done():
                        f.cancel()
                break
            
            # Check for exceptions
            try:
                future.result()
            except Exception as e:
                logger.error(f"Worker failed with error: {str(e)}")
    
    # Mark test as complete
    results.complete()
    
    return results


def ramp_load_test(
    url: str,
    max_concurrency: int,
    num_requests: int,
    duration: int,
    error_rate: float = 0.0
) -> StressTestResults:
    """
    Run a ramp load test with gradually increasing concurrency.
    
    Args:
        url: Convex URL
        max_concurrency: Maximum number of concurrent workers
        num_requests: Total number of requests to make
        duration: Maximum test duration in seconds
        error_rate: Percentage of simulated errors
        
    Returns:
        StressTestResults: Test results
    """
    logger.info(f"Starting ramp load test up to {max_concurrency} workers with {num_requests} requests")
    
    # Initialize results
    results = StressTestResults()
    
    # Calculate concurrency steps
    steps = 5
    step_size = max(1, max_concurrency // steps)
    step_duration = duration // steps
    
    # Calculate requests per step
    requests_per_step = num_requests // steps
    
    start_time = time.time()
    step_concurrency = step_size
    
    # Create threadpool with max concurrency
    with concurrent.futures.ThreadPoolExecutor(max_workers=max_concurrency) as executor:
        # Run each step
        for step in range(steps):
            if time.time() - start_time >= duration:
                logger.warning("Test duration exceeded, stopping test")
                break
            
            # Adjust concurrency for this step
            step_concurrency = min(step_size * (step + 1), max_concurrency)
            logger.info(f"Step {step+1}/{steps}: Running with {step_concurrency} workers")
            
            # Calculate requests per worker for this step
            requests_per_worker = requests_per_step // step_concurrency
            remainder = requests_per_step % step_concurrency
            
            # Submit tasks for this step
            futures = []
            for i in range(step_concurrency):
                worker_requests = requests_per_worker + (1 if i < remainder else 0)
                
                future = executor.submit(
                    run_worker,
                    (step * max_concurrency) + i + 1,
                    url,
                    worker_requests,
                    results,
                    True,
                    error_rate,
                    (0.01, 0.05)
                )
                futures.append(future)
            
            # Wait for step completion or timeout
            step_deadline = time.time() + step_duration
            for future in concurrent.futures.as_completed(futures):
                if time.time() > step_deadline:
                    logger.warning(f"Step {step+1} duration exceeded, moving to next step")
                    break
                
                # Check for exceptions
                try:
                    future.result()
                except Exception as e:
                    logger.error(f"Worker failed with error: {str(e)}")
    
    # Mark test as complete
    results.complete()
    
    return results


def spike_load_test(
    url: str,
    base_concurrency: int,
    peak_concurrency: int,
    num_requests: int,
    duration: int,
    error_rate: float = 0.0
) -> StressTestResults:
    """
    Run a spike load test with periods of high concurrency.
    
    Args:
        url: Convex URL
        base_concurrency: Baseline number of concurrent workers
        peak_concurrency: Peak number of concurrent workers during spikes
        num_requests: Total number of requests to make
        duration: Maximum test duration in seconds
        error_rate: Percentage of simulated errors
        
    Returns:
        StressTestResults: Test results
    """
    logger.info(f"Starting spike load test with baseline {base_concurrency} and peak {peak_concurrency} workers")
    
    # Initialize results
    results = StressTestResults()
    
    # Calculate spike parameters
    num_spikes = 3
    spike_duration = duration // (num_spikes * 3)  # Each cycle is baseline-spike-baseline
    
    # Calculate requests per phase
    total_phases = num_spikes * 2  # alternating baseline and spike phases
    requests_per_phase = num_requests // total_phases
    
    start_time = time.time()
    
    # Create threadpool with max concurrency
    with concurrent.futures.ThreadPoolExecutor(max_workers=peak_concurrency) as executor:
        phase = 0
        
        # Run baseline and spike phases
        while phase < total_phases and time.time() - start_time < duration:
            # Determine if this is a baseline or spike phase
            is_spike = phase % 2 == 1
            phase_concurrency = peak_concurrency if is_spike else base_concurrency
            phase_name = "Spike" if is_spike else "Baseline"
            
            logger.info(f"Phase {phase+1}/{total_phases}: {phase_name} with {phase_concurrency} workers")
            
            # Calculate requests per worker
            requests_per_worker = requests_per_phase // phase_concurrency
            remainder = requests_per_phase % phase_concurrency
            
            # Submit tasks for this phase
            futures = []
            for i in range(phase_concurrency):
                worker_requests = requests_per_worker + (1 if i < remainder else 0)
                
                future = executor.submit(
                    run_worker,
                    (phase * peak_concurrency) + i + 1,
                    url,
                    worker_requests,
                    results,
                    True,
                    error_rate,
                    (0.01, 0.05)
                )
                futures.append(future)
            
            # Wait for phase completion or timeout
            phase_deadline = time.time() + (spike_duration if is_spike else spike_duration * 2)
            for future in concurrent.futures.as_completed(futures):
                if time.time() > phase_deadline:
                    logger.warning(f"Phase {phase+1} duration exceeded, moving to next phase")
                    break
                
                # Check for exceptions
                try:
                    future.result()
                except Exception as e:
                    logger.error(f"Worker failed with error: {str(e)}")
            
            # Move to next phase
            phase += 1
    
    # Mark test as complete
    results.complete()
    
    return results


def oscillating_load_test(
    url: str,
    min_concurrency: int,
    max_concurrency: int,
    num_requests: int,
    duration: int,
    error_rate: float = 0.0
) -> StressTestResults:
    """
    Run an oscillating load test with gradually changing concurrency.
    
    Args:
        url: Convex URL
        min_concurrency: Minimum number of concurrent workers
        max_concurrency: Maximum number of concurrent workers
        num_requests: Total number of requests to make
        duration: Maximum test duration in seconds
        error_rate: Percentage of simulated errors
        
    Returns:
        StressTestResults: Test results
    """
    logger.info(f"Starting oscillating load test between {min_concurrency} and {max_concurrency} workers")
    
    # Initialize results
    results = StressTestResults()
    
    # Calculate oscillation parameters
    num_oscillations = 3
    steps_per_oscillation = 8
    total_steps = num_oscillations * steps_per_oscillation
    step_duration = duration // total_steps
    
    # Calculate requests per step
    requests_per_step = num_requests // total_steps
    
    start_time = time.time()
    
    # Create threadpool with max concurrency
    with concurrent.futures.ThreadPoolExecutor(max_workers=max_concurrency) as executor:
        # Run each step
        for step in range(total_steps):
            if time.time() - start_time >= duration:
                logger.warning("Test duration exceeded, stopping test")
                break
            
            # Calculate concurrency for this step using sine wave pattern
            progress = (step % steps_per_oscillation) / steps_per_oscillation
            oscillation_factor = (math.sin(progress * 2 * math.pi) + 1) / 2  # 0 to 1
            step_concurrency = min_concurrency + int(oscillation_factor * (max_concurrency - min_concurrency))
            
            logger.info(f"Step {step+1}/{total_steps}: Running with {step_concurrency} workers")
            
            # Calculate requests per worker
            requests_per_worker = requests_per_step // step_concurrency
            remainder = requests_per_step % step_concurrency
            
            # Submit tasks for this step
            futures = []
            for i in range(step_concurrency):
                worker_requests = requests_per_worker + (1 if i < remainder else 0)
                
                future = executor.submit(
                    run_worker,
                    (step * max_concurrency) + i + 1,
                    url,
                    worker_requests,
                    results,
                    True,
                    error_rate,
                    (0.01, 0.05)
                )
                futures.append(future)
            
            # Wait for step completion or timeout
            step_deadline = time.time() + step_duration
            for future in concurrent.futures.as_completed(futures):
                if time.time() > step_deadline:
                    logger.warning(f"Step {step+1} duration exceeded, moving to next step")
                    break
                
                # Check for exceptions
                try:
                    future.result()
                except Exception as e:
                    logger.error(f"Worker failed with error: {str(e)}")
    
    # Mark test as complete
    results.complete()
    
    return results


def generate_pool_config_report(pool_config):
    """
    Generate a report about the pool configuration.
    
    Args:
        pool_config: Pool configuration dictionary
        
    Returns:
        Dict: Pool configuration report
    """
    return {
        "convex_url": pool_config.get("url", os.environ.get("CONVEX_URL", "unknown")),
        "min_connections": pool_config.get("min_connections", 3),
        "max_connections": pool_config.get("max_connections", 20),
        "connection_timeout": pool_config.get("connection_timeout", 30),
        "connection_lifetime": pool_config.get("connection_lifetime", 3600),
        "idle_timeout": pool_config.get("idle_timeout", 300),
        "retry_attempts": pool_config.get("retry_attempts", 3),
        "health_check_interval": pool_config.get("health_check_interval", 60),
        "circuit_breaker": {
            "failure_threshold": pool_config.get("circuit_breaker_threshold", 5),
            "recovery_timeout": pool_config.get("circuit_breaker_timeout", 30),
            "recovery_success": pool_config.get("circuit_breaker_success", 2)
        }
    }


def main():
    """Main entry point for the stress tester."""
    parser = argparse.ArgumentParser(description="Convex Connection Pool Stress Tester")
    
    parser.add_argument(
        "--url",
        help="Convex URL (defaults to CONVEX_URL environment variable)",
        default=os.environ.get("CONVEX_URL", "")
    )
    
    parser.add_argument(
        "--concurrency",
        type=int,
        help=f"Number of concurrent workers (default: {DEFAULT_CONCURRENCY})",
        default=DEFAULT_CONCURRENCY
    )
    
    parser.add_argument(
        "--requests",
        type=int,
        help=f"Total number of requests to make (default: {DEFAULT_REQUESTS})",
        default=DEFAULT_REQUESTS
    )
    
    parser.add_argument(
        "--duration",
        type=int,
        help=f"Maximum test duration in seconds (default: {DEFAULT_DURATION})",
        default=DEFAULT_DURATION
    )
    
    parser.add_argument(
        "--pattern",
        choices=["steady", "ramp", "spike", "oscillate"],
        help=f"Load pattern (default: {DEFAULT_PATTERN})",
        default=DEFAULT_PATTERN
    )
    
    parser.add_argument(
        "--error-rate",
        type=float,
        help=f"Percentage of simulated errors (default: {DEFAULT_ERROR_RATE})",
        default=DEFAULT_ERROR_RATE
    )
    
    parser.add_argument(
        "--output",
        help="Output file for test results (JSON format)",
        default=None
    )
    
    parser.add_argument(
        "--min-connections",
        type=int,
        help="Minimum connections in pool",
        default=None
    )
    
    parser.add_argument(
        "--max-connections",
        type=int,
        help="Maximum connections in pool",
        default=None
    )
    
    args = parser.parse_args()
    
    # Validate URL
    if not args.url:
        print("Error: Convex URL is required. Provide with --url or set CONVEX_URL environment variable.")
        sys.exit(1)
    
    # Configure pool
    pool_config = {
        "url": args.url,
    }
    
    if args.min_connections is not None:
        pool_config["min_connections"] = args.min_connections
    
    if args.max_connections is not None:
        pool_config["max_connections"] = args.max_connections
    
    # Initialize pool
    try:
        pool = get_convex_pool(pool_config)
        logger.info("Convex connection pool initialized successfully")
    except Exception as e:
        logger.error(f"Failed to initialize Convex connection pool: {str(e)}")
        sys.exit(1)
    
    # Run the selected test pattern
    try:
        if args.pattern == "steady":
            results = steady_load_test(
                args.url,
                args.concurrency,
                args.requests,
                args.duration,
                args.error_rate
            )
        elif args.pattern == "ramp":
            results = ramp_load_test(
                args.url,
                args.concurrency,
                args.requests,
                args.duration,
                args.error_rate
            )
        elif args.pattern == "spike":
            base_concurrency = max(1, args.concurrency // 3)
            results = spike_load_test(
                args.url,
                base_concurrency,
                args.concurrency,
                args.requests,
                args.duration,
                args.error_rate
            )
        elif args.pattern == "oscillate":
            min_concurrency = max(1, args.concurrency // 3)
            results = oscillating_load_test(
                args.url,
                min_concurrency,
                args.concurrency,
                args.requests,
                args.duration,
                args.error_rate
            )
        else:
            logger.error(f"Unknown pattern: {args.pattern}")
            sys.exit(1)
        
        # Get final pool stats
        pool_stats = pool.get_stats()
        
        # Generate report
        report = results.generate_report()
        
        # Add test configuration
        report["configuration"] = {
            "pattern": args.pattern,
            "concurrency": args.concurrency,
            "requests": args.requests,
            "duration": args.duration,
            "error_rate": args.error_rate
        }
        
        # Add pool configuration and stats
        report["pool_configuration"] = generate_pool_config_report(pool_config)
        report["pool_stats"] = pool_stats
        
        # Print report
        print("\nTest Results:")
        print(f"  Pattern: {args.pattern}")
        print(f"  Concurrency: {args.concurrency}")
        print(f"  Total Requests: {report['summary']['total_requests']}")
        print(f"  Successful Requests: {report['summary']['successful_requests']}")
        print(f"  Failed Requests: {report['summary']['failed_requests']}")
        print(f"  Test Duration: {report['summary']['duration_seconds']:.2f} seconds")
        print(f"  Requests/second: {report['summary']['requests_per_second']:.2f}")
        print(f"  Success Rate: {report['summary']['success_rate_percent']:.2f}%")
        print(f"  Avg Request Time: {report['timing']['average_request_time']:.3f} seconds")
        print(f"  Min Request Time: {report['timing']['min_request_time']:.3f} seconds")
        print(f"  Max Request Time: {report['timing']['max_request_time']:.3f} seconds")
        
        # Save report to file if specified
        if args.output:
            with open(args.output, 'w') as f:
                json.dump(report, f, indent=2)
            print(f"\nDetailed report saved to {args.output}")
    
    except KeyboardInterrupt:
        print("\nTest interrupted by user")
    except Exception as e:
        logger.error(f"Error during test: {str(e)}")
        logger.debug(traceback.format_exc())
    finally:
        # Close the pool
        try:
            pool.close()
            logger.info("Convex connection pool closed")
        except Exception as e:
            logger.error(f"Error closing connection pool: {str(e)}")


if __name__ == "__main__":
    # Import math here to avoid polluting the global namespace
    import math
    import traceback
    
    main()