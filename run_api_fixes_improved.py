#!/usr/bin/env python3
"""
CryoProtect v2 - Improved API Fixes Runner

This script runs the necessary fixes to make the API fully functional with improved:
1. Monitoring and progress tracking
2. Timeout handling
3. Error recovery
4. Detailed logging

Usage:
    python run_api_fixes_improved.py
"""

import os
import sys
import json
import time
import logging
import threading
import subprocess
from datetime import datetime

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler("api_fixes_improved.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

# Add a performance logger
perf_logger = logging.getLogger("performance")
perf_handler = logging.FileHandler("performance.log")
perf_handler.setFormatter(logging.Formatter("%(asctime)s [PERF] %(message)s"))
perf_logger.addHandler(perf_handler)
perf_logger.setLevel(logging.INFO)

class ProgressMonitor:
    """Monitors and displays progress of the API fixes process."""
    
    def __init__(self):
        self.current_step = ""
        self.progress = 0
        self.total_steps = 3  # Total steps in our process
        self.running = True
        self._lock = threading.Lock()
        
    def update(self, step, progress=None):
        """Update the current step and progress."""
        with self._lock:
            self.current_step = step
            if progress is not None:
                self.progress = progress
            logger.info(f"Progress: {self.progress}/{self.total_steps} - {self.current_step}")
    
    def display_progress(self):
        """Display progress in the console."""
        last_progress = 0
        while self.running:
            with self._lock:
                if self.progress > last_progress:
                    print(f"Step {self.progress}/{self.total_steps} completed: {self.current_step}")
                    last_progress = self.progress
            time.sleep(1.0)
                
    def start(self):
        """Start the progress monitor thread."""
        self.thread = threading.Thread(target=self.display_progress)
        self.thread.daemon = True
        self.thread.start()
        
    def stop(self):
        """Stop the progress monitor thread."""
        self.running = False
        if hasattr(self, 'thread'):
            self.thread.join(timeout=1.0)

def log_performance(func):
    """Decorator to log performance of functions."""
    def wrapper(*args, **kwargs):
        start_time = time.time()
        result = func(*args, **kwargs)
        elapsed = time.time() - start_time
        perf_logger.info(f"{func.__name__} took {elapsed:.2f} seconds")
        return result
    return wrapper

def run_command_with_timeout(command, description, timeout=60):
    """Run a command with a timeout."""
    logger.info(f"Running {description}...")
    
    result = {"success": False, "output": ""}
    process = None
    
    def target():
        nonlocal process
        try:
            process = subprocess.Popen(
                command,
                shell=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True
            )
            stdout, stderr = process.communicate()
            result["success"] = process.returncode == 0
            result["output"] = stdout if process.returncode == 0 else stderr
        except Exception as e:
            result["output"] = str(e)
    
    thread = threading.Thread(target=target)
    thread.start()
    thread.join(timeout)
    
    if thread.is_alive():
        logger.error(f"{description} timed out after {timeout} seconds")
        if process:
            try:
                process.terminate()
                process.wait(timeout=5)
            except:
                process.kill()
        thread.join(5)
        return False, f"{description} timed out after {timeout} seconds"
    
    if result["success"]:
        logger.info(f"{description} completed successfully")
    else:
        logger.error(f"{description} failed: {result['output']}")
        
    return result["success"], result["output"]

@log_performance
def fix_database_tables():
    """Run the database table fix script with timeout."""
    return run_command_with_timeout("python fix_database_tables.py", "Database table fix script", timeout=120)

@log_performance
def run_api_verification():
    """Run the API verification script with timeout."""
    return run_command_with_timeout("python tests/verify_api_endpoints.py", "API verification script", timeout=180)

def create_watchdog(timeout=300):
    """Create a watchdog to terminate the process if it hangs."""
    def watchdog():
        time.sleep(timeout)
        logger.error(f"Watchdog timeout reached after {timeout} seconds. Terminating process.")
        os._exit(1)
    
    watchdog_thread = threading.Thread(target=watchdog)
    watchdog_thread.daemon = True
    watchdog_thread.start()
    return watchdog_thread

def checkpoint(name):
    """Create a checkpoint file to track progress."""
    checkpoint_dir = "checkpoints"
    os.makedirs(checkpoint_dir, exist_ok=True)
    
    with open(os.path.join(checkpoint_dir, f"{name}.txt"), "w") as f:
        f.write(str(time.time()))
    logger.info(f"Checkpoint: {name}")

def diagnose_system():
    """Diagnose the system and log information."""
    try:
        import psutil
        
        # Log CPU and memory usage
        logger.info(f"CPU usage: {psutil.cpu_percent(interval=1.0)}%")
        logger.info(f"Memory usage: {psutil.virtual_memory().percent}%")
        
        # Log disk usage
        disk = psutil.disk_usage('/')
        logger.info(f"Disk usage: {disk.percent}% (Free: {disk.free / (1024**3):.2f} GB)")
        
        # Log network information
        net_io = psutil.net_io_counters()
        logger.info(f"Network: Sent {net_io.bytes_sent / (1024**2):.2f} MB, Received {net_io.bytes_recv / (1024**2):.2f} MB")
        
        return True
    except Exception as e:
        logger.error(f"Error diagnosing system: {e}")
        return False

def main():
    """Main function to run the API fixes with improved error handling and monitoring."""
    print("\n" + "=" * 80)
    print("CryoProtect v2 - Improved API Fixes Runner")
    print("=" * 80)
    
    # Create a watchdog to terminate if the process hangs
    watchdog = create_watchdog(timeout=600)  # 10 minutes timeout
    
    # Create progress monitor
    monitor = ProgressMonitor()
    monitor.start()
    
    try:
        # Diagnose system
        logger.info("Diagnosing system...")
        diagnose_system()
        
        # Create checkpoint
        checkpoint("start")
        
        # Step 1: Fix database tables
        monitor.update("Fixing database tables", 0)
        success, output = fix_database_tables()
        if not success:
            print("Failed to fix database tables. Check the logs for details.")
            return 1
        
        # Create checkpoint
        checkpoint("database_fixed")
        monitor.update("Database tables fixed", 1)
        
        # Step 2: Run API verification
        monitor.update("Verifying API endpoints", 1)
        success, output = run_api_verification()
        if not success:
            print("API verification failed. Check the logs for details.")
            return 1
        
        # Create checkpoint
        checkpoint("api_verified")
        monitor.update("API endpoints verified", 2)
        
        # Step 3: Generate summary report
        monitor.update("Generating summary report", 2)
        
        # Try to parse the verification report
        report_path = None
        for line in output.splitlines():
            if "Report saved to:" in line:
                report_path = line.split("Report saved to:")[1].strip()
                break
        
        if report_path and os.path.exists(report_path):
            try:
                with open(report_path, 'r') as f:
                    report = json.load(f)
                
                print("\n" + "=" * 60)
                print("API Fixes Summary Report")
                print("=" * 60)
                print(f"Status: {report['status']}")
                print(f"Total Endpoints: {report['summary']['total']}")
                print(f"Implemented: {report['summary']['implemented']} ({report['summary']['implemented']/report['summary']['total']*100:.1f}%)")
                print(f"Functional: {report['summary']['functional']} ({report['summary']['functional']/report['summary']['total']*100:.1f}%)")
                print(f"Discrepancies: {report['summary']['discrepancies']}")
                
                if report['status'] == "SUCCESS":
                    print("\nAll API endpoints are now functional!")
                else:
                    print("\nSome API endpoints still have issues. Check the verification report for details.")
                
                # Save a copy of the report with a timestamp
                timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
                report_copy_path = f"api_verification_report_{timestamp}.json"
                with open(report_copy_path, 'w') as f:
                    json.dump(report, f, indent=2)
                print(f"\nA copy of the report has been saved to: {report_copy_path}")
            except Exception as e:
                logger.error(f"Error parsing verification report: {e}")
                print("\nCould not parse verification report. Check the logs for details.")
        else:
            print("\nCould not find verification report. Check the logs for details.")
        
        # Create checkpoint
        checkpoint("complete")
        monitor.update("Summary report generated", 3)
        
        print("\n" + "=" * 60)
        print("API Fixes Complete")
        print("=" * 60)
        
        return 0 if success else 1
    
    except KeyboardInterrupt:
        logger.info("Process interrupted by user")
        print("\nProcess interrupted by user")
        return 130
    except Exception as e:
        logger.error(f"Unexpected error: {e}")
        print(f"\nUnexpected error: {e}")
        return 1
    finally:
        monitor.stop()

if __name__ == "__main__":
    sys.exit(main())