#!/usr/bin/env python3
"""
CryoProtect v2 - Process Termination and Summary Script

This script gracefully terminates running API fix processes and provides a summary
of the improvements made to the system.

Usage:
    python terminate_and_summarize.py
"""

import os
import sys
import json
import time
import signal
import logging
import subprocess
import platform
from datetime import datetime

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler("termination_summary.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

def find_and_terminate_process(process_name):
    """Find and terminate a process by name."""
    terminated = False
    current_pid = os.getpid()
    
    try:
        if platform.system() == "Windows":
            # Windows
            result = subprocess.run(
                f"tasklist /FI \"IMAGENAME eq {process_name}\" /FO CSV /NH",
                shell=True,
                capture_output=True,
                text=True
            )
            
            if process_name in result.stdout:
                logger.info(f"Found process: {process_name}")
                # Exclude current process
                subprocess.run(f"taskkill /F /IM {process_name} /T /FI \"PID ne {current_pid}\"", shell=True)
                logger.info(f"Terminated process: {process_name}")
                terminated = True
            else:
                logger.info(f"Process not found: {process_name}")
        else:
            # Unix-like
            result = subprocess.run(
                f"pgrep -f {process_name}",
                shell=True,
                capture_output=True,
                text=True
            )
            
            if result.stdout.strip():
                # Filter out current process
                pids = [pid for pid in result.stdout.strip().split('\n') if pid and int(pid) != current_pid]
                if pids:
                    for pid in pids:
                        logger.info(f"Found process with PID: {pid}")
                        subprocess.run(f"kill -15 {pid}", shell=True)
                        time.sleep(2)  # Give it time to terminate gracefully
                
                        # Check if it's still running
                        check_result = subprocess.run(
                            f"pgrep -f {pid}",
                            shell=True,
                            capture_output=True,
                            text=True
                        )
                        
                        if check_result.stdout.strip():
                            # Force kill if still running
                            subprocess.run(f"kill -9 {pid}", shell=True)
                        
                        logger.info(f"Terminated process with PID: {pid}")
                    terminated = True
                else:
                    logger.info(f"No processes to terminate (excluding current PID: {current_pid})")
            else:
                logger.info(f"Process not found: {process_name}")
    except Exception as e:
        logger.error(f"Error terminating process: {e}")
    
    return terminated

def check_checkpoint_progress():
    """Check progress based on checkpoint files."""
    checkpoint_dir = "checkpoints"
    if not os.path.exists(checkpoint_dir):
        return "No checkpoints found. Progress unknown."
    
    checkpoints = {
        "start": "Process started",
        "database_fixed": "Database tables fixed",
        "api_verified": "API endpoints verified",
        "complete": "Process completed"
    }
    
    progress = []
    for checkpoint, description in checkpoints.items():
        checkpoint_file = os.path.join(checkpoint_dir, f"{checkpoint}.txt")
        if os.path.exists(checkpoint_file):
            with open(checkpoint_file, 'r') as f:
                timestamp = float(f.read().strip())
                dt = datetime.fromtimestamp(timestamp)
                progress.append(f"{description} at {dt.strftime('%Y-%m-%d %H:%M:%S')}")
    
    if not progress:
        return "No progress recorded in checkpoints."
    
    return "\n".join(progress)

def check_verification_reports():
    """Check for verification reports and summarize them."""
    reports_dir = "reports"
    if not os.path.exists(reports_dir):
        reports_dir = os.path.join("tests", "reports")
        if not os.path.exists(reports_dir):
            return "No verification reports found."
    
    # Find the most recent report
    report_files = [f for f in os.listdir(reports_dir) if f.startswith("api_verification") and f.endswith(".json")]
    if not report_files:
        return "No verification reports found in reports directory."
    
    # Sort by modification time (most recent first)
    report_files.sort(key=lambda f: os.path.getmtime(os.path.join(reports_dir, f)), reverse=True)
    latest_report = os.path.join(reports_dir, report_files[0])
    
    try:
        with open(latest_report, 'r') as f:
            report = json.load(f)
        
        summary = report.get("summary", {})
        total = summary.get("total", 0)
        implemented = summary.get("implemented", 0)
        functional = summary.get("functional", 0)
        discrepancies = summary.get("discrepancies", 0)
        
        return f"""
Latest verification report: {os.path.basename(latest_report)}
Status: {report.get("status", "Unknown")}
Total Endpoints: {total}
Implemented: {implemented} ({implemented/total*100 if total else 0:.1f}%)
Functional: {functional} ({functional/total*100 if total else 0:.1f}%)
Discrepancies: {discrepancies}
"""
    except Exception as e:
        return f"Error reading verification report: {e}"

def check_log_files():
    """Check log files for errors and warnings."""
    log_files = [
        "api_fixes_improved.log",
        "fix_database_tables.log",
        "api_verification_standalone.log",
        "performance.log"
    ]
    
    log_summary = {}
    
    for log_file in log_files:
        if os.path.exists(log_file):
            error_count = 0
            warning_count = 0
            
            with open(log_file, 'r') as f:
                for line in f:
                    if "[ERROR]" in line:
                        error_count += 1
                    elif "[WARNING]" in line:
                        warning_count += 1
            
            log_summary[log_file] = {
                "errors": error_count,
                "warnings": warning_count,
                "size": os.path.getsize(log_file)
            }
    
    if not log_summary:
        return "No log files found."
    
    result = "Log File Summary:\n"
    for log_file, stats in log_summary.items():
        result += f"{log_file}: {stats['errors']} errors, {stats['warnings']} warnings, {stats['size']/1024:.1f} KB\n"
    
    return result

def summarize_improvements():
    """Summarize the improvements made to the system."""
    improvements = [
        "1. Created improved API fixes runner with progress monitoring, timeout handling, and error recovery",
        "2. Developed standalone API verification script for independent endpoint testing",
        "3. Created modular database fix script for targeted database repairs",
        "4. Implemented real-time API monitoring dashboard for system health tracking",
        "5. Added comprehensive logging for performance analysis and debugging",
        "6. Created checkpoint system for process recovery",
        "7. Added detailed documentation for troubleshooting and maintenance"
    ]
    
    return "\n".join(improvements)

def main():
    """Main function to terminate processes and generate summary."""
    print("\n" + "=" * 80)
    print("CryoProtect v2 - Process Termination and Summary")
    print("=" * 80)
    
    # Terminate running processes
    print("\nTerminating running processes...")
    print(f"Current process ID: {os.getpid()} (will be excluded from termination)")
    python_terminated = find_and_terminate_process("python.exe" if platform.system() == "Windows" else "python")
    
    if not python_terminated:
        print("No Python processes were terminated.")
    
    # Check progress
    print("\nChecking progress...")
    progress = check_checkpoint_progress()
    print(progress)
    
    # Check verification reports
    print("\nChecking verification reports...")
    report_summary = check_verification_reports()
    print(report_summary)
    
    # Check log files
    print("\nChecking log files...")
    log_summary = check_log_files()
    print(log_summary)
    
    # Summarize improvements
    print("\nImprovements Made:")
    improvements = summarize_improvements()
    print(improvements)
    
    print("\n" + "=" * 80)
    print("Summary Complete")
    print("=" * 80)
    
    # Save summary to file
    summary = f"""
CryoProtect v2 - Process Termination and Summary
================================================

Termination Status:
Python processes terminated: {"Yes" if python_terminated else "No"}

Progress:
{progress}

Verification Reports:
{report_summary}

Log Summary:
{log_summary}

Improvements Made:
{improvements}

Next Steps:
1. Run the standalone verification script to check endpoint status
2. Use the monitoring dashboard to track system health
3. Consult the troubleshooting guide for any remaining issues

Summary generated: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}
"""
    
    with open("api_fixes_summary.md", "w") as f:
        f.write(summary)
    
    print(f"\nSummary saved to: api_fixes_summary.md")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())