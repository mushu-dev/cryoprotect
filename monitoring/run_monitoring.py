#!/usr/bin/env python3
"""
CryoProtect v2 - Run Monitoring Stack

This script provides commands to start, stop, and manage the monitoring stack
using Docker Compose.
"""

import os
import sys
import subprocess
import argparse
import time
import logging

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

# Get the directory of this script
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

def run_command(command, cwd=None):
    """
    Run a shell command and return the output.
    
    Args:
        command: Command to run
        cwd: Working directory
        
    Returns:
        Tuple of (return_code, stdout, stderr)
    """
    if cwd is None:
        cwd = SCRIPT_DIR
        
    logger.info(f"Running command: {command}")
    
    process = subprocess.Popen(
        command,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        shell=True,
        cwd=cwd
    )
    
    stdout, stderr = process.communicate()
    return_code = process.returncode
    
    stdout = stdout.decode('utf-8')
    stderr = stderr.decode('utf-8')
    
    if return_code != 0:
        logger.error(f"Command failed with return code {return_code}")
        logger.error(f"stderr: {stderr}")
    
    return return_code, stdout, stderr

def start_monitoring():
    """Start the monitoring stack."""
    logger.info("Starting monitoring stack...")
    
    # Check if Docker is running
    return_code, _, _ = run_command("docker info")
    if return_code != 0:
        logger.error("Docker is not running. Please start Docker and try again.")
        return False
    
    # Start the monitoring stack
    return_code, stdout, stderr = run_command("docker-compose up -d")
    if return_code != 0:
        logger.error("Failed to start monitoring stack.")
        return False
    
    logger.info("Monitoring stack started successfully.")
    logger.info("Grafana: http://localhost:3000 (admin/cryoprotect)")
    logger.info("Prometheus: http://localhost:9090")
    logger.info("Alertmanager: http://localhost:9093")
    
    return True

def stop_monitoring():
    """Stop the monitoring stack."""
    logger.info("Stopping monitoring stack...")
    
    # Stop the monitoring stack
    return_code, stdout, stderr = run_command("docker-compose down")
    if return_code != 0:
        logger.error("Failed to stop monitoring stack.")
        return False
    
    logger.info("Monitoring stack stopped successfully.")
    return True

def restart_monitoring():
    """Restart the monitoring stack."""
    logger.info("Restarting monitoring stack...")
    
    # Stop the monitoring stack
    stop_success = stop_monitoring()
    if not stop_success:
        logger.error("Failed to stop monitoring stack. Cannot restart.")
        return False
    
    # Wait a moment for containers to stop
    time.sleep(2)
    
    # Start the monitoring stack
    start_success = start_monitoring()
    if not start_success:
        logger.error("Failed to start monitoring stack. Restart failed.")
        return False
    
    logger.info("Monitoring stack restarted successfully.")
    return True

def status_monitoring():
    """Check the status of the monitoring stack."""
    logger.info("Checking monitoring stack status...")
    
    # Check the status of the monitoring stack
    return_code, stdout, stderr = run_command("docker-compose ps")
    if return_code != 0:
        logger.error("Failed to check monitoring stack status.")
        return False
    
    logger.info("Monitoring stack status:")
    print(stdout)
    
    return True

def logs_monitoring(service=None, follow=False):
    """
    Show logs for the monitoring stack.
    
    Args:
        service: Service to show logs for (prometheus, grafana, alertmanager, node-exporter)
        follow: Whether to follow the logs
    """
    if service:
        logger.info(f"Showing logs for {service}...")
        command = f"docker-compose logs {'-f' if follow else ''} {service}"
    else:
        logger.info("Showing logs for all services...")
        command = f"docker-compose logs {'-f' if follow else ''}"
    
    # Show the logs
    return_code, stdout, stderr = run_command(command)
    if return_code != 0:
        logger.error("Failed to show logs.")
        return False
    
    print(stdout)
    
    return True

def main():
    """Main function."""
    parser = argparse.ArgumentParser(description="CryoProtect v2 - Run Monitoring Stack")
    
    # Add subparsers for commands
    subparsers = parser.add_subparsers(dest="command", help="Command to run")
    
    # Start command
    start_parser = subparsers.add_parser("start", help="Start the monitoring stack")
    
    # Stop command
    stop_parser = subparsers.add_parser("stop", help="Stop the monitoring stack")
    
    # Restart command
    restart_parser = subparsers.add_parser("restart", help="Restart the monitoring stack")
    
    # Status command
    status_parser = subparsers.add_parser("status", help="Check the status of the monitoring stack")
    
    # Logs command
    logs_parser = subparsers.add_parser("logs", help="Show logs for the monitoring stack")
    logs_parser.add_argument("--service", help="Service to show logs for (prometheus, grafana, alertmanager, node-exporter)")
    logs_parser.add_argument("--follow", "-f", action="store_true", help="Follow the logs")
    
    # Parse arguments
    args = parser.parse_args()
    
    # Change to the script directory
    os.chdir(SCRIPT_DIR)
    
    # Run the appropriate command
    if args.command == "start":
        start_monitoring()
    elif args.command == "stop":
        stop_monitoring()
    elif args.command == "restart":
        restart_monitoring()
    elif args.command == "status":
        status_monitoring()
    elif args.command == "logs":
        logs_monitoring(args.service, args.follow)
    else:
        parser.print_help()

if __name__ == "__main__":
    main()