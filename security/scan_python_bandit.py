#!/usr/bin/env python3
"""
scan_python_bandit.py - Python security scanner using Bandit

This script scans Python code for security vulnerabilities using Bandit.
It can be run as a standalone script or integrated into CI/CD pipelines.

Usage:
    python scan_python_bandit.py [options]

Options:
    --path PATH           Path to scan (default: current directory)
    --format FORMAT       Output format (json, yaml, csv, txt) (default: txt)
    --output FILE         Output file (default: bandit-results.{format})
    --severity SEVERITY   Minimum severity to report (low, medium, high) (default: low)
    --confidence CONF     Minimum confidence to report (low, medium, high) (default: low)
    --exit-on-critical    Exit with code 1 if critical issues found
    --exclude PATTERN     Comma-separated list of paths to exclude
    --help                Show this help message and exit
"""

import argparse
import json
import os
import subprocess
import sys
from datetime import datetime


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="Python security scanner using Bandit")
    parser.add_argument("--path", default=".", help="Path to scan (default: current directory)")
    parser.add_argument("--format", default="txt", choices=["json", "yaml", "csv", "txt"],
                        help="Output format (json, yaml, csv, txt) (default: txt)")
    parser.add_argument("--output", default=None, help="Output file (default: bandit-results.{format})")
    parser.add_argument("--severity", default="low", choices=["low", "medium", "high"],
                        help="Minimum severity to report (low, medium, high) (default: low)")
    parser.add_argument("--confidence", default="low", choices=["low", "medium", "high"],
                        help="Minimum confidence to report (low, medium, high) (default: low)")
    parser.add_argument("--exit-on-critical", action="store_true", help="Exit with code 1 if critical issues found")
    parser.add_argument("--exclude", default="", help="Comma-separated list of paths to exclude")
    return parser.parse_args()


def check_bandit_installed():
    """Check if Bandit is installed."""
    try:
        subprocess.run(["bandit", "--version"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
        return True
    except (subprocess.SubprocessError, FileNotFoundError):
        return False


def install_bandit():
    """Install Bandit if not already installed."""
    print("Bandit not found. Installing...")
    try:
        subprocess.run([sys.executable, "-m", "pip", "install", "bandit"], check=True)
        print("Bandit installed successfully.")
        return True
    except subprocess.SubprocessError:
        print("Failed to install Bandit. Please install it manually: pip install bandit")
        return False


def run_bandit_scan(args):
    """Run Bandit scan with the specified arguments."""
    if not check_bandit_installed() and not install_bandit():
        return False, None

    # Set default output file if not provided
    output_file = args.output or f"bandit-results.{args.format}"

    # Build command
    cmd = ["bandit", "-r", args.path, "-f", args.format, "-o", output_file,
           "--severity-level", args.severity, "--confidence-level", args.confidence]

    # Add exclude paths if provided
    if args.exclude:
        cmd.extend(["--exclude", args.exclude])

    print(f"Running Bandit scan on {args.path}...")
    try:
        subprocess.run(cmd, check=True)
        print(f"Scan completed successfully. Results saved to {output_file}")
        return True, output_file
    except subprocess.SubprocessError as e:
        print(f"Error running Bandit scan: {e}")
        return False, None


def check_for_critical_issues(output_file, output_format):
    """Check if there are any critical issues in the scan results."""
    if not os.path.exists(output_file):
        print(f"Output file {output_file} not found.")
        return False

    if output_format == "json":
        try:
            with open(output_file, 'r') as f:
                results = json.load(f)
            
            high_severity_count = 0
            for result in results.get("results", []):
                if result.get("issue_severity") == "HIGH":
                    high_severity_count += 1
            
            if high_severity_count > 0:
                print(f"Found {high_severity_count} high severity issues!")
                return True
        except (json.JSONDecodeError, KeyError) as e:
            print(f"Error parsing JSON results: {e}")
    
    # For other formats, we'll need to parse the file differently or rely on Bandit's exit code
    return False


def generate_summary(output_file, output_format, scan_success):
    """Generate a summary of the scan results."""
    summary = {
        "scanner": "bandit",
        "timestamp": datetime.utcnow().isoformat(),
        "success": scan_success,
        "output_file": output_file,
        "output_format": output_format
    }

    if scan_success and output_format == "json" and os.path.exists(output_file):
        try:
            with open(output_file, 'r') as f:
                results = json.load(f)
            
            # Extract metrics
            metrics = results.get("metrics", {})
            summary["total_files"] = metrics.get("_totals", {}).get("loc", 0)
            
            # Count issues by severity
            severity_counts = {"HIGH": 0, "MEDIUM": 0, "LOW": 0}
            for result in results.get("results", []):
                severity = result.get("issue_severity")
                if severity in severity_counts:
                    severity_counts[severity] += 1
            
            summary["issues"] = {
                "high": severity_counts["HIGH"],
                "medium": severity_counts["MEDIUM"],
                "low": severity_counts["LOW"],
                "total": sum(severity_counts.values())
            }
        except (json.JSONDecodeError, KeyError) as e:
            print(f"Error generating summary: {e}")

    # Save summary to file
    summary_file = f"bandit-summary-{datetime.utcnow().strftime('%Y%m%d%H%M%S')}.json"
    with open(summary_file, 'w') as f:
        json.dump(summary, f, indent=2)
    
    print(f"Summary saved to {summary_file}")
    return summary


def main():
    """Main function."""
    args = parse_args()
    scan_success, output_file = run_bandit_scan(args)
    
    if scan_success and output_file:
        summary = generate_summary(output_file, args.format, scan_success)
        
        if args.exit_on_critical and check_for_critical_issues(output_file, args.format):
            print("Exiting with code 1 due to critical issues found.")
            sys.exit(1)
    
    sys.exit(0 if scan_success else 1)


if __name__ == "__main__":
    main()