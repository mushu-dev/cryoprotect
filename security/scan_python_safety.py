#!/usr/bin/env python3
"""
scan_python_safety.py - Python dependency vulnerability scanner using Safety

This script scans Python dependencies for known security vulnerabilities using Safety.
It can be run as a standalone script or integrated into CI/CD pipelines.

Usage:
    python scan_python_safety.py [options]

Options:
    --requirements FILE    Requirements file to scan (default: requirements.txt)
    --format FORMAT        Output format (json, text, bare, full) (default: text)
    --output FILE          Output file (default: safety-results.{format})
    --exit-on-critical     Exit with code 1 if critical vulnerabilities found
    --api-key KEY          Safety API key for full database access
    --help                 Show this help message and exit
"""

import argparse
import json
import os
import subprocess
import sys
from datetime import datetime


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="Python dependency vulnerability scanner using Safety")
    parser.add_argument("--requirements", default="requirements.txt", 
                        help="Requirements file to scan (default: requirements.txt)")
    parser.add_argument("--format", default="text", choices=["json", "text", "bare", "full"],
                        help="Output format (json, text, bare, full) (default: text)")
    parser.add_argument("--output", default=None, 
                        help="Output file (default: safety-results.{format})")
    parser.add_argument("--exit-on-critical", action="store_true", 
                        help="Exit with code 1 if critical vulnerabilities found")
    parser.add_argument("--api-key", default=None, 
                        help="Safety API key for full database access")
    return parser.parse_args()


def check_safety_installed():
    """Check if Safety is installed."""
    try:
        subprocess.run(["safety", "--version"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
        return True
    except (subprocess.SubprocessError, FileNotFoundError):
        return False


def install_safety():
    """Install Safety if not already installed."""
    print("Safety not found. Installing...")
    try:
        subprocess.run([sys.executable, "-m", "pip", "install", "safety"], check=True)
        print("Safety installed successfully.")
        return True
    except subprocess.SubprocessError:
        print("Failed to install Safety. Please install it manually: pip install safety")
        return False


def run_safety_scan(args):
    """Run Safety scan with the specified arguments."""
    if not check_safety_installed() and not install_safety():
        return False, None

    # Check if requirements file exists
    if not os.path.exists(args.requirements):
        print(f"Requirements file {args.requirements} not found.")
        return False, None

    # Set default output file if not provided
    output_format = "json" if args.format == "json" else "txt"
    output_file = args.output or f"safety-results.{output_format}"

    # Build command
    cmd = ["safety", "check", "-r", args.requirements, "--output", args.format]
    
    # Add API key if provided
    if args.api_key:
        cmd.extend(["--key", args.api_key])

    print(f"Running Safety scan on {args.requirements}...")
    try:
        # Redirect output to file
        with open(output_file, 'w') as f:
            subprocess.run(cmd, stdout=f, stderr=subprocess.PIPE, check=True)
        
        print(f"Scan completed successfully. Results saved to {output_file}")
        return True, output_file
    except subprocess.SubprocessError as e:
        # Safety returns non-zero exit code when vulnerabilities are found
        # We'll capture the output and continue
        with open(output_file, 'w') as f:
            result = subprocess.run(cmd, stdout=f, stderr=subprocess.PIPE)
        
        print(f"Scan completed with vulnerabilities found. Results saved to {output_file}")
        return True, output_file


def check_for_critical_vulnerabilities(output_file, output_format):
    """Check if there are any critical vulnerabilities in the scan results."""
    if not os.path.exists(output_file):
        print(f"Output file {output_file} not found.")
        return False

    if output_format == "json":
        try:
            with open(output_file, 'r') as f:
                results = json.load(f)
            
            critical_count = 0
            for vuln in results:
                # Safety doesn't have a standard severity rating, but we can check for critical keywords
                if any(kw in vuln.get("advisory", "").lower() for kw in ["critical", "high", "severe"]):
                    critical_count += 1
            
            if critical_count > 0:
                print(f"Found {critical_count} potentially critical vulnerabilities!")
                return True
        except (json.JSONDecodeError, KeyError, TypeError) as e:
            print(f"Error parsing JSON results: {e}")
    
    # For other formats, we'll need to parse the file differently
    # or rely on Safety's exit code (which is non-zero when vulnerabilities are found)
    return False


def generate_summary(output_file, output_format, scan_success):
    """Generate a summary of the scan results."""
    summary = {
        "scanner": "safety",
        "timestamp": datetime.utcnow().isoformat(),
        "success": scan_success,
        "output_file": output_file,
        "output_format": output_format
    }

    if scan_success and output_format == "json" and os.path.exists(output_file):
        try:
            with open(output_file, 'r') as f:
                results = json.load(f)
            
            # Count vulnerabilities
            summary["vulnerabilities"] = {
                "total": len(results)
            }
            
            # Try to categorize by severity (Safety doesn't provide standard severity ratings)
            severity_counts = {"critical": 0, "high": 0, "medium": 0, "low": 0}
            for vuln in results:
                advisory = vuln.get("advisory", "").lower()
                if "critical" in advisory:
                    severity_counts["critical"] += 1
                elif "high" in advisory:
                    severity_counts["high"] += 1
                elif "medium" in advisory:
                    severity_counts["medium"] += 1
                else:
                    severity_counts["low"] += 1
            
            summary["vulnerabilities"].update(severity_counts)
        except (json.JSONDecodeError, KeyError, TypeError) as e:
            print(f"Error generating summary: {e}")

    # Save summary to file
    summary_file = f"safety-summary-{datetime.utcnow().strftime('%Y%m%d%H%M%S')}.json"
    with open(summary_file, 'w') as f:
        json.dump(summary, f, indent=2)
    
    print(f"Summary saved to {summary_file}")
    return summary


def main():
    """Main function."""
    args = parse_args()
    scan_success, output_file = run_safety_scan(args)
    
    if scan_success and output_file:
        summary = generate_summary(output_file, args.format, scan_success)
        
        if args.exit_on_critical and check_for_critical_vulnerabilities(output_file, args.format):
            print("Exiting with code 1 due to critical vulnerabilities found.")
            sys.exit(1)
    
    sys.exit(0 if scan_success else 1)


if __name__ == "__main__":
    main()