#!/usr/bin/env python3
"""
Run Comprehensive ChEMBL Verification

This script executes the comprehensive ChEMBL verification process and outputs
a user-friendly summary of the results. It's designed to be integrated into
the ChEMBL import workflow or run independently to validate the ChEMBL data
in the CryoProtect database.

Usage:
    python run_comprehensive_chembl_verification.py [--with-api-check] [--with-visualizations]
"""

import os
import sys
import json
import argparse
import subprocess
import webbrowser
from datetime import datetime

def run_verification(api_check=False, visualizations=False, sample_size=50):
    """
    Run the verification process.
    
    Args:
        api_check: Whether to check against the ChEMBL API
        visualizations: Whether to generate visualizations
        sample_size: Sample size for data quality checks
        
    Returns:
        Tuple of (success_flag, report_path)
    """
    cmd = [sys.executable, "comprehensive_chembl_verification.py"]
    
    # Add optional arguments
    if api_check:
        cmd.append("--check_with_chembl_api")
    
    if visualizations:
        cmd.append("--generate_visualizations")
    
    if sample_size:
        cmd.extend(["--sample_size", str(sample_size)])
    
    # Run the verification script
    print(f"Running: {' '.join(cmd)}")
    print("=" * 80)
    
    try:
        result = subprocess.run(cmd, check=True)
        
        # Find the latest HTML report
        reports_dir = "reports"
        html_reports = [f for f in os.listdir(reports_dir) if f.startswith("chembl_verification_") and f.endswith(".html")]
        if html_reports:
            # Sort by filename (which includes timestamp) to get the latest
            latest_report = sorted(html_reports)[-1]
            report_path = os.path.join(reports_dir, latest_report)
            return True, report_path
        else:
            print("No HTML report found")
            return True, None
    
    except subprocess.CalledProcessError as e:
        print(f"Verification process failed with exit code {e.returncode}")
        return False, None
    
    except Exception as e:
        print(f"Error running verification: {str(e)}")
        return False, None

def open_report(report_path):
    """
    Open the HTML report in the default web browser.
    
    Args:
        report_path: Path to the HTML report
    """
    if report_path and os.path.exists(report_path):
        try:
            print(f"Opening report: {report_path}")
            webbrowser.open(f"file://{os.path.abspath(report_path)}")
        except Exception as e:
            print(f"Error opening report: {str(e)}")
            print(f"Please open the report manually: {report_path}")
    else:
        print("No report to open")

def print_summary(success_flag, report_path):
    """
    Print a summary of the verification results.
    
    Args:
        success_flag: Whether the verification was successful
        report_path: Path to the verification report
    """
    # Try to read the report for detailed summary
    if report_path and os.path.exists(report_path.replace(".html", ".json")):
        try:
            with open(report_path.replace(".html", ".json"), "r") as f:
                report_data = json.load(f)
            
            # Extract key metrics
            molecule_count = report_data.get("data_completeness", {}).get("molecule_count", {}).get("with_chembl_id", 0)
            reference_count = len(report_data.get("data_completeness", {}).get("reference_compounds", {}).get("found", []))
            total_reference = len(report_data.get("data_completeness", {}).get("reference_compounds", {}).get("required", []))
            
            # If we have property coverage data, calculate the average
            property_coverage = 0
            property_count = 0
            properties = report_data.get("data_completeness", {}).get("property_coverage", {}).get("properties", {})
            for prop_name, data in properties.items():
                property_coverage += data.get("percentage", 0)
                property_count += 1
            
            avg_property_coverage = property_coverage / property_count if property_count > 0 else 0
            
            # Get cross-reference coverage
            pubchem_percentage = report_data.get("cross_references", {}).get("pubchem_coverage", {}).get("percentage", 0)
            inchikey_percentage = report_data.get("cross_references", {}).get("inchikey_coverage", {}).get("percentage", 0)
            
            # Print detailed summary
            print("\n" + "=" * 80)
            print("CHEMBL VERIFICATION SUMMARY")
            print("=" * 80)
            print(f"Status: {report_data.get('summary', {}).get('overall_result', 'UNKNOWN')}")
            print(f"ChEMBL Molecules: {molecule_count}")
            print(f"Reference Compounds: {reference_count}/{total_reference}")
            print(f"Average Property Coverage: {avg_property_coverage:.1f}%")
            print(f"PubChem Cross-Reference: {pubchem_percentage:.1f}%")
            print(f"InChIKey Coverage: {inchikey_percentage:.1f}%")
            print("=" * 80)
            print(f"Detailed report: {report_path}")
            print("=" * 80)
            
            return
        
        except Exception as e:
            print(f"Error reading report data: {str(e)}")
    
    # Fallback summary
    print("\n" + "=" * 80)
    print("CHEMBL VERIFICATION SUMMARY")
    print("=" * 80)
    print(f"Status: {'PASSED' if success_flag else 'FAILED'}")
    print(f"Report: {report_path if report_path else 'Not available'}")
    print("=" * 80)

def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(description="Run Comprehensive ChEMBL Verification")
    parser.add_argument("--with-api-check", action="store_true", help="Check data against ChEMBL API")
    parser.add_argument("--with-visualizations", action="store_true", help="Generate visualizations")
    parser.add_argument("--sample-size", type=int, default=50, help="Sample size for data quality checks")
    parser.add_argument("--open-report", action="store_true", help="Open HTML report after verification")
    args = parser.parse_args()
    
    print(f"Starting ChEMBL verification at {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    
    # Create necessary directories
    os.makedirs("reports", exist_ok=True)
    os.makedirs("logs", exist_ok=True)
    
    # Run verification
    success_flag, report_path = run_verification(
        api_check=args.with_api_check,
        visualizations=args.with_visualizations,
        sample_size=args.sample_size
    )
    
    # Print summary
    print_summary(success_flag, report_path)
    
    # Open report if requested
    if args.open_report and report_path:
        open_report(report_path)
    
    # Exit with appropriate code
    return 0 if success_flag else 1

if __name__ == "__main__":
    sys.exit(main())