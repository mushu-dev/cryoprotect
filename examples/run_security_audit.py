"""
Example script demonstrating how to use the CryoProtect v2 Security Audit System.

This script runs a security audit on the project and generates both JSON and Markdown reports.
It demonstrates how to use the security auditor programmatically and how to interpret the results.

Usage:
    python examples/run_security_audit.py
"""

import os
import sys
import json
from pathlib import Path

# Add the project root directory to the Python path
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

# Import the security auditor
from security.auditor import run_audit

def main():
    """Run a security audit and display the results."""
    print("Running CryoProtect v2 Security Audit...")
    
    # Create reports directory if it doesn't exist
    reports_dir = project_root / "reports"
    reports_dir.mkdir(exist_ok=True)
    
    # Define output file paths
    json_report_path = reports_dir / "security_audit_report.json"
    md_report_path = reports_dir / "security_audit_report.md"
    
    # Run the security audit
    report = run_audit(
        base_dir=str(project_root),
        output_json=str(json_report_path),
        output_md=str(md_report_path)
    )
    
    # Display a summary of the results
    print("\nSecurity Audit Summary:")
    print(f"Total findings: {report['summary']['total_findings']}")
    print(f"Risk level: {report['summary']['risk_level']} (Score: {report['summary']['risk_score']})")
    
    print("\nFindings by severity:")
    for severity, count in report['summary']['severity_counts'].items():
        print(f"  {severity}: {count}")
    
    print("\nFindings by category:")
    for category, count in sorted(report['summary']['category_counts'].items(), key=lambda x: x[1], reverse=True):
        print(f"  {category}: {count}")
    
    print("\nTop 5 critical/high severity findings:")
    critical_high_findings = [f for f in report['findings'] if f['severity'] in ('CRITICAL', 'HIGH')]
    for i, finding in enumerate(critical_high_findings[:5], 1):
        print(f"  {i}. {finding['title']} ({finding['severity']})")
        print(f"     Location: {finding['location']}")
    
    print("\nReports generated:")
    print(f"  JSON report: {json_report_path}")
    print(f"  Markdown report: {md_report_path}")
    
    print("\nNext steps:")
    print("1. Review the detailed findings in the Markdown report")
    print("2. Prioritize remediation based on severity and impact")
    print("3. Implement the recommended fixes")
    print("4. Re-run the audit to verify that issues have been resolved")

if __name__ == "__main__":
    main()