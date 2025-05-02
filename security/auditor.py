"""
CryoProtect v2 Security Auditor

This module provides the core functionality for running security audits on the CryoProtect v2 system.
It coordinates the execution of various audit modules and generates comprehensive reports in both
JSON and Markdown formats.

Usage:
    python security/auditor.py --dir /path/to/project --output-json reports/security_audit_report.json --output-md reports/security_audit_report.md
"""

import os
import sys
import json
import logging
import argparse
import datetime
from typing import Dict, List, Any, Optional, Tuple
from enum import Enum
from pathlib import Path
import importlib.util

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Import audit modules
try:
    from security.audit_modules.authentication import run_audit as run_auth_audit
    from security.audit_modules.api import run_audit as run_api_audit
    from security.audit_modules.data_protection import run_audit as run_data_protection_audit
except ImportError as e:
    logger.error(f"Error importing audit modules: {e}")
    logger.error("Make sure you're running the script from the project root directory.")
    sys.exit(1)

class SeverityLevel(Enum):
    """Severity levels for security findings"""
    CRITICAL = 5
    HIGH = 4
    MEDIUM = 3
    LOW = 2
    INFO = 1
    
    def __str__(self):
        return self.name

class SecurityAuditor:
    """
    Core security auditor that coordinates the execution of various audit modules
    and generates comprehensive reports.
    """
    
    def __init__(self, base_dir: str = None):
        """
        Initialize the security auditor.
        
        Args:
            base_dir: Base directory of the project to audit
        """
        self.base_dir = base_dir or os.getcwd()
        self.findings = []
        self.audit_summary = {}
        logger.info(f"Initialized SecurityAuditor with base directory: {self.base_dir}")
    
    def run_all_audits(self) -> Dict[str, Any]:
        """
        Run all security audit modules and aggregate the findings.
        
        Returns:
            Dictionary containing all findings and summary information
        """
        logger.info("Running all security audits")
        
        # Clear previous findings
        self.findings = []
        
        # Run all audit modules
        auth_findings = run_auth_audit(self.base_dir)
        api_findings = run_api_audit(self.base_dir)
        data_protection_findings = run_data_protection_audit(self.base_dir)
        
        # Aggregate findings
        self.findings = auth_findings + api_findings + data_protection_findings
        
        # Sort findings by severity (highest first)
        self.findings.sort(key=lambda x: x['severity_value'], reverse=True)
        
        # Generate audit summary
        self._generate_audit_summary()
        
        # Create the final report
        report = {
            "audit_info": {
                "timestamp": datetime.datetime.now().isoformat(),
                "project_directory": self.base_dir,
                "modules_executed": ["authentication", "api", "data_protection"]
            },
            "summary": self.audit_summary,
            "findings": self.findings
        }
        
        logger.info(f"Completed all security audits. Found {len(self.findings)} issues.")
        return report
    
    def _generate_audit_summary(self) -> None:
        """Generate a summary of the audit findings."""
        # Count findings by severity
        severity_counts = {level.name: 0 for level in SeverityLevel}
        for finding in self.findings:
            severity_counts[finding['severity']] += 1
        
        # Count findings by category
        category_counts = {}
        for finding in self.findings:
            category = finding['category']
            if category not in category_counts:
                category_counts[category] = 0
            category_counts[category] += 1
        
        # Count findings by module
        module_counts = {}
        for finding in self.findings:
            module = finding['module']
            if module not in module_counts:
                module_counts[module] = 0
            module_counts[module] += 1
        
        # Generate overall risk score (weighted average of severity counts)
        total_findings = len(self.findings)
        if total_findings > 0:
            risk_score = (
                severity_counts['CRITICAL'] * 10 +
                severity_counts['HIGH'] * 5 +
                severity_counts['MEDIUM'] * 2 +
                severity_counts['LOW'] * 1
            ) / total_findings
        else:
            risk_score = 0
        
        # Determine risk level
        if risk_score >= 7:
            risk_level = "Critical"
        elif risk_score >= 4:
            risk_level = "High"
        elif risk_score >= 2:
            risk_level = "Medium"
        elif risk_score > 0:
            risk_level = "Low"
        else:
            risk_level = "None"
        
        # Create summary
        self.audit_summary = {
            "total_findings": total_findings,
            "severity_counts": severity_counts,
            "category_counts": category_counts,
            "module_counts": module_counts,
            "risk_score": round(risk_score, 2),
            "risk_level": risk_level
        }
    
    def generate_json_report(self, output_file: str) -> None:
        """
        Generate a JSON report of the audit findings.
        
        Args:
            output_file: Path to the output JSON file
        """
        logger.info(f"Generating JSON report: {output_file}")
        
        # Run audits if not already run
        if not self.findings:
            self.run_all_audits()
        
        # Create the report
        report = {
            "audit_info": {
                "timestamp": datetime.datetime.now().isoformat(),
                "project_directory": self.base_dir,
                "modules_executed": ["authentication", "api", "data_protection"]
            },
            "summary": self.audit_summary,
            "findings": self.findings
        }
        
        # Write the report to file
        os.makedirs(os.path.dirname(output_file), exist_ok=True)
        with open(output_file, 'w') as f:
            json.dump(report, f, indent=2)
        
        logger.info(f"JSON report generated: {output_file}")
    
    def generate_markdown_report(self, output_file: str) -> None:
        """
        Generate a Markdown report of the audit findings.
        
        Args:
            output_file: Path to the output Markdown file
        """
        logger.info(f"Generating Markdown report: {output_file}")
        
        # Run audits if not already run
        if not self.findings:
            self.run_all_audits()
        
        # Create the report
        timestamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        
        md_content = [
            "# CryoProtect v2 Security Audit Report",
            "",
            f"**Date:** {timestamp}",
            f"**Project Directory:** {self.base_dir}",
            f"**Modules Executed:** authentication, api, data_protection",
            "",
            "## Executive Summary",
            "",
            f"**Total Findings:** {self.audit_summary['total_findings']}",
            f"**Risk Level:** {self.audit_summary['risk_level']} (Score: {self.audit_summary['risk_score']})",
            "",
            "### Findings by Severity",
            ""
        ]
        
        # Add severity counts
        for severity, count in self.audit_summary['severity_counts'].items():
            md_content.append(f"- **{severity}:** {count}")
        
        md_content.extend([
            "",
            "### Findings by Category",
            ""
        ])
        
        # Add category counts
        for category, count in sorted(self.audit_summary['category_counts'].items(), key=lambda x: x[1], reverse=True):
            md_content.append(f"- **{category}:** {count}")
        
        md_content.extend([
            "",
            "## Critical Findings",
            ""
        ])
        
        # Add critical findings
        critical_findings = [f for f in self.findings if f['severity'] == 'CRITICAL']
        if critical_findings:
            for i, finding in enumerate(critical_findings, 1):
                md_content.extend([
                    f"### {i}. {finding['title']}",
                    "",
                    f"**Category:** {finding['category']}",
                    f"**Location:** {finding['location']}",
                    "",
                    f"{finding['description']}",
                    "",
                    "**Recommendation:**",
                    f"{finding['recommendation']}",
                    ""
                ])
        else:
            md_content.append("No critical findings.")
            md_content.append("")
        
        md_content.extend([
            "## High Severity Findings",
            ""
        ])
        
        # Add high severity findings
        high_findings = [f for f in self.findings if f['severity'] == 'HIGH']
        if high_findings:
            for i, finding in enumerate(high_findings, 1):
                md_content.extend([
                    f"### {i}. {finding['title']}",
                    "",
                    f"**Category:** {finding['category']}",
                    f"**Location:** {finding['location']}",
                    "",
                    f"{finding['description']}",
                    "",
                    "**Recommendation:**",
                    f"{finding['recommendation']}",
                    ""
                ])
        else:
            md_content.append("No high severity findings.")
            md_content.append("")
        
        md_content.extend([
            "## Medium Severity Findings",
            ""
        ])
        
        # Add medium severity findings
        medium_findings = [f for f in self.findings if f['severity'] == 'MEDIUM']
        if medium_findings:
            for i, finding in enumerate(medium_findings, 1):
                md_content.extend([
                    f"### {i}. {finding['title']}",
                    "",
                    f"**Category:** {finding['category']}",
                    f"**Location:** {finding['location']}",
                    "",
                    f"{finding['description']}",
                    "",
                    "**Recommendation:**",
                    f"{finding['recommendation']}",
                    ""
                ])
        else:
            md_content.append("No medium severity findings.")
            md_content.append("")
        
        md_content.extend([
            "## Low Severity Findings",
            ""
        ])
        
        # Add low severity findings
        low_findings = [f for f in self.findings if f['severity'] == 'LOW']
        if low_findings:
            for i, finding in enumerate(low_findings, 1):
                md_content.extend([
                    f"### {i}. {finding['title']}",
                    "",
                    f"**Category:** {finding['category']}",
                    f"**Location:** {finding['location']}",
                    "",
                    f"{finding['description']}",
                    "",
                    "**Recommendation:**",
                    f"{finding['recommendation']}",
                    ""
                ])
        else:
            md_content.append("No low severity findings.")
            md_content.append("")
        
        md_content.extend([
            "## Informational Findings",
            ""
        ])
        
        # Add informational findings
        info_findings = [f for f in self.findings if f['severity'] == 'INFO']
        if info_findings:
            for i, finding in enumerate(info_findings, 1):
                md_content.extend([
                    f"### {i}. {finding['title']}",
                    "",
                    f"**Category:** {finding['category']}",
                    f"**Location:** {finding['location']}",
                    "",
                    f"{finding['description']}",
                    "",
                    "**Recommendation:**",
                    f"{finding['recommendation']}",
                    ""
                ])
        else:
            md_content.append("No informational findings.")
            md_content.append("")
        
        md_content.extend([
            "## Appendix: OWASP Top 10 Coverage",
            "",
            "This security audit covers the following OWASP Top 10 (2017) categories:",
            "",
            "1. **A1:2017 Injection** - Checks for SQL, NoSQL, and command injection vulnerabilities",
            "2. **A2:2017 Broken Authentication** - Checks for authentication and session management issues",
            "3. **A3:2017 Sensitive Data Exposure** - Checks for proper encryption and data protection",
            "4. **A4:2017 XML External Entities (XXE)** - Checks for XXE vulnerabilities in XML processing",
            "5. **A5:2017 Broken Access Control** - Checks for proper access control mechanisms",
            "6. **A6:2017 Security Misconfiguration** - Checks for secure configuration of components",
            "7. **A7:2017 Cross-Site Scripting (XSS)** - Checks for XSS vulnerabilities",
            "8. **A8:2017 Insecure Deserialization** - Checks for insecure deserialization",
            "9. **A9:2017 Using Components with Known Vulnerabilities** - Checks for vulnerable dependencies",
            "10. **A10:2017 Insufficient Logging & Monitoring** - Checks for proper logging and monitoring",
            "",
            "## Methodology",
            "",
            "This security audit was performed using static code analysis techniques to identify potential security vulnerabilities in the codebase. The audit focused on identifying common security issues based on the OWASP Top 10 and other security best practices.",
            "",
            "The audit process involved:",
            "",
            "1. **Code scanning** - Automated scanning of the codebase for security patterns",
            "2. **Pattern matching** - Identification of common security anti-patterns",
            "3. **Configuration analysis** - Review of configuration files for security issues",
            "4. **Dependency analysis** - Review of dependencies for known vulnerabilities",
            "",
            "## Disclaimer",
            "",
            "This security audit report is based on automated analysis and may contain false positives or miss certain vulnerabilities. It should be used as a starting point for a more comprehensive security review, including manual code review and penetration testing.",
            "",
            f"*Report generated on {timestamp}*"
        ])
        
        # Write the report to file
        os.makedirs(os.path.dirname(output_file), exist_ok=True)
        with open(output_file, 'w') as f:
            f.write('\n'.join(md_content))
        
        logger.info(f"Markdown report generated: {output_file}")


def run_audit(base_dir: str = None, output_json: str = None, output_md: str = None) -> Dict[str, Any]:
    """
    Run the security audit and generate reports.
    
    Args:
        base_dir: Base directory of the project to audit
        output_json: Path to the output JSON report file
        output_md: Path to the output Markdown report file
        
    Returns:
        Dictionary containing all findings and summary information
    """
    auditor = SecurityAuditor(base_dir)
    report = auditor.run_all_audits()
    
    if output_json:
        auditor.generate_json_report(output_json)
    
    if output_md:
        auditor.generate_markdown_report(output_md)
    
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="CryoProtect v2 Security Auditor")
    parser.add_argument("--dir", help="Base directory to audit", default=os.getcwd())
    parser.add_argument("--output-json", help="Output file for JSON report")
    parser.add_argument("--output-md", help="Output file for Markdown report")
    args = parser.parse_args()
    
    if not args.output_json and not args.output_md:
        parser.error("At least one of --output-json or --output-md must be specified")
    
    run_audit(args.dir, args.output_json, args.output_md)