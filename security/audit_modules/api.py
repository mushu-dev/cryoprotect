"""
API Security Audit Module

This module provides functionality to audit API security configurations and vulnerabilities
in the CryoProtect v2 system, focusing on OWASP Top 10 vulnerabilities related to APIs, such as:
- A1:2017 Injection
- A3:2017 Sensitive Data Exposure
- A4:2017 XML External Entities (XXE)
- A6:2017 Security Misconfiguration
- A7:2017 Cross-Site Scripting (XSS)
"""

import os
import re
import json
import logging
from enum import Enum
from typing import Dict, List, Any, Optional, Tuple
from pathlib import Path
import configparser

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

class SeverityLevel(Enum):
    """Severity levels for security findings"""
    CRITICAL = 5
    HIGH = 4
    MEDIUM = 3
    LOW = 2
    INFO = 1
    
    def __str__(self):
        return self.name

class ApiAuditor:
    """
    Auditor for API security configurations and vulnerabilities.
    Checks for common security vulnerabilities related to APIs.
    """
    
    def __init__(self, base_dir: str = None):
        """
        Initialize the API auditor.
        
        Args:
            base_dir: Base directory of the project to audit
        """
        self.base_dir = base_dir or os.getcwd()
        self.findings = []
        logger.info(f"Initialized ApiAuditor with base directory: {self.base_dir}")
    
    def run_all_checks(self) -> List[Dict[str, Any]]:
        """
        Run all API security checks.
        
        Returns:
            List of findings from all checks
        """
        logger.info("Running all API security checks")
        
        # Clear previous findings
        self.findings = []
        
        # Run all checks
        self.check_sql_injection()
        self.check_nosql_injection()
        self.check_command_injection()
        self.check_xxe_vulnerabilities()
        self.check_cors_configuration()
        self.check_rate_limiting()
        self.check_api_documentation()
        self.check_input_validation()
        self.check_error_handling()
        self.check_api_versioning()
        self.check_content_type()
        
        logger.info(f"Completed API security checks. Found {len(self.findings)} issues.")
        return self.findings
    
    def check_sql_injection(self) -> None:
        """Check for SQL injection vulnerabilities"""
        logger.info("Checking for SQL injection vulnerabilities")
        
        # Look for database query code
        sql_files = self._find_files_with_pattern(r'execute.*sql|cursor\.execute|db\.query|query\(|raw.*sql')
        
        if not sql_files:
            logger.info("No SQL query code found to check for SQL injection")
            return
        
        # Check for SQL injection in found files
        for file_path in sql_files:
            with open(file_path, 'r', encoding='utf-8', errors='ignore') as f:
                content = f.read()
                
                # Check for string concatenation in SQL queries
                if re.search(r'execute\(.*\+.*\)|execute\(.*%.*\)|execute\(.*\.format\(|execute\(.*f[\'"]', content):
                    self._add_finding(
                        title="Potential SQL Injection",
                        description=f"SQL queries in {file_path} may be vulnerable to SQL injection due to string concatenation or formatting.",
                        severity=SeverityLevel.CRITICAL,
                        category="Injection",
                        location=file_path,
                        recommendation="Use parameterized queries or prepared statements instead of string concatenation or formatting."
                    )
                
                # Check for lack of parameterized queries
                if re.search(r'execute\(', content) and not re.search(r'execute\(.*\%s.*\)|execute\(.*\?.*\)|execute\(.*:.*\)', content):
                    self._add_finding(
                        title="Lack of parameterized queries",
                        description=f"SQL queries in {file_path} may not be using parameterized queries, which increases the risk of SQL injection.",
                        severity=SeverityLevel.HIGH,
                        category="Injection",
                        location=file_path,
                        recommendation="Use parameterized queries or prepared statements with placeholders for all user-supplied input."
                    )
    
    def check_nosql_injection(self) -> None:
        """Check for NoSQL injection vulnerabilities"""
        logger.info("Checking for NoSQL injection vulnerabilities")
        
        # Look for NoSQL database code
        nosql_files = self._find_files_with_pattern(r'mongodb|mongoose|firestore|dynamodb|cosmosdb')
        
        if not nosql_files:
            logger.info("No NoSQL database code found to check for NoSQL injection")
            return
        
        # Check for NoSQL injection in found files
        for file_path in nosql_files:
            with open(file_path, 'r', encoding='utf-8', errors='ignore') as f:
                content = f.read()
                
                # Check for string concatenation in NoSQL queries
                if re.search(r'find\(.*\+.*\)|find\(.*%.*\)|find\(.*\.format\(|find\(.*f[\'"]', content):
                    self._add_finding(
                        title="Potential NoSQL Injection",
                        description=f"NoSQL queries in {file_path} may be vulnerable to NoSQL injection due to string concatenation or formatting.",
                        severity=SeverityLevel.CRITICAL,
                        category="Injection",
                        location=file_path,
                        recommendation="Use parameterized queries and proper input validation for NoSQL databases."
                    )
                
                # Check for user input directly in query objects
                if re.search(r'find\(\s*{.*request\..*}', content):
                    self._add_finding(
                        title="User input in NoSQL query objects",
                        description=f"NoSQL queries in {file_path} may be using user input directly in query objects, which increases the risk of NoSQL injection.",
                        severity=SeverityLevel.HIGH,
                        category="Injection",
                        location=file_path,
                        recommendation="Validate and sanitize user input before using it in NoSQL query objects."
                    )
    
    def check_command_injection(self) -> None:
        """Check for command injection vulnerabilities"""
        logger.info("Checking for command injection vulnerabilities")
        
        # Look for command execution code
        cmd_files = self._find_files_with_pattern(r'subprocess|os\.system|os\.popen|exec|eval|child_process|spawn|shell_exec')
        
        if not cmd_files:
            logger.info("No command execution code found to check for command injection")
            return
        
        # Check for command injection in found files
        for file_path in cmd_files:
            with open(file_path, 'r', encoding='utf-8', errors='ignore') as f:
                content = f.read()
                
                # Check for string concatenation in command execution
                if re.search(r'(subprocess\.call|os\.system|os\.popen|exec|eval)\(.*\+.*\)', content) or \
                   re.search(r'(subprocess\.call|os\.system|os\.popen|exec|eval)\(.*%.*\)', content) or \
                   re.search(r'(subprocess\.call|os\.system|os\.popen|exec|eval)\(.*\.format\(', content) or \
                   re.search(r'(subprocess\.call|os\.system|os\.popen|exec|eval)\(.*f[\'"]', content):
                    self._add_finding(
                        title="Potential Command Injection",
                        description=f"Command execution in {file_path} may be vulnerable to command injection due to string concatenation or formatting.",
                        severity=SeverityLevel.CRITICAL,
                        category="Injection",
                        location=file_path,
                        recommendation="Use subprocess module with shell=False and pass arguments as a list instead of a string."
                    )
                
                # Check for shell=True in subprocess calls
                if re.search(r'subprocess\..*shell\s*=\s*True', content):
                    self._add_finding(
                        title="Unsafe subprocess shell=True",
                        description=f"Subprocess calls in {file_path} use shell=True, which increases the risk of command injection.",
                        severity=SeverityLevel.HIGH,
                        category="Injection",
                        location=file_path,
                        recommendation="Avoid using shell=True in subprocess calls. Pass arguments as a list instead of a string."
                    )
    
    def check_xxe_vulnerabilities(self) -> None:
        """Check for XML External Entity (XXE) vulnerabilities"""
        logger.info("Checking for XXE vulnerabilities")
        
        # Look for XML parsing code
        xml_files = self._find_files_with_pattern(r'xml|ElementTree|minidom|lxml|etree|parseXml|DOMParser')
        
        if not xml_files:
            logger.info("No XML parsing code found to check for XXE vulnerabilities")
            return
        
        # Check for XXE vulnerabilities in found files
        for file_path in xml_files:
            with open(file_path, 'r', encoding='utf-8', errors='ignore') as f:
                content = f.read()
                
                # Check for disabling of external entity processing
                if re.search(r'xml', content, re.IGNORECASE) and not re.search(r'resolve_entities\s*=\s*False|no_network\s*=\s*True|FEATURE.*external-general-entities.*false', content, re.IGNORECASE):
                    self._add_finding(
                        title="Potential XXE Vulnerability",
                        description=f"XML parsing in {file_path} may be vulnerable to XXE attacks due to not disabling external entity processing.",
                        severity=SeverityLevel.HIGH,
                        category="XML External Entities (XXE)",
                        location=file_path,
                        recommendation="Disable external entity processing in XML parsers by setting appropriate parser options."
                    )
    
    def check_cors_configuration(self) -> None:
        """Check for insecure CORS configuration"""
        logger.info("Checking CORS configuration")
        
        # Look for CORS configuration code
        cors_files = self._find_files_with_pattern(r'cors|Access-Control-Allow|CORS_ORIGIN')
        
        if not cors_files:
            self._add_finding(
                title="No CORS configuration found",
                description="Could not find evidence of CORS configuration in the codebase.",
                severity=SeverityLevel.MEDIUM,
                category="Security Misconfiguration",
                location="N/A",
                recommendation="Implement proper CORS configuration to restrict cross-origin requests to trusted domains."
            )
            return
        
        # Check CORS configuration in found files
        for file_path in cors_files:
            with open(file_path, 'r', encoding='utf-8', errors='ignore') as f:
                content = f.read()
                
                # Check for wildcard origin
                if re.search(r'Access-Control-Allow-Origin\s*:\s*\*|origin\s*:\s*\*|origins\s*=\s*\*', content, re.IGNORECASE):
                    self._add_finding(
                        title="Wildcard CORS origin",
                        description=f"CORS configuration in {file_path} allows requests from any origin using a wildcard (*).",
                        severity=SeverityLevel.HIGH,
                        category="Security Misconfiguration",
                        location=file_path,
                        recommendation="Restrict CORS to specific trusted domains instead of using a wildcard."
                    )
                
                # Check for overly permissive origins
                if re.search(r'origin.*null|null.*origin', content, re.IGNORECASE):
                    self._add_finding(
                        title="CORS allows null origin",
                        description=f"CORS configuration in {file_path} allows requests from 'null' origin, which can be spoofed.",
                        severity=SeverityLevel.MEDIUM,
                        category="Security Misconfiguration",
                        location=file_path,
                        recommendation="Remove 'null' from allowed origins in CORS configuration."
                    )
    
    def check_rate_limiting(self) -> None:
        """Check for API rate limiting"""
        logger.info("Checking for API rate limiting")
        
        # Look for rate limiting code
        rate_files = self._find_files_with_pattern(r'rate.*limit|throttle|limiter')
        
        if not rate_files:
            self._add_finding(
                title="No API rate limiting found",
                description="Could not find evidence of API rate limiting in the codebase.",
                severity=SeverityLevel.MEDIUM,
                category="Security Misconfiguration",
                location="N/A",
                recommendation="Implement API rate limiting to prevent abuse and denial of service attacks."
            )
            return
        
        # If rate limiting is implemented, check for common issues
        for file_path in rate_files:
            with open(file_path, 'r', encoding='utf-8', errors='ignore') as f:
                content = f.read()
                
                # Check for high rate limits
                rate_match = re.search(r'rate.*?(\d+).*?(?:minute|min|second|sec)', content, re.IGNORECASE)
                if rate_match and int(rate_match.group(1)) > 100:
                    self._add_finding(
                        title="High API rate limit",
                        description=f"API rate limit in {file_path} is set to {rate_match.group(1)} requests, which may be too high.",
                        severity=SeverityLevel.LOW,
                        category="Security Misconfiguration",
                        location=file_path,
                        recommendation="Consider reducing the API rate limit to prevent abuse."
                    )
    
    def check_api_documentation(self) -> None:
        """Check for API documentation"""
        logger.info("Checking for API documentation")
        
        # Look for API documentation
        doc_files = self._find_files_with_pattern(r'swagger|openapi|api.*doc|apidoc')
        
        if not doc_files:
            self._add_finding(
                title="No API documentation found",
                description="Could not find evidence of API documentation in the codebase.",
                severity=SeverityLevel.LOW,
                category="Documentation",
                location="N/A",
                recommendation="Implement API documentation using tools like Swagger or OpenAPI to improve security awareness and testing."
            )
            return
    
    def check_input_validation(self) -> None:
        """Check for input validation"""
        logger.info("Checking for input validation")
        
        # Look for input validation code
        validation_files = self._find_files_with_pattern(r'validate|validation|sanitize|sanitization|escape|filter')
        
        if not validation_files:
            self._add_finding(
                title="No input validation found",
                description="Could not find evidence of input validation in the codebase.",
                severity=SeverityLevel.HIGH,
                category="Injection",
                location="N/A",
                recommendation="Implement proper input validation for all user-supplied data to prevent injection attacks."
            )
            return
        
        # Check input validation in found files
        for file_path in validation_files:
            with open(file_path, 'r', encoding='utf-8', errors='ignore') as f:
                content = f.read()
                
                # Check for client-side only validation
                if 'javascript' in file_path.lower() or 'js' in file_path.lower() or '.ts' in file_path.lower():
                    if re.search(r'validate|validation', content, re.IGNORECASE) and not self._find_server_validation():
                        self._add_finding(
                            title="Client-side only validation",
                            description=f"Input validation in {file_path} appears to be implemented only on the client side.",
                            severity=SeverityLevel.MEDIUM,
                            category="Injection",
                            location=file_path,
                            recommendation="Implement server-side validation in addition to client-side validation."
                        )
    
    def _find_server_validation(self) -> bool:
        """
        Check if server-side validation exists.
        
        Returns:
            True if server-side validation is found, False otherwise
        """
        server_files = []
        for root, _, files in os.walk(self.base_dir):
            for file in files:
                if file.endswith(('.py', '.rb', '.php', '.java', '.cs')):
                    file_path = os.path.join(root, file)
                    server_files.append(file_path)
        
        for file_path in server_files:
            try:
                with open(file_path, 'r', encoding='utf-8', errors='ignore') as f:
                    content = f.read()
                    if re.search(r'validate|validation|sanitize|sanitization', content, re.IGNORECASE):
                        return True
            except Exception as e:
                logger.warning(f"Error reading file {file_path}: {e}")
        
        return False
    
    def check_error_handling(self) -> None:
        """Check for proper error handling"""
        logger.info("Checking for proper error handling")
        
        # Look for error handling code
        error_files = self._find_files_with_pattern(r'try.*catch|except|error.*handler|handle.*error|error.*response')
        
        if not error_files:
            self._add_finding(
                title="No error handling found",
                description="Could not find evidence of proper error handling in the codebase.",
                severity=SeverityLevel.MEDIUM,
                category="Security Misconfiguration",
                location="N/A",
                recommendation="Implement proper error handling to prevent information leakage and improve security."
            )
            return
        
        # Check error handling in found files
        for file_path in error_files:
            with open(file_path, 'r', encoding='utf-8', errors='ignore') as f:
                content = f.read()
                
                # Check for sensitive information in error messages
                if re.search(r'(error|exception).*message.*\+|stack.*trace', content, re.IGNORECASE):
                    self._add_finding(
                        title="Potential sensitive information in error messages",
                        description=f"Error handling in {file_path} may expose sensitive information such as stack traces or detailed error messages.",
                        severity=SeverityLevel.MEDIUM,
                        category="Sensitive Data Exposure",
                        location=file_path,
                        recommendation="Implement proper error handling that does not expose sensitive information to users."
                    )
    
    def check_api_versioning(self) -> None:
        """Check for API versioning"""
        logger.info("Checking for API versioning")
        
        # Look for API versioning code
        version_files = self._find_files_with_pattern(r'api.*v\d|v\d.*api|version')
        
        if not version_files:
            self._add_finding(
                title="No API versioning found",
                description="Could not find evidence of API versioning in the codebase.",
                severity=SeverityLevel.LOW,
                category="API Design",
                location="N/A",
                recommendation="Implement API versioning to ensure backward compatibility and security updates."
            )
            return
    
    def check_content_type(self) -> None:
        """Check for proper Content-Type headers"""
        logger.info("Checking for proper Content-Type headers")
        
        # Look for Content-Type headers
        content_type_files = self._find_files_with_pattern(r'Content-Type|content.*type|MIME')
        
        if not content_type_files:
            self._add_finding(
                title="No Content-Type headers found",
                description="Could not find evidence of Content-Type headers in the codebase.",
                severity=SeverityLevel.LOW,
                category="Security Misconfiguration",
                location="N/A",
                recommendation="Set proper Content-Type headers for all API responses to prevent content sniffing attacks."
            )
            return
        
        # Check Content-Type headers in found files
        for file_path in content_type_files:
            with open(file_path, 'r', encoding='utf-8', errors='ignore') as f:
                content = f.read()
                
                # Check for X-Content-Type-Options header
                if 'Content-Type' in content and not re.search(r'X-Content-Type-Options\s*:\s*nosniff', content, re.IGNORECASE):
                    self._add_finding(
                        title="Missing X-Content-Type-Options header",
                        description=f"Content-Type headers in {file_path} are not accompanied by X-Content-Type-Options: nosniff header.",
                        severity=SeverityLevel.LOW,
                        category="Security Misconfiguration",
                        location=file_path,
                        recommendation="Add X-Content-Type-Options: nosniff header to prevent content sniffing attacks."
                    )
    
    def _find_files_with_pattern(self, pattern: str) -> List[str]:
        """
        Find files containing the given pattern.
        
        Args:
            pattern: Regular expression pattern to search for
            
        Returns:
            List of file paths containing the pattern
        """
        result = []
        pattern_re = re.compile(pattern, re.IGNORECASE)
        
        for root, _, files in os.walk(self.base_dir):
            for file in files:
                if file.endswith(('.py', '.js', '.ts', '.html', '.config', '.json', '.xml')):
                    file_path = os.path.join(root, file)
                    try:
                        with open(file_path, 'r', encoding='utf-8', errors='ignore') as f:
                            content = f.read()
                            if pattern_re.search(content):
                                result.append(file_path)
                    except Exception as e:
                        logger.warning(f"Error reading file {file_path}: {e}")
        
        return result
    
    def _add_finding(self, title: str, description: str, severity: SeverityLevel,
                    category: str, location: str, recommendation: str) -> None:
        """
        Add a security finding to the list.
        
        Args:
            title: Short title of the finding
            description: Detailed description of the finding
            severity: Severity level of the finding
            category: Category of the finding (e.g., "Injection")
            location: File or component where the finding was detected
            recommendation: Recommended action to address the finding
        """
        finding = {
            "title": title,
            "description": description,
            "severity": severity.name,
            "severity_value": severity.value,
            "category": category,
            "location": location,
            "recommendation": recommendation,
            "module": "api"
        }
        
        self.findings.append(finding)
        logger.info(f"Added finding: {title} ({severity.name})")


def run_audit(base_dir: str = None) -> List[Dict[str, Any]]:
    """
    Run the API security audit.
    
    Args:
        base_dir: Base directory of the project to audit
        
    Returns:
        List of findings from the audit
    """
    auditor = ApiAuditor(base_dir)
    return auditor.run_all_checks()


if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="API Security Audit")
    parser.add_argument("--dir", help="Base directory to audit", default=os.getcwd())
    parser.add_argument("--output", help="Output file for findings (JSON format)")
    args = parser.parse_args()
    
    findings = run_audit(args.dir)
    
    if args.output:
        with open(args.output, 'w') as f:
            json.dump(findings, f, indent=2)
        print(f"Findings written to {args.output}")
    else:
        print(json.dumps(findings, indent=2))