"""
Data Protection Security Audit Module

This module provides functionality to audit data protection and privacy controls
in the CryoProtect v2 system, focusing on OWASP Top 10 vulnerabilities related to
data protection, such as:
- A3:2017 Sensitive Data Exposure
- A9:2017 Using Components with Known Vulnerabilities
- A10:2017 Insufficient Logging & Monitoring
"""

import os
import re
import json
import logging
from enum import Enum
from typing import Dict, List, Any, Optional, Tuple
from pathlib import Path
import configparser
import hashlib
import base64

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

class DataProtectionAuditor:
    """
    Auditor for data protection and privacy controls.
    Checks for common security vulnerabilities related to data protection.
    """
    
    def __init__(self, base_dir: str = None):
        """
        Initialize the data protection auditor.
        
        Args:
            base_dir: Base directory of the project to audit
        """
        self.base_dir = base_dir or os.getcwd()
        self.findings = []
        logger.info(f"Initialized DataProtectionAuditor with base directory: {self.base_dir}")
    
    def run_all_checks(self) -> List[Dict[str, Any]]:
        """
        Run all data protection security checks.
        
        Returns:
            List of findings from all checks
        """
        logger.info("Running all data protection security checks")
        
        # Clear previous findings
        self.findings = []
        
        # Run all checks
        self.check_sensitive_data_storage()
        self.check_data_encryption()
        self.check_pii_handling()
        self.check_logging_sensitive_data()
        self.check_data_retention()
        self.check_secure_file_uploads()
        self.check_database_backups()
        self.check_data_masking()
        self.check_dependency_vulnerabilities()
        self.check_data_validation()
        self.check_secure_headers()
        
        logger.info(f"Completed data protection security checks. Found {len(self.findings)} issues.")
        return self.findings
    
    def check_sensitive_data_storage(self) -> None:
        """Check for sensitive data storage issues"""
        logger.info("Checking sensitive data storage")
        
        # Look for sensitive data patterns
        sensitive_files = self._find_files_with_pattern(r'password|secret|key|token|credential|api.*key|private.*key')
        
        if not sensitive_files:
            logger.info("No files with sensitive data patterns found")
            return
        
        # Check for sensitive data in found files
        for file_path in sensitive_files:
            with open(file_path, 'r', encoding='utf-8', errors='ignore') as f:
                content = f.read()
                
                # Check for hardcoded credentials
                if re.search(r'password\s*=\s*[\'"][^\'"]+[\'"]|api.?key\s*=\s*[\'"][^\'"]+[\'"]|secret\s*=\s*[\'"][^\'"]+[\'"]', content, re.IGNORECASE):
                    self._add_finding(
                        title="Hardcoded credentials",
                        description=f"Hardcoded credentials found in {file_path}, which poses a security risk.",
                        severity=SeverityLevel.CRITICAL,
                        category="Sensitive Data Exposure",
                        location=file_path,
                        recommendation="Remove hardcoded credentials and use environment variables, "
                                      "secure credential storage, or a secrets management system."
                    )
                
                # Check for sensitive data in comments
                if re.search(r'#.*password.*=|//.*password.*=|/\*.*password.*=', content, re.IGNORECASE):
                    self._add_finding(
                        title="Sensitive data in comments",
                        description=f"Sensitive data found in comments in {file_path}, which poses a security risk.",
                        severity=SeverityLevel.MEDIUM,
                        category="Sensitive Data Exposure",
                        location=file_path,
                        recommendation="Remove sensitive data from comments in the codebase."
                    )
                
                # Check for sensitive data in configuration files
                if file_path.endswith(('.json', '.xml', '.config', '.ini', '.yml', '.yaml')) and \
                   re.search(r'password|secret|key|token|credential|api.*key', content, re.IGNORECASE):
                    self._add_finding(
                        title="Sensitive data in configuration files",
                        description=f"Sensitive data found in configuration file {file_path}, which poses a security risk.",
                        severity=SeverityLevel.HIGH,
                        category="Sensitive Data Exposure",
                        location=file_path,
                        recommendation="Remove sensitive data from configuration files and use environment variables "
                                      "or a secrets management system."
                    )
    
    def check_data_encryption(self) -> None:
        """Check for data encryption issues"""
        logger.info("Checking data encryption")
        
        # Look for encryption code
        encryption_files = self._find_files_with_pattern(r'encrypt|decrypt|cipher|aes|rsa|crypto')
        
        if not encryption_files:
            self._add_finding(
                title="No encryption implementation found",
                description="Could not find evidence of data encryption in the codebase.",
                severity=SeverityLevel.HIGH,
                category="Sensitive Data Exposure",
                location="N/A",
                recommendation="Implement data encryption for sensitive data at rest and in transit."
            )
            return
        
        # Check encryption in found files
        for file_path in encryption_files:
            with open(file_path, 'r', encoding='utf-8', errors='ignore') as f:
                content = f.read()
                
                # Check for weak encryption algorithms
                if re.search(r'des|rc4|md5|sha1', content, re.IGNORECASE):
                    self._add_finding(
                        title="Weak encryption algorithm",
                        description=f"Weak encryption algorithm (DES, RC4, MD5, or SHA1) found in {file_path}.",
                        severity=SeverityLevel.HIGH,
                        category="Sensitive Data Exposure",
                        location=file_path,
                        recommendation="Use strong encryption algorithms such as AES-256, RSA-2048, or SHA-256."
                    )
                
                # Check for ECB mode
                if re.search(r'ecb|electronic.*code.*book', content, re.IGNORECASE):
                    self._add_finding(
                        title="Insecure ECB mode",
                        description=f"Insecure ECB mode for encryption found in {file_path}.",
                        severity=SeverityLevel.HIGH,
                        category="Sensitive Data Exposure",
                        location=file_path,
                        recommendation="Use a more secure encryption mode such as CBC, GCM, or CTR."
                    )
                
                # Check for static initialization vectors
                if re.search(r'iv\s*=\s*[\'"][^\'"]+[\'"]', content, re.IGNORECASE):
                    self._add_finding(
                        title="Static initialization vector",
                        description=f"Static initialization vector for encryption found in {file_path}.",
                        severity=SeverityLevel.HIGH,
                        category="Sensitive Data Exposure",
                        location=file_path,
                        recommendation="Use a random initialization vector for each encryption operation."
                    )
    
    def check_pii_handling(self) -> None:
        """Check for personally identifiable information (PII) handling"""
        logger.info("Checking PII handling")
        
        # Look for PII patterns
        pii_files = self._find_files_with_pattern(r'name|email|address|phone|ssn|social.*security|credit.*card|passport|driver.*license|dob|birth.*date')
        
        if not pii_files:
            logger.info("No files with PII patterns found")
            return
        
        # Check for PII handling in found files
        for file_path in pii_files:
            with open(file_path, 'r', encoding='utf-8', errors='ignore') as f:
                content = f.read()
                
                # Check for PII in logs
                if re.search(r'log.*email|log.*phone|log.*ssn|log.*credit.*card', content, re.IGNORECASE):
                    self._add_finding(
                        title="PII in logs",
                        description=f"Personally identifiable information (PII) may be logged in {file_path}.",
                        severity=SeverityLevel.HIGH,
                        category="Sensitive Data Exposure",
                        location=file_path,
                        recommendation="Ensure that PII is not logged or is properly masked in logs."
                    )
                
                # Check for unencrypted PII storage
                if re.search(r'(email|phone|ssn|credit.*card).*=', content, re.IGNORECASE) and not re.search(r'encrypt', content, re.IGNORECASE):
                    self._add_finding(
                        title="Unencrypted PII storage",
                        description=f"Personally identifiable information (PII) may be stored unencrypted in {file_path}.",
                        severity=SeverityLevel.HIGH,
                        category="Sensitive Data Exposure",
                        location=file_path,
                        recommendation="Encrypt PII at rest and in transit."
                    )
    
    def check_logging_sensitive_data(self) -> None:
        """Check for sensitive data in logs"""
        logger.info("Checking for sensitive data in logs")
        
        # Look for logging code
        log_files = self._find_files_with_pattern(r'log\.|logger|logging|console\.log|print|debug')
        
        if not log_files:
            logger.info("No logging code found")
            return
        
        # Check for sensitive data in logs
        for file_path in log_files:
            with open(file_path, 'r', encoding='utf-8', errors='ignore') as f:
                content = f.read()
                
                # Check for sensitive data in logs
                if re.search(r'log.*password|log.*token|log.*key|log.*secret', content, re.IGNORECASE):
                    self._add_finding(
                        title="Sensitive data in logs",
                        description=f"Sensitive data may be logged in {file_path}.",
                        severity=SeverityLevel.HIGH,
                        category="Sensitive Data Exposure",
                        location=file_path,
                        recommendation="Ensure that sensitive data is not logged or is properly masked in logs."
                    )
                
                # Check for exception details in logs
                if re.search(r'log.*exception|log.*error|log.*stack.*trace', content, re.IGNORECASE):
                    self._add_finding(
                        title="Exception details in logs",
                        description=f"Exception details may be logged in {file_path}, which could expose sensitive information.",
                        severity=SeverityLevel.MEDIUM,
                        category="Sensitive Data Exposure",
                        location=file_path,
                        recommendation="Ensure that exception details are not logged in production or are properly sanitized."
                    )
    
    def check_data_retention(self) -> None:
        """Check for data retention policies"""
        logger.info("Checking data retention policies")
        
        # Look for data retention code
        retention_files = self._find_files_with_pattern(r'retention|expire|ttl|time.*live|delete.*data|purge.*data')
        
        if not retention_files:
            self._add_finding(
                title="No data retention policy found",
                description="Could not find evidence of data retention policies in the codebase.",
                severity=SeverityLevel.MEDIUM,
                category="Data Privacy",
                location="N/A",
                recommendation="Implement data retention policies to ensure that data is not kept longer than necessary."
            )
            return
    
    def check_secure_file_uploads(self) -> None:
        """Check for secure file upload handling"""
        logger.info("Checking secure file upload handling")
        
        # Look for file upload code
        upload_files = self._find_files_with_pattern(r'upload|file.*upload|multipart|form.*data')
        
        if not upload_files:
            logger.info("No file upload code found")
            return
        
        # Check for secure file upload handling
        for file_path in upload_files:
            with open(file_path, 'r', encoding='utf-8', errors='ignore') as f:
                content = f.read()
                
                # Check for file type validation
                if not re.search(r'content.*type|mime.*type|file.*type|extension', content, re.IGNORECASE):
                    self._add_finding(
                        title="No file type validation",
                        description=f"File upload handling in {file_path} may not validate file types, which could lead to security vulnerabilities.",
                        severity=SeverityLevel.HIGH,
                        category="Insecure File Upload",
                        location=file_path,
                        recommendation="Implement file type validation for all file uploads."
                    )
                
                # Check for file size validation
                if not re.search(r'size|length|maxsize|max.*size', content, re.IGNORECASE):
                    self._add_finding(
                        title="No file size validation",
                        description=f"File upload handling in {file_path} may not validate file sizes, which could lead to denial of service attacks.",
                        severity=SeverityLevel.MEDIUM,
                        category="Insecure File Upload",
                        location=file_path,
                        recommendation="Implement file size validation for all file uploads."
                    )
                
                # Check for file content validation
                if not re.search(r'validate.*content|scan|antivirus|malware', content, re.IGNORECASE):
                    self._add_finding(
                        title="No file content validation",
                        description=f"File upload handling in {file_path} may not validate file contents, which could lead to security vulnerabilities.",
                        severity=SeverityLevel.MEDIUM,
                        category="Insecure File Upload",
                        location=file_path,
                        recommendation="Implement file content validation for all file uploads."
                    )
    
    def check_database_backups(self) -> None:
        """Check for database backup security"""
        logger.info("Checking database backup security")
        
        # Look for database backup code
        backup_files = self._find_files_with_pattern(r'backup|dump|export|restore')
        
        if not backup_files:
            logger.info("No database backup code found")
            return
        
        # Check for database backup security
        for file_path in backup_files:
            with open(file_path, 'r', encoding='utf-8', errors='ignore') as f:
                content = f.read()
                
                # Check for encrypted backups
                if not re.search(r'encrypt|password|secure', content, re.IGNORECASE):
                    self._add_finding(
                        title="Unencrypted database backups",
                        description=f"Database backup handling in {file_path} may not encrypt backups, which could expose sensitive data.",
                        severity=SeverityLevel.HIGH,
                        category="Sensitive Data Exposure",
                        location=file_path,
                        recommendation="Encrypt all database backups."
                    )
                
                # Check for backup access controls
                if not re.search(r'access.*control|permission|restrict|acl', content, re.IGNORECASE):
                    self._add_finding(
                        title="No backup access controls",
                        description=f"Database backup handling in {file_path} may not implement access controls for backups, which could expose sensitive data.",
                        severity=SeverityLevel.MEDIUM,
                        category="Sensitive Data Exposure",
                        location=file_path,
                        recommendation="Implement access controls for all database backups."
                    )
    
    def check_data_masking(self) -> None:
        """Check for data masking implementation"""
        logger.info("Checking data masking implementation")
        
        # Look for data masking code
        masking_files = self._find_files_with_pattern(r'mask|redact|anonymize|pseudonymize|obfuscate')
        
        if not masking_files:
            self._add_finding(
                title="No data masking implementation found",
                description="Could not find evidence of data masking implementation in the codebase.",
                severity=SeverityLevel.MEDIUM,
                category="Sensitive Data Exposure",
                location="N/A",
                recommendation="Implement data masking for sensitive data in logs, error messages, and user interfaces."
            )
            return
    
    def check_dependency_vulnerabilities(self) -> None:
        """Check for vulnerable dependencies"""
        logger.info("Checking for vulnerable dependencies")
        
        # Look for dependency files
        dependency_files = []
        for root, _, files in os.walk(self.base_dir):
            for file in files:
                if file in ['requirements.txt', 'package.json', 'Gemfile', 'pom.xml', 'build.gradle']:
                    dependency_files.append(os.path.join(root, file))
        
        if not dependency_files:
            logger.info("No dependency files found")
            return
        
        # Check for dependency locking
        for file_path in dependency_files:
            with open(file_path, 'r', encoding='utf-8', errors='ignore') as f:
                content = f.read()
                
                # Check for version pinning
                if file_path.endswith('requirements.txt') and re.search(r'[a-zA-Z0-9_-]+[^=><~]$', content, re.MULTILINE):
                    self._add_finding(
                        title="Unpinned Python dependencies",
                        description=f"Python dependencies in {file_path} are not pinned to specific versions, which could lead to security vulnerabilities.",
                        severity=SeverityLevel.MEDIUM,
                        category="Using Components with Known Vulnerabilities",
                        location=file_path,
                        recommendation="Pin all Python dependencies to specific versions."
                    )
                
                # Check for package-lock.json
                if file_path.endswith('package.json') and not os.path.exists(os.path.join(os.path.dirname(file_path), 'package-lock.json')):
                    self._add_finding(
                        title="Missing package-lock.json",
                        description=f"No package-lock.json found for {file_path}, which could lead to inconsistent dependency versions.",
                        severity=SeverityLevel.LOW,
                        category="Using Components with Known Vulnerabilities",
                        location=file_path,
                        recommendation="Use package-lock.json to lock dependency versions."
                    )
    
    def check_data_validation(self) -> None:
        """Check for data validation"""
        logger.info("Checking data validation")
        
        # Look for data validation code
        validation_files = self._find_files_with_pattern(r'validate|validation|sanitize|sanitization|escape|filter')
        
        if not validation_files:
            self._add_finding(
                title="No data validation found",
                description="Could not find evidence of data validation in the codebase.",
                severity=SeverityLevel.HIGH,
                category="Injection",
                location="N/A",
                recommendation="Implement proper data validation for all user-supplied data to prevent injection attacks."
            )
            return
    
    def check_secure_headers(self) -> None:
        """Check for secure HTTP headers"""
        logger.info("Checking secure HTTP headers")
        
        # Look for HTTP header code
        header_files = self._find_files_with_pattern(r'header|response.*header|set.*header')
        
        if not header_files:
            logger.info("No HTTP header code found")
            return
        
        # Check for secure HTTP headers
        for file_path in header_files:
            with open(file_path, 'r', encoding='utf-8', errors='ignore') as f:
                content = f.read()
                
                # Check for Content-Security-Policy header
                if not re.search(r'Content-Security-Policy', content, re.IGNORECASE):
                    self._add_finding(
                        title="Missing Content-Security-Policy header",
                        description=f"HTTP header handling in {file_path} may not set the Content-Security-Policy header, which could lead to security vulnerabilities.",
                        severity=SeverityLevel.MEDIUM,
                        category="Security Misconfiguration",
                        location=file_path,
                        recommendation="Implement the Content-Security-Policy header to prevent XSS attacks."
                    )
                
                # Check for X-Content-Type-Options header
                if not re.search(r'X-Content-Type-Options', content, re.IGNORECASE):
                    self._add_finding(
                        title="Missing X-Content-Type-Options header",
                        description=f"HTTP header handling in {file_path} may not set the X-Content-Type-Options header, which could lead to security vulnerabilities.",
                        severity=SeverityLevel.LOW,
                        category="Security Misconfiguration",
                        location=file_path,
                        recommendation="Implement the X-Content-Type-Options header to prevent MIME type sniffing."
                    )
                
                # Check for X-Frame-Options header
                if not re.search(r'X-Frame-Options', content, re.IGNORECASE):
                    self._add_finding(
                        title="Missing X-Frame-Options header",
                        description=f"HTTP header handling in {file_path} may not set the X-Frame-Options header, which could lead to clickjacking attacks.",
                        severity=SeverityLevel.MEDIUM,
                        category="Security Misconfiguration",
                        location=file_path,
                        recommendation="Implement the X-Frame-Options header to prevent clickjacking attacks."
                    )
                
                # Check for Strict-Transport-Security header
                if not re.search(r'Strict-Transport-Security', content, re.IGNORECASE):
                    self._add_finding(
                        title="Missing Strict-Transport-Security header",
                        description=f"HTTP header handling in {file_path} may not set the Strict-Transport-Security header, which could lead to protocol downgrade attacks.",
                        severity=SeverityLevel.MEDIUM,
                        category="Security Misconfiguration",
                        location=file_path,
                        recommendation="Implement the Strict-Transport-Security header to enforce HTTPS."
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
            category: Category of the finding (e.g., "Sensitive Data Exposure")
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
            "module": "data_protection"
        }
        
        self.findings.append(finding)
        logger.info(f"Added finding: {title} ({severity.name})")


def run_audit(base_dir: str = None) -> List[Dict[str, Any]]:
    """
    Run the data protection security audit.
    
    Args:
        base_dir: Base directory of the project to audit
        
    Returns:
        List of findings from the audit
    """
    auditor = DataProtectionAuditor(base_dir)
    return auditor.run_all_checks()


if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="Data Protection Security Audit")
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