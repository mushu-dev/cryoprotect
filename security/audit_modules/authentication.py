"""
Authentication Security Audit Module

This module provides functionality to audit authentication and authorization mechanisms
in the CryoProtect v2 system, focusing on OWASP Top 10 vulnerabilities related to
authentication, such as:
- A2:2017 Broken Authentication
- A5:2017 Broken Access Control
- A6:2017 Security Misconfiguration
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

class AuthenticationAuditor:
    """
    Auditor for authentication and authorization mechanisms.
    Checks for common security vulnerabilities related to authentication.
    """
    
    def __init__(self, base_dir: str = None):
        """
        Initialize the authentication auditor.
        
        Args:
            base_dir: Base directory of the project to audit
        """
        self.base_dir = base_dir or os.getcwd()
        self.findings = []
        logger.info(f"Initialized AuthenticationAuditor with base directory: {self.base_dir}")
    
    def run_all_checks(self) -> List[Dict[str, Any]]:
        """
        Run all authentication security checks.
        
        Returns:
            List of findings from all checks
        """
        logger.info("Running all authentication security checks")
        
        # Clear previous findings
        self.findings = []
        
        # Run all checks
        self.check_password_policies()
        self.check_session_management()
        self.check_access_control()
        self.check_auth_logging()
        self.check_mfa_implementation()
        self.check_password_storage()
        self.check_account_lockout()
        self.check_csrf_protection()
        
        logger.info(f"Completed authentication security checks. Found {len(self.findings)} issues.")
        return self.findings
    
    def check_password_policies(self) -> None:
        """Check for weak password policies"""
        logger.info("Checking password policies")
        
        # Look for password validation code
        password_files = self._find_files_with_pattern(r'password.*validation|validate.*password')
        
        if not password_files:
            self._add_finding(
                title="No password validation found",
                description="Could not find evidence of password validation in the codebase.",
                severity=SeverityLevel.HIGH,
                category="Broken Authentication",
                location="N/A",
                recommendation="Implement strong password validation that enforces minimum length, "
                              "complexity requirements, and checks against common passwords."
            )
            return
        
        # Check password validation in found files
        for file_path in password_files:
            with open(file_path, 'r', encoding='utf-8', errors='ignore') as f:
                content = f.read()
                
                # Check for minimum password length
                length_check = re.search(r'(min.*length|length.*min).*?(\d+)', content, re.IGNORECASE)
                if not length_check or int(length_check.group(2)) < 8:
                    self._add_finding(
                        title="Weak password length requirement",
                        description=f"Password validation in {file_path} does not enforce a minimum length of 8 characters.",
                        severity=SeverityLevel.MEDIUM,
                        category="Broken Authentication",
                        location=file_path,
                        recommendation="Enforce a minimum password length of at least 8 characters, "
                                      "preferably 12 or more for sensitive systems."
                    )
                
                # Check for complexity requirements
                if not re.search(r'(upper.*lower|lower.*upper)', content, re.IGNORECASE):
                    self._add_finding(
                        title="Insufficient password complexity",
                        description=f"Password validation in {file_path} does not enforce the use of both uppercase and lowercase letters.",
                        severity=SeverityLevel.MEDIUM,
                        category="Broken Authentication",
                        location=file_path,
                        recommendation="Enforce password complexity requirements including uppercase, "
                                      "lowercase, numbers, and special characters."
                    )
                
                # Check for common password check
                if not re.search(r'common.*password|dictionary.*attack', content, re.IGNORECASE):
                    self._add_finding(
                        title="No check against common passwords",
                        description=f"Password validation in {file_path} does not check against commonly used or breached passwords.",
                        severity=SeverityLevel.MEDIUM,
                        category="Broken Authentication",
                        location=file_path,
                        recommendation="Implement checks against lists of commonly used or breached passwords."
                    )
    
    def check_session_management(self) -> None:
        """Check for secure session management"""
        logger.info("Checking session management")
        
        # Look for session management code
        session_files = self._find_files_with_pattern(r'session|cookie|jwt|token')
        
        if not session_files:
            self._add_finding(
                title="No session management found",
                description="Could not find evidence of session management in the codebase.",
                severity=SeverityLevel.MEDIUM,
                category="Broken Authentication",
                location="N/A",
                recommendation="Implement secure session management with proper timeout, "
                              "secure cookie flags, and token validation."
            )
            return
        
        # Check session management in found files
        for file_path in session_files:
            with open(file_path, 'r', encoding='utf-8', errors='ignore') as f:
                content = f.read()
                
                # Check for secure cookie flag
                if 'cookie' in content.lower() and not re.search(r'secure\s*[=:]\s*true', content, re.IGNORECASE):
                    self._add_finding(
                        title="Missing Secure flag on cookies",
                        description=f"Cookies in {file_path} are not set with the Secure flag, allowing transmission over unencrypted connections.",
                        severity=SeverityLevel.HIGH,
                        category="Broken Authentication",
                        location=file_path,
                        recommendation="Set the Secure flag on all sensitive cookies to ensure they are only transmitted over HTTPS."
                    )
                
                # Check for HttpOnly flag
                if 'cookie' in content.lower() and not re.search(r'httponly\s*[=:]\s*true', content, re.IGNORECASE):
                    self._add_finding(
                        title="Missing HttpOnly flag on cookies",
                        description=f"Cookies in {file_path} are not set with the HttpOnly flag, making them accessible to client-side scripts.",
                        severity=SeverityLevel.HIGH,
                        category="Broken Authentication",
                        location=file_path,
                        recommendation="Set the HttpOnly flag on all sensitive cookies to prevent access from client-side scripts."
                    )
                
                # Check for SameSite attribute
                if 'cookie' in content.lower() and not re.search(r'samesite\s*[=:]\s*(strict|lax)', content, re.IGNORECASE):
                    self._add_finding(
                        title="Missing SameSite attribute on cookies",
                        description=f"Cookies in {file_path} are not set with the SameSite attribute, making them vulnerable to CSRF attacks.",
                        severity=SeverityLevel.MEDIUM,
                        category="Broken Authentication",
                        location=file_path,
                        recommendation="Set the SameSite attribute to 'Strict' or 'Lax' on all cookies to prevent CSRF attacks."
                    )
                
                # Check for session timeout
                if 'session' in content.lower() and not re.search(r'(timeout|expire|ttl|maxage)', content, re.IGNORECASE):
                    self._add_finding(
                        title="No session timeout",
                        description=f"Session management in {file_path} does not implement session timeout.",
                        severity=SeverityLevel.MEDIUM,
                        category="Broken Authentication",
                        location=file_path,
                        recommendation="Implement session timeout to automatically invalidate sessions after a period of inactivity."
                    )
    
    def check_access_control(self) -> None:
        """Check for proper access control mechanisms"""
        logger.info("Checking access control mechanisms")
        
        # Look for access control code
        access_files = self._find_files_with_pattern(r'auth|permission|role|access|rbac')
        
        if not access_files:
            self._add_finding(
                title="No access control found",
                description="Could not find evidence of access control mechanisms in the codebase.",
                severity=SeverityLevel.HIGH,
                category="Broken Access Control",
                location="N/A",
                recommendation="Implement proper access control mechanisms such as RBAC or ABAC "
                              "to restrict access to resources based on user roles and permissions."
            )
            return
        
        # Check access control in found files
        for file_path in access_files:
            with open(file_path, 'r', encoding='utf-8', errors='ignore') as f:
                content = f.read()
                
                # Check for hardcoded roles or permissions
                if re.search(r'role\s*[=:]\s*[\'"]admin[\'"]', content, re.IGNORECASE):
                    self._add_finding(
                        title="Hardcoded roles",
                        description=f"Hardcoded roles found in {file_path}, which may lead to inflexible access control.",
                        severity=SeverityLevel.LOW,
                        category="Broken Access Control",
                        location=file_path,
                        recommendation="Use a configurable role-based access control system instead of hardcoded roles."
                    )
                
                # Check for insecure direct object references
                if re.search(r'id\s*=\s*request\.(get|post|param)', content, re.IGNORECASE):
                    self._add_finding(
                        title="Potential Insecure Direct Object Reference (IDOR)",
                        description=f"Potential IDOR vulnerability found in {file_path}, where object IDs are directly used from user input without proper authorization checks.",
                        severity=SeverityLevel.HIGH,
                        category="Broken Access Control",
                        location=file_path,
                        recommendation="Implement proper authorization checks before accessing objects based on user-supplied IDs."
                    )
    
    def check_auth_logging(self) -> None:
        """Check for authentication event logging"""
        logger.info("Checking authentication event logging")
        
        # Look for authentication logging code
        log_files = self._find_files_with_pattern(r'log.*auth|auth.*log|login.*log|log.*login')
        
        if not log_files:
            self._add_finding(
                title="No authentication logging found",
                description="Could not find evidence of authentication event logging in the codebase.",
                severity=SeverityLevel.MEDIUM,
                category="Insufficient Logging & Monitoring",
                location="N/A",
                recommendation="Implement logging for all authentication events including successes, "
                              "failures, password changes, and privilege changes."
            )
            return
        
        # Check logging in found files
        for file_path in log_files:
            with open(file_path, 'r', encoding='utf-8', errors='ignore') as f:
                content = f.read()
                
                # Check for login failure logging
                if not re.search(r'log.*fail|fail.*log', content, re.IGNORECASE):
                    self._add_finding(
                        title="No login failure logging",
                        description=f"Authentication logging in {file_path} does not appear to log login failures.",
                        severity=SeverityLevel.MEDIUM,
                        category="Insufficient Logging & Monitoring",
                        location=file_path,
                        recommendation="Implement logging for authentication failures to detect potential brute force attacks."
                    )
                
                # Check for sensitive data in logs
                if re.search(r'log.*password|password.*log', content, re.IGNORECASE):
                    self._add_finding(
                        title="Potential sensitive data in logs",
                        description=f"Authentication logging in {file_path} may be logging sensitive data such as passwords.",
                        severity=SeverityLevel.HIGH,
                        category="Sensitive Data Exposure",
                        location=file_path,
                        recommendation="Ensure that sensitive data such as passwords are never logged, even in hashed form."
                    )
    
    def check_mfa_implementation(self) -> None:
        """Check for multi-factor authentication implementation"""
        logger.info("Checking multi-factor authentication")
        
        # Look for MFA code
        mfa_files = self._find_files_with_pattern(r'mfa|2fa|two.*factor|multi.*factor|otp|totp')
        
        if not mfa_files:
            self._add_finding(
                title="No multi-factor authentication found",
                description="Could not find evidence of multi-factor authentication in the codebase.",
                severity=SeverityLevel.MEDIUM,
                category="Broken Authentication",
                location="N/A",
                recommendation="Implement multi-factor authentication for all user accounts, "
                              "especially for administrative or privileged accounts."
            )
            return
        
        # If MFA is implemented, check for common issues
        for file_path in mfa_files:
            with open(file_path, 'r', encoding='utf-8', errors='ignore') as f:
                content = f.read()
                
                # Check for SMS-based MFA (less secure)
                if re.search(r'sms|text.*message', content, re.IGNORECASE):
                    self._add_finding(
                        title="SMS-based MFA",
                        description=f"SMS-based MFA found in {file_path}, which is vulnerable to SIM swapping attacks.",
                        severity=SeverityLevel.LOW,
                        category="Broken Authentication",
                        location=file_path,
                        recommendation="Consider using more secure MFA methods such as TOTP apps or hardware tokens instead of SMS."
                    )
    
    def check_password_storage(self) -> None:
        """Check for secure password storage"""
        logger.info("Checking password storage")
        
        # Look for password storage code
        password_files = self._find_files_with_pattern(r'password.*hash|hash.*password|encrypt.*password|password.*encrypt|bcrypt|scrypt|pbkdf2')
        
        if not password_files:
            self._add_finding(
                title="No secure password storage found",
                description="Could not find evidence of secure password hashing in the codebase.",
                severity=SeverityLevel.CRITICAL,
                category="Sensitive Data Exposure",
                location="N/A",
                recommendation="Implement secure password storage using modern hashing algorithms "
                              "such as bcrypt, Argon2, or PBKDF2 with sufficient work factors."
            )
            return
        
        # Check password storage in found files
        for file_path in password_files:
            with open(file_path, 'r', encoding='utf-8', errors='ignore') as f:
                content = f.read()
                
                # Check for weak hashing algorithms
                if re.search(r'md5|sha1', content, re.IGNORECASE):
                    self._add_finding(
                        title="Weak password hashing algorithm",
                        description=f"Weak password hashing algorithm (MD5 or SHA1) found in {file_path}.",
                        severity=SeverityLevel.CRITICAL,
                        category="Sensitive Data Exposure",
                        location=file_path,
                        recommendation="Replace weak hashing algorithms with modern algorithms "
                                      "such as bcrypt, Argon2, or PBKDF2 with sufficient work factors."
                    )
                
                # Check for missing salt
                if not re.search(r'salt', content, re.IGNORECASE) and not re.search(r'bcrypt|argon2', content, re.IGNORECASE):
                    self._add_finding(
                        title="Missing salt in password hashing",
                        description=f"No evidence of salt usage in password hashing found in {file_path}.",
                        severity=SeverityLevel.HIGH,
                        category="Sensitive Data Exposure",
                        location=file_path,
                        recommendation="Use a modern hashing algorithm that automatically handles salting, "
                                      "such as bcrypt or Argon2."
                    )
    
    def check_account_lockout(self) -> None:
        """Check for account lockout mechanisms"""
        logger.info("Checking account lockout mechanisms")
        
        # Look for account lockout code
        lockout_files = self._find_files_with_pattern(r'lockout|bruteforce|brute.*force|rate.*limit|attempt.*limit|max.*attempt')
        
        if not lockout_files:
            self._add_finding(
                title="No account lockout mechanism found",
                description="Could not find evidence of account lockout or rate limiting mechanisms in the codebase.",
                severity=SeverityLevel.HIGH,
                category="Broken Authentication",
                location="N/A",
                recommendation="Implement account lockout or rate limiting mechanisms to prevent "
                              "brute force attacks on authentication endpoints."
            )
            return
        
        # If lockout mechanisms are found, check for common issues
        for file_path in lockout_files:
            with open(file_path, 'r', encoding='utf-8', errors='ignore') as f:
                content = f.read()
                
                # Check for too many allowed attempts
                attempts_match = re.search(r'max.*attempt.*?(\d+)|attempt.*limit.*?(\d+)', content, re.IGNORECASE)
                if attempts_match:
                    attempts = int(attempts_match.group(1) or attempts_match.group(2))
                    if attempts > 10:
                        self._add_finding(
                            title="High account lockout threshold",
                            description=f"Account lockout threshold in {file_path} is set to {attempts} attempts, which is too high.",
                            severity=SeverityLevel.MEDIUM,
                            category="Broken Authentication",
                            location=file_path,
                            recommendation="Reduce the account lockout threshold to 5 or fewer consecutive failed attempts."
                        )
    
    def check_csrf_protection(self) -> None:
        """Check for CSRF protection"""
        logger.info("Checking CSRF protection")
        
        # Look for CSRF protection code
        csrf_files = self._find_files_with_pattern(r'csrf|xsrf|cross.*site.*request|token')
        
        if not csrf_files:
            self._add_finding(
                title="No CSRF protection found",
                description="Could not find evidence of CSRF protection in the codebase.",
                severity=SeverityLevel.HIGH,
                category="Cross-Site Request Forgery",
                location="N/A",
                recommendation="Implement CSRF protection using tokens for all state-changing operations."
            )
            return
        
        # Check CSRF protection in found files
        for file_path in csrf_files:
            with open(file_path, 'r', encoding='utf-8', errors='ignore') as f:
                content = f.read()
                
                # Check for CSRF token validation
                if not re.search(r'validate.*token|verify.*token|check.*token', content, re.IGNORECASE):
                    self._add_finding(
                        title="Incomplete CSRF protection",
                        description=f"CSRF tokens may be generated but not validated in {file_path}.",
                        severity=SeverityLevel.MEDIUM,
                        category="Cross-Site Request Forgery",
                        location=file_path,
                        recommendation="Ensure that CSRF tokens are both generated and validated for all state-changing operations."
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
            category: Category of the finding (e.g., "Broken Authentication")
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
            "module": "authentication"
        }
        
        self.findings.append(finding)
        logger.info(f"Added finding: {title} ({severity.name})")


def run_audit(base_dir: str = None) -> List[Dict[str, Any]]:
    """
    Run the authentication security audit.
    
    Args:
        base_dir: Base directory of the project to audit
        
    Returns:
        List of findings from the audit
    """
    auditor = AuthenticationAuditor(base_dir)
    return auditor.run_all_checks()


if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="Authentication Security Audit")
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