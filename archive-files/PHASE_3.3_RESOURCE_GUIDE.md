# Phase 3.3: Security Implementation Resource Guide

This guide provides resources, patterns, and examples for implementing the security phase of CryoProtect v2. It complements the Implementation Directive and Line References documents.

## Key Security Components

### 1. Security Audit Framework

The security audit should comprehensively check all aspects of the application:

- **Authentication**: JWT implementation, token validation, user session management
- **Authorization**: RLS policies, role-based access control, API permissions
- **Data Protection**: Encryption, secure storage, database security
- **API Security**: Input validation, rate limiting, CORS configuration
- **Infrastructure**: Docker security, environment configuration, secrets management
- **Web Security**: CSP, security headers, cookie protection

**Pattern for Security Auditor Implementation**:
- Create a modular audit framework with individual audit modules
- Implement scoring system for severity classification
- Use fixture data to test security controls
- Generate actionable reports with remediation steps

### 2. Vulnerability Scanning Tools

Key tools to integrate:

1. **Dependency Scanning**
   - Python: Safety, PyUp, Snyk
   - JavaScript: npm audit, Snyk
   - Docker: Trivy, Clair

2. **Code Scanning**
   - Python: Bandit
   - JavaScript: ESLint with security plugins
   - General: Semgrep, OWASP Dependency Check

3. **CI/CD Integration**
   - GitHub Action workflows for automated scanning
   - Fail build on critical vulnerabilities
   - Regular scheduled scans

### 3. Data Encryption Implementation

Effective encryption strategies:

1. **Field-Level Encryption**
   - Implement transparent encryption/decryption
   - Only encrypt sensitive fields (PII, credentials, etc.)
   - Use strong encryption algorithms (AES-256-GCM)

2. **Key Management**
   - Secure key storage (environment variables, Vault)
   - Key rotation mechanism
   - Backup and recovery procedures

3. **Implementation Approaches**
   - Database-level encryption (when possible)
   - Application-level encryption (for cross-platform)
   - API-level encryption (for external services)

### 4. Security Headers and Best Practices

Essential security headers:

1. **Content-Security-Policy (CSP)**
   ```
   Content-Security-Policy: default-src 'self'; script-src 'self' 'unsafe-inline'; style-src 'self' 'unsafe-inline';
   ```

2. **Strict-Transport-Security (HSTS)**
   ```
   Strict-Transport-Security: max-age=31536000; includeSubDomains
   ```

3. **X-Content-Type-Options**
   ```
   X-Content-Type-Options: nosniff
   ```

4. **X-Frame-Options**
   ```
   X-Frame-Options: DENY
   ```

5. **Referrer-Policy**
   ```
   Referrer-Policy: strict-origin-when-cross-origin
   ```

6. **Feature-Policy**
   ```
   Feature-Policy: camera 'none'; microphone 'none'; geolocation 'none'
   ```

### 5. Secure Cookie Implementation

Essential cookie security attributes:

1. **Secure Flag**
   - Ensures cookies are only sent over HTTPS

2. **HttpOnly Flag**
   - Prevents JavaScript access to cookies

3. **SameSite Attribute**
   - Options: Strict, Lax, None
   - Helps prevent CSRF attacks

4. **Max-Age/Expires**
   - Set short lifetimes for sensitive cookies

## Implementation Resources

### Security Audit Tools

1. **OWASP ZAP**
   - Web application security scanner
   - Automated scanning and reporting
   - API endpoint testing

2. **SQLMap**
   - Database security testing
   - SQL injection vulnerability detection

3. **Trufflehog**
   - Secret scanning in codebase

### Encryption Libraries

1. **cryptography** (Python)
   - Modern cryptography library
   - Fernet symmetric encryption
   - PBKDF2 key derivation

2. **node-forge** (JavaScript)
   - JavaScript cryptography library
   - AES encryption

### Security Headers Middleware

1. **Flask-Talisman**
   - Flask extension for security headers
   - CSP configuration
   - HTTPS enforcement

2. **Helmet** (for Node.js/Express)
   - Reference for security header implementations

### Cookie Security

1. **Flask-Session**
   - Server-side session management
   - Secure cookie handling

2. **itsdangerous**
   - Secure cookie signing
   - Tamper-proof data

## Implementation Examples

### Security Audit Script

```python
# security/auditor.py
import os
import json
import sys
import requests
from datetime import datetime
import logging

class SecurityAuditor:
    """Comprehensive security audit framework."""
    
    def __init__(self, base_url="http://localhost:5000"):
        self.base_url = base_url
        self.issues = []
        self.timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        self.logger = logging.getLogger("security_audit")
    
    def check_security_headers(self):
        """Check for proper security headers."""
        self.logger.info("Checking security headers")
        headers_to_check = {
            'Content-Security-Policy': None,
            'Strict-Transport-Security': None,
            'X-Content-Type-Options': 'nosniff',
            'X-Frame-Options': ['DENY', 'SAMEORIGIN'],
            'Referrer-Policy': None
        }
        
        try:
            response = requests.get(f"{self.base_url}/")
            
            for header, expected_value in headers_to_check.items():
                if header not in response.headers:
                    self.issues.append({
                        'severity': 'high',
                        'category': 'headers',
                        'title': f"Missing {header} header",
                        'description': f"The {header} security header is missing from HTTP responses.",
                        'remediation': f"Add the {header} header to all HTTP responses."
                    })
                elif expected_value:
                    if isinstance(expected_value, list):
                        if response.headers[header] not in expected_value:
                            self.issues.append({
                                'severity': 'medium',
                                'category': 'headers',
                                'title': f"Invalid {header} value",
                                'description': f"The {header} header has an invalid value: {response.headers[header]}",
                                'remediation': f"Set {header} to one of: {', '.join(expected_value)}"
                            })
                    elif response.headers[header] != expected_value:
                        self.issues.append({
                            'severity': 'medium',
                            'category': 'headers',
                            'title': f"Invalid {header} value",
                            'description': f"The {header} header has an invalid value: {response.headers[header]}",
                            'remediation': f"Set {header} to: {expected_value}"
                        })
        except Exception as e:
            self.logger.error(f"Error checking security headers: {str(e)}")
            self.issues.append({
                'severity': 'error',
                'category': 'general',
                'title': "Error checking security headers",
                'description': f"Exception occurred: {str(e)}",
                'remediation': "Verify the application is running and accessible."
            })
            
    # Additional audit methods would be implemented here...
    
    def run_full_audit(self):
        """Execute complete security audit."""
        self.check_security_headers()
        # Run other audit checks...
        self.generate_report()
    
    def generate_report(self):
        """Generate security audit report."""
        report = {
            'timestamp': self.timestamp,
            'target': self.base_url,
            'issues': self.issues,
            'summary': {
                'critical': len([i for i in self.issues if i['severity'] == 'critical']),
                'high': len([i for i in self.issues if i['severity'] == 'high']),
                'medium': len([i for i in self.issues if i['severity'] == 'medium']),
                'low': len([i for i in self.issues if i['severity'] == 'low']),
                'info': len([i for i in self.issues if i['severity'] == 'info']),
                'total': len(self.issues)
            }
        }
        
        # Ensure output directory exists
        os.makedirs('reports/security', exist_ok=True)
        
        # Write JSON report
        report_file = f"reports/security/audit_{self.timestamp}.json"
        with open(report_file, 'w') as f:
            json.dump(report, f, indent=2)
        
        # Write human-readable report
        report_md = f"reports/security/audit_{self.timestamp}.md"
        with open(report_md, 'w') as f:
            f.write(f"# Security Audit Report\n\n")
            f.write(f"**Date:** {self.timestamp}\n")
            f.write(f"**Target:** {self.base_url}\n\n")
            
            f.write("## Summary\n\n")
            f.write(f"- Critical: {report['summary']['critical']}\n")
            f.write(f"- High: {report['summary']['high']}\n")
            f.write(f"- Medium: {report['summary']['medium']}\n")
            f.write(f"- Low: {report['summary']['low']}\n")
            f.write(f"- Info: {report['summary']['info']}\n")
            f.write(f"- Total: {report['summary']['total']}\n\n")
            
            f.write("## Issues\n\n")
            for i, issue in enumerate(self.issues):
                f.write(f"### {i+1}. {issue['title']} ({issue['severity'].upper()})\n\n")
                f.write(f"**Category:** {issue['category']}\n\n")
                f.write(f"**Description:** {issue['description']}\n\n")
                f.write(f"**Remediation:** {issue['remediation']}\n\n")
            
            f.write("\n---\n")
            f.write("Generated with SecurityAuditor")
        
        self.logger.info(f"Security audit report generated: {report_file}")
        return report_file, report_md
```

### Encryption Service Implementation

```python
# security/encryption.py
from cryptography.fernet import Fernet
from cryptography.hazmat.primitives import hashes
from cryptography.hazmat.primitives.kdf.pbkdf2 import PBKDF2HMAC
import os
import base64
import json
import logging

class EncryptionService:
    """Service for encrypting and decrypting sensitive data."""
    
    _instance = None
    
    @classmethod
    def get_instance(cls, key_file=None, config=None):
        """Get singleton instance of EncryptionService."""
        if cls._instance is None:
            cls._instance = EncryptionService(key_file, config)
        return cls._instance
    
    def __init__(self, key_file=None, config=None):
        """Initialize encryption service."""
        self.logger = logging.getLogger("encryption_service")
        self.config = config or {}
        
        # Key file location
        self.key_file = key_file or os.environ.get(
            'ENCRYPTION_KEY_FILE',
            self.config.get('ENCRYPTION_KEY_FILE', 'keys/encryption.key')
        )
        
        # Load or generate key
        self.key = self._load_or_generate_key()
        self.fernet = Fernet(self.key)
        
        self.logger.info("Encryption service initialized")
    
    def _load_or_generate_key(self):
        """Load existing key or generate new one."""
        try:
            if os.path.exists(self.key_file):
                with open(self.key_file, 'rb') as f:
                    key = f.read()
                    self.logger.info(f"Loaded encryption key from {self.key_file}")
                    return key
            
            # Generate new key if none exists
            key = Fernet.generate_key()
            
            # Ensure directory exists
            os.makedirs(os.path.dirname(self.key_file), exist_ok=True)
            
            # Save key
            with open(self.key_file, 'wb') as f:
                f.write(key)
            
            self.logger.info(f"Generated new encryption key and saved to {self.key_file}")
            return key
            
        except Exception as e:
            self.logger.error(f"Error in key management: {str(e)}")
            # Fall back to environment variable in emergency
            env_key = os.environ.get('ENCRYPTION_KEY')
            if env_key:
                self.logger.warning("Using encryption key from environment variable")
                return env_key.encode()
            
            # Last resort: generate temporary key (will not persist)
            self.logger.warning("Using temporary encryption key - data cannot be decrypted in future sessions")
            return Fernet.generate_key()
    
    def encrypt(self, data):
        """Encrypt data."""
        if data is None:
            return None
            
        try:
            # Convert to JSON and then to bytes if not already bytes
            if not isinstance(data, bytes):
                if isinstance(data, str):
                    data_bytes = data.encode('utf-8')
                else:
                    data_bytes = json.dumps(data).encode('utf-8')
            else:
                data_bytes = data
                
            # Encrypt the data
            encrypted = self.fernet.encrypt(data_bytes)
            
            # Return base64 string representation
            return base64.b64encode(encrypted).decode('utf-8')
            
        except Exception as e:
            self.logger.error(f"Encryption error: {str(e)}")
            raise
    
    def decrypt(self, encrypted_data):
        """Decrypt data."""
        if encrypted_data is None:
            return None
            
        try:
            # Convert from base64 string to bytes if needed
            if isinstance(encrypted_data, str):
                encrypted_bytes = base64.b64decode(encrypted_data.encode('utf-8'))
            else:
                encrypted_bytes = encrypted_data
                
            # Decrypt the data
            decrypted_bytes = self.fernet.decrypt(encrypted_bytes)
            
            # Try to parse as JSON, fall back to string
            try:
                return json.loads(decrypted_bytes.decode('utf-8'))
            except json.JSONDecodeError:
                return decrypted_bytes.decode('utf-8')
                
        except Exception as e:
            self.logger.error(f"Decryption error: {str(e)}")
            raise
    
    def rotate_key(self):
        """Generate new key and return migration helper for re-encryption."""
        old_key = self.key
        old_fernet = self.fernet
        
        # Generate new key
        self.key = Fernet.generate_key()
        self.fernet = Fernet(self.key)
        
        # Save new key
        with open(self.key_file, 'wb') as f:
            f.write(self.key)
        
        # Create helper for data migration
        migration_helper = {
            'old_fernet': old_fernet,
            'new_fernet': self.fernet,
            'reencrypt': lambda data: self.fernet.encrypt(old_fernet.decrypt(data))
        }
        
        self.logger.info("Encryption key rotated successfully")
        return migration_helper
```

### Security Headers Middleware

```python
# api/security_middleware.py
from flask import Flask, request, Response
from functools import wraps
import secrets
import logging

class SecurityHeadersMiddleware:
    """Middleware to add security headers to HTTP responses."""
    
    def __init__(self, app=None, csp_policy=None, hsts_age=31536000):
        self.app = app
        self.csp_policy = csp_policy or self._default_csp_policy()
        self.hsts_age = hsts_age
        self.logger = logging.getLogger("security_headers")
        
        if app is not None:
            self.init_app(app)
    
    def init_app(self, app):
        """Initialize middleware with Flask app."""
        self.app = app
        
        # Add after_request handler
        app.after_request(self.add_security_headers)
        self.logger.info("Security headers middleware initialized")
    
    def _default_csp_policy(self):
        """Default Content Security Policy."""
        return {
            'default-src': ["'self'"],
            'script-src': ["'self'", "'unsafe-inline'"],
            'style-src': ["'self'", "'unsafe-inline'"],
            'img-src': ["'self'", "data:"],
            'font-src': ["'self'"],
            'connect-src': ["'self'"],
            'frame-ancestors': ["'none'"],
            'form-action': ["'self'"]
        }
    
    def _build_csp_header(self):
        """Build Content-Security-Policy header value."""
        csp_parts = []
        for directive, sources in self.csp_policy.items():
            csp_parts.append(f"{directive} {' '.join(sources)}")
        return "; ".join(csp_parts)
    
    def add_security_headers(self, response):
        """Add security headers to response."""
        try:
            # Content Security Policy
            response.headers['Content-Security-Policy'] = self._build_csp_header()
            
            # HTTP Strict Transport Security
            if not request.is_secure and not self.app.debug:
                response.headers['Strict-Transport-Security'] = f"max-age={self.hsts_age}; includeSubDomains"
            
            # Prevent MIME type sniffing
            response.headers['X-Content-Type-Options'] = 'nosniff'
            
            # Control frame embedding
            response.headers['X-Frame-Options'] = 'DENY'
            
            # Enable browser XSS protection
            response.headers['X-XSS-Protection'] = '1; mode=block'
            
            # Control referrer information
            response.headers['Referrer-Policy'] = 'strict-origin-when-cross-origin'
            
            # Restrict browser features
            response.headers['Feature-Policy'] = "geolocation 'none'; microphone 'none'; camera 'none'"
            
            return response
            
        except Exception as e:
            self.logger.error(f"Error adding security headers: {str(e)}")
            return response

# Usage:
# app = Flask(__name__)
# security_headers = SecurityHeadersMiddleware(app)
```

These implementation examples provide a starting point for each security component. Adapt and integrate them into the CryoProtect v2 codebase as needed.