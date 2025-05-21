# Security Module Reference

## Overview

The security module implements comprehensive security features for CryoProtect, including authentication, authorization, data protection, and security best practices throughout the application.

## Key Components

### Authentication System
- **JWT Authentication**: `jwt_auth.py` for token-based authentication
- **Enhanced JWT**: `enhanced_jwt_auth.py` with additional security features
- **Service Role Auth**: Service account authentication for system operations
- **Session Management**: `session_management.py` for secure session handling

### Authorization Framework
- **Role-Based Access Control**: `rbac.py` implementation
- **Permission System**: Granular permission definitions
- **Policy Enforcement**: Authorization policy validation
- **Row-Level Security**: Integration with database RLS policies

### Data Protection
- **CSRF Protection**: Cross-Site Request Forgery prevention in `csrf.py`
- **Input Validation**: Request validation with schemas
- **Output Sanitization**: Response sanitization
- **Secure Headers**: Security headers implementation in `security_headers.py`

## Authentication Flows

### User Authentication Flow
1. **Login Request**: User credentials submitted to auth endpoint
2. **Credential Verification**: Validation against Supabase Auth
3. **JWT Generation**: Creation of signed JWT with user claims
4. **Token Response**: Return of access and refresh tokens
5. **Subsequent Requests**: JWT included in Authorization header
6. **Token Verification**: Validation of signature and claims
7. **Token Refresh**: Obtaining new tokens with refresh token

### Service Role Authentication
1. **Service Credentials**: Service API key and secret
2. **Service JWT**: Special JWT with service role claims
3. **RLS Bypass**: Ability to bypass row-level security
4. **Permission Elevation**: Access to administrative operations

## Row-Level Security

The system implements comprehensive Row-Level Security (RLS) in the database:

- **User-Based Policies**: Records accessible only to owning users
- **Organization Policies**: Data segregation by organization
- **Role-Based Policies**: Access based on user roles
- **Public Data Policies**: Controlled public data access
- **Security Definer Functions**: Optimized policy enforcement

## Security Headers

Security headers implemented in `security_headers.py`:

- **Content-Security-Policy**: Restrict resource loading
- **X-Content-Type-Options**: Prevent MIME type sniffing
- **X-Frame-Options**: Control iframe embedding
- **Strict-Transport-Security**: Enforce HTTPS
- **X-XSS-Protection**: XSS filtering (legacy browsers)
- **Referrer-Policy**: Control referrer information
- **Permissions-Policy**: Restrict browser features

## Best Practices

1. **Defense in Depth**: Multiple security layers
2. **Least Privilege**: Minimal permissions needed
3. **Secure by Default**: Security enabled without configuration
4. **Input Validation**: Validate all user input
5. **Secure Communications**: Encrypted data in transit
6. **Audit Logging**: Log security-relevant events
7. **Regular Testing**: Security testing and review

## Common Vulnerabilities Addressed

1. **SQL Injection**: Parameterized queries and ORM
2. **XSS (Cross-Site Scripting)**: Content sanitization and CSP
3. **CSRF (Cross-Site Request Forgery)**: Anti-CSRF tokens
4. **SSRF (Server-Side Request Forgery)**: URL validation
5. **Broken Authentication**: Secure session management
6. **Sensitive Data Exposure**: Proper data encryption
7. **Insecure Deserialization**: Safe deserialization practices

## Security Configuration

Security settings are managed in `config.py` with these key parameters:

- **JWT_SECRET_KEY**: Secret for JWT signing
- **JWT_ALGORITHM**: Algorithm used for JWT (RS256/HS256)
- **TOKEN_EXPIRY**: Token lifetime configuration
- **CSRF_SECRET**: Secret for CSRF token generation
- **PASSWORD_POLICY**: Password strength requirements
- **RATE_LIMIT_CONFIG**: API rate limiting settings
- **CSP_POLICY**: Content Security Policy configuration

## Usage Examples

```python
# Protecting an endpoint with authentication
@jwt_required
def protected_resource():
    # Only authenticated users reach here
    pass

# Role-based authorization
@requires_role('admin')
def admin_resource():
    # Only admin users reach here
    pass

# Permission-based authorization
@requires_permission('write:data')
def write_data_resource():
    # Only users with write permission reach here
    pass

# CSRF protection
@csrf_protect
def form_submission():
    # Protected against CSRF attacks
    pass

# Add security headers to response
@app.after_request
def add_security_headers(response):
    return set_security_headers(response)
```