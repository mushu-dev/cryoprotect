# Phase 3.3: Security Implementation Line References

This document provides specific file and line references for implementing the Security phase of CryoProtect v2. Use these references to quickly locate the relevant code sections without extensive searching.

## 1. Comprehensive Security Audit

### Key Files to Examine
- **app.py**: Main Flask application with route registrations and middleware
- **api/__init__.py**: API route registrations
- **api/resources.py**: Base resource class definitions
- **api/jwt_auth.py**: JWT authentication implementation
- **api/enhanced_jwt_auth.py**: Enhanced JWT authentication
- **api/session_management.py**: User session handling
- **api/rbac.py**: Role-based access control implementation
- **migrations/006_rls_policies.sql**: RLS policy definitions
- **migrations/007_service_role_rls.sql**: Service role RLS policies
- **config.py**: Application configuration with security settings

### Implementation Locations
- Create new file: **security/auditor.py** - Security audit script
- Create new file: **security/scanner.py** - Vulnerability scanning tools
- Create new file: **reports/security/audit_template.md** - Audit report template

## 2. Vulnerability Scanning Integration

### Key Files to Examine
- **.github/workflows/deploy.yml**: Current CI/CD workflow
- **.github/workflows/ci-cd.yml**: Additional CI/CD workflow
- **requirements.txt**: Python dependencies
- **package.json**: Node.js dependencies

### Implementation Locations
- Update file: **.github/workflows/deploy.yml** - Add security scanning steps
- Update file: **.github/workflows/ci-cd.yml** - Add vulnerability scanning
- Create new file: **security/scan_dependencies.py** - Dependency scanning script
- Create new file: **security/vulnerability_manager.py** - Vulnerability management

## 3. Enhanced Data Encryption

### Key Files to Examine
- **config.py**: Application configuration
- **api/models.py**: Database models
- **api/resources.py**: API resources that handle sensitive data
- **api/user_profile_resources.py**: User profile data handling

### Implementation Locations
- Create new file: **security/encryption.py** - Encryption service implementation
- Update file: **api/models.py** - Add encryption to sensitive fields
- Update file: **config.py** - Add encryption configuration
- Create new file: **scripts/rotate_encryption_keys.py** - Key rotation script
- Create new file: **tests/test_encryption.py** - Tests for encryption service

## 4. Security Best Practices Implementation

### Key Files to Examine
- **app.py**: Main Flask application
- **api/__init__.py**: API initialization
- **api/utils.py**: Utility functions
- **templates/base.html**: Base HTML template
- **static/js/app.js**: Frontend JavaScript

### Implementation Locations
- Create new file: **api/security_middleware.py** - Security headers middleware
- Update file: **app.py** - Apply security middleware
- Update file: **api/session_management.py** - Enhance cookie security
- Update file: **templates/base.html** - Add CSP nonce integration
- Update file: **static/js/api.js** - Add CSRF protection

## 5. Documentation and Checklist

### Implementation Locations
- Create new file: **docs/security_implementation.md** - Security implementation documentation
- Create new file: **docs/security_checklist.md** - Security checklist for new features
- Create new file: **docs/security_response.md** - Security incident response procedures
- Create new file: **docs/security_training.md** - Security awareness training

## Code References for Key Security Components

### JWT Authentication
```python
# api/jwt_auth.py - JWT Authentication Integration

def generate_token(user_id, role='user'):
    """Generate a JWT token for a user with specific role."""
    payload = {
        'exp': datetime.utcnow() + timedelta(days=1),
        'iat': datetime.utcnow(),
        'sub': str(user_id),
        'role': role
    }
    return jwt.encode(
        payload,
        current_app.config.get('SECRET_KEY'),
        algorithm='HS256'
    )
```

### Role-Based Access Control
```python
# api/rbac.py - Role-Based Access Control

def role_required(roles):
    """Decorator to restrict access based on user role."""
    def decorator(f):
        @wraps(f)
        def decorated_function(*args, **kwargs):
            token = extract_token_from_request()
            if not token:
                return jsonify({'message': 'Authentication required'}), 401
            
            try:
                payload = verify_token(token)
                role = payload.get('role', 'user')
                if role not in roles:
                    return jsonify({'message': 'Insufficient permissions'}), 403
            except:
                return jsonify({'message': 'Invalid token'}), 401
            
            return f(*args, **kwargs)
        return decorated_function
    return decorator
```

### Row Level Security
```sql
-- migrations/006_rls_policies.sql - Example RLS Policy

-- Enable RLS on table
ALTER TABLE molecules ENABLE ROW LEVEL SECURITY;

-- Create policy for users
CREATE POLICY user_molecules_policy ON molecules 
    USING (is_public = true OR created_by = auth.uid());
    
-- Create policy for team access
CREATE POLICY team_molecules_policy ON molecules 
    USING (EXISTS (
        SELECT 1 FROM user_teams ut 
        WHERE ut.user_id = auth.uid() 
        AND ut.team_id = molecules.team_id
    ));
```

### API Rate Limiting
```python
# api/rate_limiter.py - Rate Limiting Implementation

class RateLimiter:
    """Rate limiting implementation for API endpoints."""
    
    def __init__(self, redis_client, limit=100, window=3600):
        self.redis = redis_client
        self.limit = limit
        self.window = window
    
    def is_rate_limited(self, key):
        """Check if the given key is rate limited."""
        current = self.redis.get(key)
        if current is None:
            self.redis.set(key, 1, ex=self.window)
            return False
        
        current = int(current)
        if current >= self.limit:
            return True
            
        self.redis.incr(key)
        return False
```