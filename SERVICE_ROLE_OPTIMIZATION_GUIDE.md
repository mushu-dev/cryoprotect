# Service Role Authentication Optimization Guide

This guide provides detailed implementation steps for enhancing the service role authentication mechanism in the CryoProtect system. It serves as a companion document to the main CRYOPROTECT_OPTIMIZATION_PLAN.md and focuses specifically on service role security improvements.

## 1. Current Service Role Implementation Analysis

### Current Limitations
1. **Direct service role key usage**: The service role key is directly used in application code
2. **Limited token management**: No proper token lifecycle management
3. **Inconsistent RLS policies**: Service role bypass policies are not consistently applied
4. **No proper role claims**: JWT tokens lack proper role claims
5. **No token revocation**: No mechanism to revoke compromised tokens
6. **Limited audit trail**: Insufficient logging of service role operations

### Security and Performance Impact
- Increased risk of service role key exposure
- Difficulty tracking service role usage
- Inconsistent permissions across tables
- Complex application code to handle service role authentication
- Limited ability to restrict service role capabilities
- Potential for unintended data access

## 2. Optimization Steps

### 2.1 JWT-Based Service Role Authentication

#### Problem
The current implementation directly uses the service role key in application code, which is both a security risk and makes it difficult to control and track service role access.

#### Solution
Implement a JWT-based service role authentication mechanism that:
- Generates short-lived JWT tokens with specific permissions
- Includes proper role claims and scopes
- Has configurable expiration
- Provides token revocation capability
- Creates an audit trail of service role token usage

#### Implementation Details

**1. Service Role Token Generator Class**

```python
# service_role_token.py
import time
import uuid
import jwt
import logging
from typing import Dict, List, Optional, Any

logger = logging.getLogger(__name__)

class ServiceRoleTokenManager:
    """
    Manages JWT tokens for service role authentication.
    
    This class is responsible for generating, validating, and revoking
    service role tokens with specific permissions.
    """
    
    def __init__(self, service_role_key: str, jwt_secret: str, token_expiry: int = 3600):
        """
        Initialize the token manager.
        
        Args:
            service_role_key: The Supabase service role key
            jwt_secret: Secret used to sign the JWT tokens
            token_expiry: Token expiration time in seconds (default: 1 hour)
        """
        self.service_role_key = service_role_key
        self.jwt_secret = jwt_secret
        self.token_expiry = token_expiry
        self.revoked_tokens = set()
        
    def generate_token(self, 
                       scopes: List[str], 
                       user_id: Optional[str] = None, 
                       metadata: Optional[Dict[str, Any]] = None) -> str:
        """
        Generate a service role JWT token with specific scopes.
        
        Args:
            scopes: List of permission scopes (e.g., 'molecules:read', 'molecular_properties:write')
            user_id: Optional user ID to associate with the token
            metadata: Optional metadata to include in the token
            
        Returns:
            JWT token string
        """
        now = int(time.time())
        token_id = str(uuid.uuid4())
        
        payload = {
            'iss': 'cryoprotect',                  # Issuer
            'sub': user_id or 'service_role',      # Subject
            'role': 'service_role',                # Role claim
            'scopes': scopes,                      # Permission scopes
            'jti': token_id,                       # JWT ID for revocation
            'iat': now,                            # Issued at time
            'exp': now + self.token_expiry,        # Expiration time
            'metadata': metadata or {}             # Additional metadata
        }
        
        token = jwt.encode(payload, self.jwt_secret, algorithm='HS256')
        
        # Log token generation
        logger.info(f"Generated service role token {token_id[:8]} for {'user ' + user_id if user_id else 'service'} with scopes: {', '.join(scopes)}")
        
        return token
    
    def validate_token(self, token: str) -> Dict[str, Any]:
        """
        Validate a service role JWT token.
        
        Args:
            token: JWT token string
            
        Returns:
            Token payload if valid
            
        Raises:
            jwt.InvalidTokenError: If token is invalid
            ValueError: If token is revoked or expired
        """
        try:
            payload = jwt.decode(token, self.jwt_secret, algorithms=['HS256'])
            
            # Check if token is revoked
            if payload['jti'] in self.revoked_tokens:
                raise ValueError("Token has been revoked")
            
            return payload
        except jwt.ExpiredSignatureError:
            logger.warning(f"Attempted to use expired token")
            raise ValueError("Token has expired")
        except jwt.InvalidTokenError as e:
            logger.warning(f"Invalid token: {str(e)}")
            raise
    
    def revoke_token(self, token: str) -> bool:
        """
        Revoke a JWT token.
        
        Args:
            token: JWT token string
            
        Returns:
            True if token was revoked, False otherwise
        """
        try:
            payload = jwt.decode(token, self.jwt_secret, algorithms=['HS256'], options={'verify_exp': False})
            self.revoked_tokens.add(payload['jti'])
            
            logger.info(f"Revoked service role token {payload['jti'][:8]}")
            return True
        except Exception as e:
            logger.error(f"Failed to revoke token: {str(e)}")
            return False
    
    def has_scope(self, token: str, required_scope: str) -> bool:
        """
        Check if a token has a specific scope.
        
        Args:
            token: JWT token string
            required_scope: The scope to check for
            
        Returns:
            True if token has the required scope, False otherwise
        """
        try:
            payload = self.validate_token(token)
            scopes = payload.get('scopes', [])
            
            # Check for exact match
            if required_scope in scopes:
                return True
                
            # Check for wildcard matches (e.g., 'molecules:*' matches 'molecules:read')
            for scope in scopes:
                if scope.endswith(':*') and required_scope.startswith(scope[:-1]):
                    return True
            
            return False
        except:
            return False
```

**2. Enhanced JWT Authentication Middleware**

```python
# enhanced_jwt_auth.py
from functools import wraps
from flask import request, jsonify, g, current_app
import jwt
from typing import List, Optional, Callable

def jwt_required(scopes: Optional[List[str]] = None):
    """
    Decorator to protect routes with JWT authentication.
    
    Args:
        scopes: Optional list of required permission scopes
        
    Returns:
        Decorated function
    """
    def decorator(f: Callable):
        @wraps(f)
        def decorated_function(*args, **kwargs):
            auth_header = request.headers.get('Authorization', '')
            
            if not auth_header.startswith('Bearer '):
                return jsonify({"error": "Authorization header with Bearer token required"}), 401
            
            token = auth_header.split(' ')[1]
            
            try:
                # Check if it's a service role token
                if token.startswith('srt_'):
                    # Validate service role token
                    token_manager = current_app.service_role_token_manager
                    payload = token_manager.validate_token(token)
                    
                    # Check scopes if provided
                    if scopes:
                        for scope in scopes:
                            if not token_manager.has_scope(token, scope):
                                return jsonify({"error": f"Insufficient permissions. Required scope: {scope}"}), 403
                    
                    # Set token payload in g for access in routes
                    g.jwt_payload = payload
                    g.jwt_type = 'service_role'
                    
                else:
                    # Regular user token validation
                    payload = jwt.decode(token, current_app.config['JWT_SECRET_KEY'], algorithms=['HS256'])
                    
                    # Set token payload in g for access in routes
                    g.jwt_payload = payload
                    g.jwt_type = 'user'
                
                return f(*args, **kwargs)
                
            except jwt.ExpiredSignatureError:
                return jsonify({"error": "Token has expired"}), 401
            except jwt.InvalidTokenError as e:
                return jsonify({"error": f"Invalid token: {str(e)}"}), 401
            except ValueError as e:
                return jsonify({"error": str(e)}), 401
                
        return decorated_function
    return decorator

def service_role_required(scopes: Optional[List[str]] = None):
    """
    Decorator to protect routes with service role authentication.
    
    Args:
        scopes: Optional list of required permission scopes
        
    Returns:
        Decorated function
    """
    def decorator(f: Callable):
        @wraps(f)
        def decorated_function(*args, **kwargs):
            auth_header = request.headers.get('Authorization', '')
            
            if not auth_header.startswith('Bearer '):
                return jsonify({"error": "Authorization header with Bearer token required"}), 401
            
            token = auth_header.split(' ')[1]
            
            # Only allow service role tokens
            if not token.startswith('srt_'):
                return jsonify({"error": "Service role token required"}), 403
            
            try:
                # Validate service role token
                token_manager = current_app.service_role_token_manager
                payload = token_manager.validate_token(token)
                
                # Check scopes if provided
                if scopes:
                    for scope in scopes:
                        if not token_manager.has_scope(token, scope):
                            return jsonify({"error": f"Insufficient permissions. Required scope: {scope}"}), 403
                
                # Set token payload in g for access in routes
                g.jwt_payload = payload
                g.jwt_type = 'service_role'
                
                return f(*args, **kwargs)
                
            except jwt.ExpiredSignatureError:
                return jsonify({"error": "Token has expired"}), 401
            except jwt.InvalidTokenError as e:
                return jsonify({"error": f"Invalid token: {str(e)}"}), 401
            except ValueError as e:
                return jsonify({"error": str(e)}), 401
                
        return decorated_function
    return decorator
```

**3. Service Role Helper Class**

```python
# service_role_helper.py
import requests
import json
import logging
from typing import Dict, List, Any, Optional

logger = logging.getLogger(__name__)

class ServiceRoleClient:
    """
    Client for making requests using service role authentication.
    
    This class encapsulates the logic for making authenticated
    requests to Supabase using the service role.
    """
    
    def __init__(self, base_url: str, token_manager):
        """
        Initialize the service role client.
        
        Args:
            base_url: The Supabase API base URL
            token_manager: The service role token manager instance
        """
        self.base_url = base_url
        self.token_manager = token_manager
        
    def _get_token(self, scopes: List[str], user_id: Optional[str] = None) -> str:
        """Get service role token with required scopes."""
        return self.token_manager.generate_token(scopes, user_id)
        
    def execute_query(self, query: str, params: Optional[Dict[str, Any]] = None, 
                     user_id: Optional[str] = None) -> Dict[str, Any]:
        """
        Execute a SQL query with service role permissions.
        
        Args:
            query: SQL query to execute
            params: Query parameters
            user_id: Optional user ID to associate with the operation
            
        Returns:
            Query result
        """
        token = self._get_token(['sql:execute'], user_id)
        
        headers = {
            'Authorization': f'Bearer {token}',
            'Content-Type': 'application/json'
        }
        
        payload = {
            'query': query,
            'params': params or []
        }
        
        try:
            response = requests.post(
                f"{self.base_url}/rest/v1/rpc/execute_sql",
                headers=headers,
                json=payload
            )
            response.raise_for_status()
            return response.json()
        except Exception as e:
            logger.error(f"Error executing query with service role: {str(e)}")
            raise
            
    def insert_data(self, table: str, data: Dict[str, Any], 
                   user_id: Optional[str] = None) -> Dict[str, Any]:
        """
        Insert data into a table with service role permissions.
        
        Args:
            table: Table name
            data: Data to insert
            user_id: Optional user ID to associate with the operation
            
        Returns:
            Inserted record
        """
        token = self._get_token([f'{table}:insert'], user_id)
        
        headers = {
            'Authorization': f'Bearer {token}',
            'Content-Type': 'application/json',
            'Prefer': 'return=representation'
        }
        
        try:
            response = requests.post(
                f"{self.base_url}/rest/v1/{table}",
                headers=headers,
                json=data
            )
            response.raise_for_status()
            return response.json()
        except Exception as e:
            logger.error(f"Error inserting data with service role: {str(e)}")
            raise
    
    def update_data(self, table: str, id_column: str, id_value: Any, 
                   data: Dict[str, Any], user_id: Optional[str] = None) -> Dict[str, Any]:
        """
        Update data in a table with service role permissions.
        
        Args:
            table: Table name
            id_column: Primary key column name
            id_value: Primary key value
            data: Data to update
            user_id: Optional user ID to associate with the operation
            
        Returns:
            Updated record
        """
        token = self._get_token([f'{table}:update'], user_id)
        
        headers = {
            'Authorization': f'Bearer {token}',
            'Content-Type': 'application/json',
            'Prefer': 'return=representation'
        }
        
        try:
            response = requests.patch(
                f"{self.base_url}/rest/v1/{table}?{id_column}=eq.{id_value}",
                headers=headers,
                json=data
            )
            response.raise_for_status()
            return response.json()
        except Exception as e:
            logger.error(f"Error updating data with service role: {str(e)}")
            raise
    
    def delete_data(self, table: str, id_column: str, id_value: Any, 
                   user_id: Optional[str] = None) -> Dict[str, Any]:
        """
        Delete data from a table with service role permissions.
        
        Args:
            table: Table name
            id_column: Primary key column name
            id_value: Primary key value
            user_id: Optional user ID to associate with the operation
            
        Returns:
            Deleted record
        """
        token = self._get_token([f'{table}:delete'], user_id)
        
        headers = {
            'Authorization': f'Bearer {token}',
            'Content-Type': 'application/json',
            'Prefer': 'return=representation'
        }
        
        try:
            response = requests.delete(
                f"{self.base_url}/rest/v1/{table}?{id_column}=eq.{id_value}",
                headers=headers
            )
            response.raise_for_status()
            return response.json()
        except Exception as e:
            logger.error(f"Error deleting data with service role: {str(e)}")
            raise
```

**4. Application Integration**

```python
# app.py
from flask import Flask
from config import get_config
from service_role_token import ServiceRoleTokenManager
from service_role_helper import ServiceRoleClient

def create_app():
    app = Flask(__name__)
    
    # Load configuration
    config = get_config()
    app.config.from_object(config)
    
    # Initialize service role token manager
    app.service_role_token_manager = ServiceRoleTokenManager(
        service_role_key=app.config['SUPABASE_SERVICE_KEY'],
        jwt_secret=app.config['JWT_SECRET_KEY'],
        token_expiry=app.config.get('SERVICE_ROLE_TOKEN_EXPIRY', 3600)
    )
    
    # Initialize service role client
    app.service_role_client = ServiceRoleClient(
        base_url=app.config['SUPABASE_URL'],
        token_manager=app.service_role_token_manager
    )
    
    # Register blueprints and routes...
    
    return app
```

### 2.2 Service Role RLS Policy Unification

#### Problem
The current implementation has inconsistent RLS policies for service role access across tables. This makes it difficult to maintain and can lead to unintended access restrictions.

#### Solution
Create a unified approach to service role RLS policies:
- Add a consistent policy for service role access to all tables
- Ensure policies are properly documented
- Implement audit logging for service role operations

#### Implementation Details

**Service Role Policy Script**

```sql
-- service_role_policies.sql

-- Function to execute this script as a migration
CREATE OR REPLACE FUNCTION apply_service_role_policies()
RETURNS void AS $$
DECLARE
    table_name text;
    policy_exists boolean;
BEGIN
    -- Create an audit log table if it doesn't exist
    CREATE TABLE IF NOT EXISTS service_role_access_log (
        id SERIAL PRIMARY KEY,
        table_name TEXT NOT NULL,
        operation TEXT NOT NULL,
        record_id UUID,
        user_id UUID,
        access_time TIMESTAMP WITH TIME ZONE DEFAULT now(),
        client_info JSONB
    );
    
    -- Enable row level security on the audit log table
    ALTER TABLE service_role_access_log ENABLE ROW LEVEL SECURITY;
    
    -- Create policy for service role to write to audit log
    DROP POLICY IF EXISTS "service_role_can_insert_logs" ON service_role_access_log;
    CREATE POLICY "service_role_can_insert_logs" ON service_role_access_log
    FOR INSERT TO service_role
    WITH CHECK (true);
    
    -- Create policy for administrators to read audit log
    DROP POLICY IF EXISTS "admins_can_read_logs" ON service_role_access_log;
    CREATE POLICY "admins_can_read_logs" ON service_role_access_log
    FOR SELECT TO authenticated
    USING (
        EXISTS (
            SELECT 1 FROM user_profile
            WHERE auth_user_id = auth.uid()
            AND is_admin = true
        )
    );

    -- Create audit function
    CREATE OR REPLACE FUNCTION log_service_role_access()
    RETURNS TRIGGER AS $$
    BEGIN
        IF current_setting('role') = 'service_role' THEN
            INSERT INTO service_role_access_log (
                table_name, 
                operation, 
                record_id,
                user_id,
                client_info
            )
            VALUES (
                TG_TABLE_NAME,
                TG_OP,
                CASE 
                    WHEN TG_OP = 'DELETE' THEN OLD.id
                    ELSE NEW.id
                END,
                CASE
                    WHEN TG_OP = 'DELETE' THEN OLD.created_by
                    ELSE NEW.created_by
                END,
                jsonb_build_object(
                    'client_ip', current_setting('request.headers', true)::jsonb->>'x-forwarded-for',
                    'user_agent', current_setting('request.headers', true)::jsonb->>'user-agent'
                )
            );
        END IF;
        
        RETURN NULL;
    END;
    $$ LANGUAGE plpgsql SECURITY DEFINER;

    -- Iterate through all public tables
    FOR table_name IN 
        SELECT tablename FROM pg_tables 
        WHERE schemaname = 'public' 
        AND tablename NOT LIKE 'pg_%'
        AND tablename NOT IN ('schema_migrations', 'schema_log', 'service_role_access_log')
    LOOP
        -- Check if service role policy exists
        SELECT EXISTS (
            SELECT 1 FROM pg_policies 
            WHERE schemaname = 'public' 
            AND tablename = table_name 
            AND policyname = 'service_role_all_access'
        ) INTO policy_exists;

        -- Create unified service role policy if it doesn't exist
        IF NOT policy_exists THEN
            EXECUTE format('
                DROP POLICY IF EXISTS "service_role_all_access" ON %I;
                CREATE POLICY "service_role_all_access" ON %I
                FOR ALL
                TO service_role
                USING (true)
                WITH CHECK (true);
            ', table_name, table_name);
            
            RAISE NOTICE 'Created service role policy for table: %', table_name;
        END IF;
        
        -- Create audit log trigger if it has an id and created_by column
        IF EXISTS (
            SELECT 1 FROM information_schema.columns 
            WHERE table_schema = 'public' 
            AND table_name = table_name 
            AND column_name = 'id'
        ) AND EXISTS (
            SELECT 1 FROM information_schema.columns 
            WHERE table_schema = 'public' 
            AND table_name = table_name 
            AND column_name = 'created_by'
        ) THEN
            EXECUTE format('
                DROP TRIGGER IF EXISTS service_role_audit_trigger ON %I;
                CREATE TRIGGER service_role_audit_trigger
                AFTER INSERT OR UPDATE OR DELETE ON %I
                FOR EACH ROW
                EXECUTE FUNCTION log_service_role_access();
            ', table_name, table_name);
            
            RAISE NOTICE 'Created audit trigger for table: %', table_name;
        END IF;
    END LOOP;
END;
$$ LANGUAGE plpgsql;

-- Execute the function
SELECT apply_service_role_policies();
```

### 2.3 Application Updates

#### Problem
The current application code directly uses the service role key, which is both insecure and inflexible.

#### Solution
Update application code to use the new JWT-based service role authentication:
- Replace direct service role key usage with service role client
- Add proper error handling for service role operations
- Implement token-based scopes for different operations

#### Implementation Examples

**1. API Resources Update**

```python
# resources.py
from flask import g, current_app, request, jsonify
from functools import wraps
from .enhanced_jwt_auth import jwt_required, service_role_required

# Service role operation decorator
def service_role_operation(scopes):
    """
    Decorator for operations that require service role access.
    
    This decorator converts a regular user operation to a service role
    operation by using the service role client with the user's ID.
    """
    def decorator(f):
        @wraps(f)
        @jwt_required()
        def decorated_function(*args, **kwargs):
            # Get the current user ID
            user_id = g.jwt_payload.get('sub')
            
            # Store original JWT type for later
            original_jwt_type = g.jwt_type
            
            try:
                # Set service role client in g for access in the route
                g.service_role_client = current_app.service_role_client
                g.service_role_user_id = user_id
                g.service_role_scopes = scopes
                
                # Mark as service role operation
                g.jwt_type = 'service_role_operation'
                
                return f(*args, **kwargs)
            finally:
                # Restore original JWT type
                g.jwt_type = original_jwt_type
                
        return decorated_function
    return decorator

# Example resource with service role operation
class MoleculeResource:
    @jwt_required()
    def get(self, molecule_id):
        # Regular user access
        return jsonify({"message": "Molecule details"})
    
    @service_role_operation(['molecules:insert'])
    def post(self):
        # This operation will use service role with user context
        client = g.service_role_client
        user_id = g.service_role_user_id
        
        # Get request data
        data = request.get_json()
        
        # Add user ID as created_by
        data['created_by'] = user_id
        
        # Insert data using service role
        result = client.insert_data('molecules', data, user_id)
        
        return jsonify(result)
        
    @service_role_required(['molecules:admin'])
    def delete(self, molecule_id):
        # This endpoint is only accessible with service role token
        # that has the molecules:admin scope
        client = g.service_role_client
        
        # Delete data using service role
        result = client.delete_data('molecules', 'id', molecule_id)
        
        return jsonify(result)
```

**2. Background Tasks Update**

```python
# background_tasks.py
import threading
import logging
from flask import current_app

logger = logging.getLogger(__name__)

def process_molecular_calculations(molecule_ids):
    """
    Background task to process molecular calculations.
    
    Args:
        molecule_ids: List of molecule IDs to process
    """
    # Get service role client from application
    service_role_client = current_app.service_role_client
    
    # Generate service role token with appropriate scopes
    scopes = [
        'molecules:read', 
        'molecular_properties:insert',
        'calculation_methods:read'
    ]
    
    try:
        for molecule_id in molecule_ids:
            # Get molecule data
            query = "SELECT * FROM molecules WHERE id = $1"
            result = service_role_client.execute_query(query, [molecule_id])
            
            if not result or not result.get('data'):
                logger.warning(f"Molecule not found: {molecule_id}")
                continue
                
            molecule = result['data'][0]
            
            # Calculate properties
            # ...
            
            # Store calculated properties
            property_data = {
                'molecule_id': molecule_id,
                'property_name': 'calculated_property',
                'property_value': '123.45',
                'property_unit': 'units',
                'calculation_method': 'background_calculation',
                'created_by': molecule['created_by']  # Preserve original creator
            }
            
            # Insert using service role with original user context
            service_role_client.insert_data(
                'molecular_properties', 
                property_data, 
                user_id=molecule['created_by']
            )
            
            logger.info(f"Processed molecule {molecule_id}")
            
    except Exception as e:
        logger.error(f"Error in background processing: {str(e)}")
```

## 3. Implementation Guide

### 3.1 Step-by-Step Approach

1. **Create JWT-Based Token Mechanism**
   - Implement the ServiceRoleTokenManager class
   - Add configuration parameters for token management
   - Implement token revocation mechanism

2. **Implement Enhanced JWT Authentication**
   - Create enhanced JWT auth middleware
   - Add scope-based permission checks
   - Update Flask application to use new authentication

3. **Create Service Role Client**
   - Implement the ServiceRoleClient class
   - Update application code to use client
   - Add proper error handling and logging

4. **Unify RLS Policies**
   - Implement service_role_policies.sql script
   - Apply to all database tables
   - Set up audit logging for service role operations

5. **Update Application Resources**
   - Update API resources to use new service role client
   - Add service_role_operation decorator
   - Implement proper error handling

### 3.2 Configuration Changes

Update the configuration to include the new service role parameters:

```python
# In config.py
class BaseConfig:
    # Existing configuration...

    # Service role configuration
    SERVICE_ROLE_TOKEN_EXPIRY: int = 3600  # 1 hour
    SERVICE_ROLE_JWT_SECRET: str = os.getenv('SERVICE_ROLE_JWT_SECRET', 'change-me-in-production')
    SERVICE_ROLE_ENABLE_AUDIT: bool = True
    SERVICE_ROLE_DEFAULT_SCOPES: List[str] = [
        'molecules:read',
        'molecular_properties:read',
        'mixtures:read',
        'experiments:read'
    ]
```

### 3.3 Testing Recommendations

1. **Authentication Testing**
   - Test token generation and validation
   - Test scope-based access control
   - Test token revocation
   - Test expired token handling

2. **Security Testing**
   - Verify audit logging works correctly
   - Test for token tampering vulnerabilities
   - Test RLS policy effectiveness
   - Test unauthorized access attempts

3. **API Testing**
   - Test all API endpoints with service role
   - Test service_role_operation decorator
   - Test background tasks with service role
   - Test error handling

## 4. Security Considerations

- Store JWT secret securely and never expose it
- Use short token expiration times (1 hour or less)
- Implement proper token revocation for security incidents
- Use specific scopes rather than wildcard permissions
- Always audit service role operations
- Encrypt sensitive data in tokens
- Regularly rotate JWT secrets

## 5. Performance Expectations

### Before Optimization
- Service role authentication overhead: Minimal
- Audit trail: None or manual
- Token management: None
- Security granularity: All-or-nothing

### After Optimization
- Service role authentication overhead: Small (JWT validation)
- Audit trail: Comprehensive
- Token management: Automated with revocation capability
- Security granularity: Fine-grained with scopes

## 6. Monitoring and Audit

### Key Audit Events
- Service role token generation
- Service role token validation failures
- Service role token revocation
- Service role operations by table and operation type
- Service role permission denials

### Audit Log Query Examples

```sql
-- Recent service role operations
SELECT 
    table_name,
    operation,
    COUNT(*) as operation_count
FROM service_role_access_log
WHERE access_time > NOW() - INTERVAL '1 day'
GROUP BY table_name, operation
ORDER BY operation_count DESC;

-- Service role operations by user
SELECT 
    u.email,
    sal.table_name,
    sal.operation,
    COUNT(*) as operation_count
FROM service_role_access_log sal
JOIN auth.users u ON sal.user_id = u.id
WHERE access_time > NOW() - INTERVAL '1 week'
GROUP BY u.email, sal.table_name, sal.operation
ORDER BY u.email, operation_count DESC;

-- Potential suspicious activity
SELECT 
    table_name,
    operation,
    record_id,
    user_id,
    access_time,
    client_info
FROM service_role_access_log
WHERE (
    -- Multiple operations in short time
    user_id IN (
        SELECT user_id
        FROM service_role_access_log
        WHERE access_time > NOW() - INTERVAL '5 minutes'
        GROUP BY user_id
        HAVING COUNT(*) > 100
    )
    -- Or unusual operations
    OR operation = 'DELETE'
    -- Or unusual times
    OR EXTRACT(HOUR FROM access_time) BETWEEN 0 AND 5
)
ORDER BY access_time DESC;
```

## 7. Implementation Risks and Mitigations

### Risks
- **Breaking changes**: Switching to JWT-based authentication may break existing code
- **Performance overhead**: JWT validation adds some overhead
- **Secret management**: JWT secret must be properly secured
- **Complex permissions**: Scopes may become overly complex to manage
- **Token leakage**: Service role tokens could be leaked

### Mitigations
- Implement the changes incrementally with backward compatibility
- Optimize JWT validation for performance
- Use proper secret management solutions
- Keep permission model simple and document well
- Use short expiration times and proper revocation for tokens

## 8. Integration with Application Code

### API Layer Updates
- Update API resources to use service_role_operation decorator
- Update API documentation to reflect new authentication
- Implement proper error handling for service role operations

### Service Updates
- Update background tasks to use service role client
- Update scheduled jobs to use service role client
- Update data processing to use service role client

## 9. Conclusion

Implementing JWT-based service role authentication with proper token management, scope-based permissions, and comprehensive audit logging will significantly improve the security and maintainability of the CryoProtect system. By creating a clear separation between direct service role key usage and application-level permissions, we can reduce security risks while improving operational visibility.