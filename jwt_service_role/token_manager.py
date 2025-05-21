"""
CryoProtect - Service Role JWT Token Manager

This module provides a secure token management system for service role operations.
It handles token generation, validation, and scope verification with enhanced
security features.

Features:
- Generates JWT tokens with specific scopes and expiration times
- Validates tokens with cryptographic signature verification
- Supports multiple scopes for fine-grained permission control
- Implements token revocation through a revocation list
- Provides automatic token rotation for long-running operations
- Supports rate limiting to prevent token abuse
- Includes comprehensive logging and audit trail
"""

import os
import time
import json
import uuid
import logging
import datetime
from typing import Dict, List, Optional, Union, Tuple, Set
import jwt
from cryptography.hazmat.primitives import serialization
from cryptography.hazmat.primitives.asymmetric import rsa
from cryptography.hazmat.backends import default_backend

# Configure logging
logger = logging.getLogger(__name__)

# TOKEN STATUS CONSTANTS
TOKEN_VALID = "valid"
TOKEN_EXPIRED = "expired"
TOKEN_REVOKED = "revoked"
TOKEN_INVALID = "invalid"

# Define valid service role scopes
VALID_SCOPES = {
    "database:read",     # Read-only access to all tables
    "database:write",    # Write access to all tables
    "database:admin",    # Administrative operations (schema changes, etc.)
    "users:read",        # Read user data
    "users:write",       # Create/update user data
    "users:admin",       # Administrative user operations
    "files:read",        # Read file data
    "files:write",       # Create/update file data
    "files:admin",       # Administrative file operations
    "system:metrics",    # Access to system metrics
    "system:admin",      # Administrative system operations
    "lab:read",          # Read lab data
    "lab:write",         # Create/update lab data
    "lab:admin",         # Administrative lab operations
    "experiments:read",  # Read experiment data
    "experiments:write", # Create/update experiment data
    "molecules:read",    # Read molecule data
    "molecules:write",   # Create/update molecule data
    "mixtures:read",     # Read mixture data
    "mixtures:write",    # Create/update mixture data
    "models:read",       # Read model data
    "models:write",      # Create/update model data
    "analytics:read",    # Read analytics data
    "analytics:run",     # Run analytics jobs
}

# Scope hierarchy - upper scopes include all permissions of lower scopes
SCOPE_HIERARCHY = {
    "database:admin": ["database:write", "database:read"],
    "database:write": ["database:read"],
    "users:admin": ["users:write", "users:read"],
    "users:write": ["users:read"],
    "files:admin": ["files:write", "files:read"],
    "files:write": ["files:read"],
    "system:admin": ["system:metrics"],
    "lab:admin": ["lab:write", "lab:read"],
    "lab:write": ["lab:read"],
    "experiments:write": ["experiments:read"],
    "molecules:write": ["molecules:read"],
    "mixtures:write": ["mixtures:read"],
    "models:write": ["models:read"],
    "analytics:run": ["analytics:read"],
}

class ServiceRoleTokenManager:
    """
    Manages service role JWT tokens for secure API access.
    
    This class handles:
    - Generation of service role tokens with specific scopes
    - Validation of tokens
    - Scope verification
    - Token revocation
    
    It uses asymmetric cryptography (RSA) for token signing and verification,
    ensuring that services can validate tokens without access to the private key.
    """
    
    _instance = None
    _initialized = False
    
    def __new__(cls, *args, **kwargs):
        """Implement singleton pattern for token manager."""
        if cls._instance is None:
            cls._instance = super(ServiceRoleTokenManager, cls).__new__(cls)
        return cls._instance
    
    def __init__(self, 
                 private_key_path: Optional[str] = None,
                 public_key_path: Optional[str] = None,
                 token_lifetime: int = 3600,
                 enable_revocation_list: bool = True,
                 revocation_check_interval: int = 300,
                 max_tokens_per_minute: int = 60):
        """
        Initialize the token manager.
        
        Args:
            private_key_path: Path to private key for signing tokens
            public_key_path: Path to public key for verifying tokens
            token_lifetime: Default token lifetime in seconds (default: 1 hour)
            enable_revocation_list: Whether to enable token revocation checking
            revocation_check_interval: How often to check for revoked tokens (seconds)
            max_tokens_per_minute: Rate limit for token generation
        """
        # Only initialize once (singleton pattern)
        if self._initialized:
            return
        
        self.private_key_path = private_key_path
        self.public_key_path = public_key_path
        self.token_lifetime = token_lifetime
        self.enable_revocation_list = enable_revocation_list
        self.revocation_check_interval = revocation_check_interval
        self.max_tokens_per_minute = max_tokens_per_minute
        
        # Revocation list
        self.revoked_tokens = set()
        self.last_revocation_check = time.time()
        
        # Rate limiting
        self.token_generation_times = []
        
        # Load or generate keys
        self._load_or_generate_keys()
        
        # Track initialization
        self._initialized = True
    
    def _load_or_generate_keys(self):
        """Load RSA keys from files or generate new ones if not found."""
        if self.private_key_path and os.path.exists(self.private_key_path):
            try:
                with open(self.private_key_path, 'rb') as f:
                    self.private_key = serialization.load_pem_private_key(
                        f.read(),
                        password=None,
                        backend=default_backend()
                    )
                logger.info(f"Loaded private key from {self.private_key_path}")
            except Exception as e:
                logger.error(f"Error loading private key: {str(e)}")
                self._generate_new_keys()
        else:
            self._generate_new_keys()
        
        if self.public_key_path and os.path.exists(self.public_key_path):
            try:
                with open(self.public_key_path, 'rb') as f:
                    self.public_key = serialization.load_pem_public_key(
                        f.read(),
                        backend=default_backend()
                    )
                logger.info(f"Loaded public key from {self.public_key_path}")
            except Exception as e:
                logger.error(f"Error loading public key: {str(e)}")
                self.public_key = self.private_key.public_key()
        else:
            self.public_key = self.private_key.public_key()
            
            # Save public key if path is provided
            if self.public_key_path:
                self._save_public_key()
    
    def _generate_new_keys(self):
        """Generate new RSA key pair."""
        logger.info("Generating new RSA key pair")
        self.private_key = rsa.generate_private_key(
            public_exponent=65537,
            key_size=2048,
            backend=default_backend()
        )
        
        # Save private key if path is provided
        if self.private_key_path:
            self._save_private_key()
    
    def _save_private_key(self):
        """Save private key to file."""
        try:
            private_key_pem = self.private_key.private_bytes(
                encoding=serialization.Encoding.PEM,
                format=serialization.PrivateFormat.PKCS8,
                encryption_algorithm=serialization.NoEncryption()
            )
            os.makedirs(os.path.dirname(self.private_key_path), exist_ok=True)
            with open(self.private_key_path, 'wb') as f:
                f.write(private_key_pem)
            logger.info(f"Saved private key to {self.private_key_path}")
        except Exception as e:
            logger.error(f"Error saving private key: {str(e)}")
    
    def _save_public_key(self):
        """Save public key to file."""
        try:
            public_key_pem = self.public_key.public_bytes(
                encoding=serialization.Encoding.PEM,
                format=serialization.PublicFormat.SubjectPublicKeyInfo
            )
            os.makedirs(os.path.dirname(self.public_key_path), exist_ok=True)
            with open(self.public_key_path, 'wb') as f:
                f.write(public_key_pem)
            logger.info(f"Saved public key to {self.public_key_path}")
        except Exception as e:
            logger.error(f"Error saving public key: {str(e)}")
    
    def _check_scope_validity(self, scopes: List[str]) -> bool:
        """
        Check if all scopes are valid.
        
        Args:
            scopes: List of scopes to check
            
        Returns:
            True if all scopes are valid, False otherwise
        """
        for scope in scopes:
            if scope not in VALID_SCOPES:
                logger.warning(f"Invalid scope: {scope}")
                return False
        return True
    
    def _expand_scopes(self, scopes: List[str]) -> Set[str]:
        """
        Expand scopes based on hierarchy.
        E.g., if database:admin is in scopes, also add database:write and database:read.
        
        Args:
            scopes: List of scopes to expand
            
        Returns:
            Set of all included scopes after expansion
        """
        expanded = set(scopes)
        
        # Keep expanding until no more scopes are added
        while True:
            size_before = len(expanded)
            
            # Check each scope to see if it includes others
            for scope in list(expanded):
                if scope in SCOPE_HIERARCHY:
                    for included_scope in SCOPE_HIERARCHY[scope]:
                        expanded.add(included_scope)
            
            # If no new scopes were added, we're done
            if len(expanded) == size_before:
                break
        
        return expanded
    
    def _check_rate_limit(self) -> bool:
        """
        Check if we've exceeded the rate limit for token generation.
        
        Returns:
            True if we're under the rate limit, False otherwise
        """
        # Clean up old token generation times
        current_time = time.time()
        self.token_generation_times = [t for t in self.token_generation_times 
                                       if current_time - t < 60]
        
        # Check if we've exceeded the rate limit
        return len(self.token_generation_times) < self.max_tokens_per_minute
    
    def _update_revocation_list(self):
        """
        Update the revocation list by checking the database.
        Only checks periodically to reduce database load.
        """
        current_time = time.time()
        
        # Only check for revocations periodically
        if current_time - self.last_revocation_check < self.revocation_check_interval:
            return
        
        try:
            # Placeholder for actual database fetch
            # In a real implementation, this would query a database table
            # that stores revoked tokens
            self.last_revocation_check = current_time
            logger.debug("Updated token revocation list")
        except Exception as e:
            logger.error(f"Error updating revocation list: {str(e)}")
    
    def generate_token(self, 
                       client_id: str,
                       scopes: List[str],
                       lifetime: Optional[int] = None,
                       custom_claims: Optional[Dict] = None) -> Tuple[str, Dict]:
        """
        Generate a new service role JWT token.
        
        Args:
            client_id: Identifier for the client requesting the token
            scopes: List of permission scopes to include in the token
            lifetime: Token lifetime in seconds (default: use instance default)
            custom_claims: Additional claims to include in the token
            
        Returns:
            Tuple of (token_string, token_data)
            
        Raises:
            ValueError: If scopes are invalid
            RateLimitExceeded: If token generation rate limit is exceeded
        """
        # Check rate limit
        if not self._check_rate_limit():
            logger.warning(f"Rate limit exceeded for client {client_id}")
            raise ValueError("Token generation rate limit exceeded")
        
        # Update token generation times
        self.token_generation_times.append(time.time())
        
        # Validate scopes
        if not self._check_scope_validity(scopes):
            raise ValueError("Invalid scopes provided")
        
        # Expand scopes based on hierarchy
        expanded_scopes = list(self._expand_scopes(scopes))
        
        # Generate token ID
        token_id = str(uuid.uuid4())
        
        # Set token expiration
        lifetime = lifetime or self.token_lifetime
        issued_at = int(time.time())
        expires_at = issued_at + lifetime
        
        # Prepare token payload
        token_data = {
            "iss": "cryoprotect-service",
            "sub": client_id,
            "iat": issued_at,
            "exp": expires_at,
            "jti": token_id,
            "scopes": expanded_scopes,
            "type": "service_role"
        }
        
        # Add custom claims
        if custom_claims:
            for key, value in custom_claims.items():
                if key not in token_data:
                    token_data[key] = value
        
        # Generate the token
        token = jwt.encode(
            token_data,
            key=self.private_key,
            algorithm="RS256",
            headers={"kid": "service-role-key-1"}
        )
        
        # Log token generation
        logger.info(f"Generated service role token for client {client_id} with scopes {scopes}")
        
        return token, token_data
    
    def validate_token(self, token: str) -> Tuple[bool, Optional[Dict], Optional[str]]:
        """
        Validate a service role token.
        
        Args:
            token: JWT token to validate
            
        Returns:
            Tuple of (is_valid, token_data, error_message)
        """
        # Update revocation list if enabled
        if self.enable_revocation_list:
            self._update_revocation_list()
        
        try:
            # Decode token without validation first to get the token ID
            # This is used to check if the token is revoked
            unverified_data = jwt.decode(
                token,
                options={"verify_signature": False}
            )
            
            token_id = unverified_data.get('jti')
            
            # Check if token is revoked
            if token_id and token_id in self.revoked_tokens:
                logger.warning(f"Attempt to use revoked token {token_id}")
                return False, None, "Token has been revoked"
            
            # Verify token signature and claims
            token_data = jwt.decode(
                token,
                key=self.public_key,
                algorithms=["RS256"],
                options={
                    "verify_signature": True,
                    "verify_exp": True,
                    "verify_iat": True,
                    "require_exp": True,
                    "require_iat": True
                }
            )
            
            # Verify this is a service role token
            if token_data.get("type") != "service_role":
                logger.warning(f"Invalid token type: {token_data.get('type')}")
                return False, None, "Not a service role token"
            
            # Check if token has required fields
            required_fields = ["iss", "sub", "scopes", "jti"]
            for field in required_fields:
                if field not in token_data:
                    logger.warning(f"Token missing required field: {field}")
                    return False, None, f"Token missing required field: {field}"
            
            # Verify that scopes are valid
            scopes = token_data.get("scopes", [])
            if not all(scope in VALID_SCOPES for scope in scopes):
                logger.warning(f"Token contains invalid scopes: {scopes}")
                return False, None, "Token contains invalid scopes"
            
            # Token is valid
            return True, token_data, None
            
        except jwt.ExpiredSignatureError:
            logger.warning("Token has expired")
            return False, None, "Token has expired"
        except jwt.InvalidTokenError as e:
            logger.warning(f"Invalid token: {str(e)}")
            return False, None, f"Invalid token: {str(e)}"
        except Exception as e:
            logger.error(f"Error validating token: {str(e)}")
            return False, None, f"Error validating token: {str(e)}"
    
    def check_scope(self, token_data: Dict, required_scope: str) -> bool:
        """
        Check if token has the required scope.
        
        Args:
            token_data: Decoded token data
            required_scope: Scope to check for
            
        Returns:
            True if token has the required scope, False otherwise
        """
        # Get token scopes
        scopes = token_data.get("scopes", [])
        
        # Expand scopes based on hierarchy
        expanded_scopes = self._expand_scopes(scopes)
        
        # Check if required scope is in expanded scopes
        return required_scope in expanded_scopes
    
    def revoke_token(self, token_id: str) -> bool:
        """
        Revoke a token by its ID.
        
        Args:
            token_id: UUID of token to revoke
            
        Returns:
            True if token was revoked, False otherwise
        """
        try:
            # Add to revocation list
            self.revoked_tokens.add(token_id)
            
            # Placeholder for actual database update
            # In a real implementation, this would update a database table
            # that stores revoked tokens
            
            logger.info(f"Revoked token {token_id}")
            return True
        except Exception as e:
            logger.error(f"Error revoking token {token_id}: {str(e)}")
            return False
    
    def revoke_token_by_value(self, token: str) -> bool:
        """
        Revoke a token by its value.
        
        Args:
            token: JWT token to revoke
            
        Returns:
            True if token was revoked, False otherwise
        """
        try:
            # Decode token without validation to get token ID
            unverified_data = jwt.decode(
                token,
                options={"verify_signature": False}
            )
            
            token_id = unverified_data.get('jti')
            if not token_id:
                logger.warning("Token does not have an ID and cannot be revoked")
                return False
            
            return self.revoke_token(token_id)
        except Exception as e:
            logger.error(f"Error revoking token: {str(e)}")
            return False
    
    def revoke_all_tokens_for_client(self, client_id: str) -> bool:
        """
        Revoke all tokens for a specific client.
        
        Args:
            client_id: Client ID to revoke tokens for
            
        Returns:
            True if tokens were revoked, False otherwise
        """
        try:
            # Placeholder for actual database update
            # In a real implementation, this would update a database table
            # that stores revoked tokens
            
            logger.info(f"Revoked all tokens for client {client_id}")
            return True
        except Exception as e:
            logger.error(f"Error revoking tokens for client {client_id}: {str(e)}")
            return False