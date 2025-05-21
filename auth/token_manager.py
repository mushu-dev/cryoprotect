#!/usr/bin/env python3
"""
Token manager for JWT-based service role authentication.

This module handles the creation, validation, and revocation of JWT tokens
for service-to-service communication in CryoProtect.
"""

import os
import time
import uuid
import json
import logging
import hashlib
import threading
from typing import Dict, List, Optional, Any, Union, Tuple
from datetime import datetime, timedelta

import jwt
from cryptography.hazmat.primitives import serialization
from cryptography.hazmat.primitives.asymmetric import rsa
from cryptography.hazmat.backends import default_backend

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class TokenManager:
    """
    Manages JWT tokens for service role authentication.
    
    Features:
    - Token creation with configurable expiration
    - Token validation with signature and claims verification
    - Token revocation list for invalidating tokens
    - Auto-rotation of signing keys
    - Thread-safe operations
    """
    
    def __init__(
        self,
        secret_key: Optional[str] = None,
        token_expiration: int = 3600,  # 1 hour
        algorithm: str = 'HS256',
        issuer: str = 'cryoprotect-auth',
        auto_rotate_keys: bool = True,
        key_rotation_interval: int = 86400,  # 24 hours
        clean_interval: int = 3600  # 1 hour
    ):
        """
        Initialize the token manager.
        
        Args:
            secret_key: Secret key for token signing (if None, a key will be generated)
            token_expiration: Token expiration time in seconds
            algorithm: JWT signing algorithm (HS256, HS384, HS512, RS256, etc.)
            issuer: Token issuer identifier
            auto_rotate_keys: Whether to automatically rotate keys
            key_rotation_interval: Interval in seconds for key rotation
            clean_interval: Interval in seconds for cleaning revoked tokens
        """
        self.token_expiration = token_expiration
        self.algorithm = algorithm
        self.issuer = issuer
        self.auto_rotate_keys = auto_rotate_keys
        self.key_rotation_interval = key_rotation_interval
        self.clean_interval = clean_interval
        
        # Use RSA algorithm if specified
        self.use_rsa = algorithm.startswith('RS')
        
        # Initialize key storage
        self.current_key_id = None
        self.keys = {}
        self.revoked_tokens = {}
        self.lock = threading.RLock()
        
        # Initialize keys
        if self.use_rsa:
            self._generate_rsa_key()
        else:
            # Use provided secret key or generate one
            if secret_key:
                self._set_hmac_key(secret_key)
            else:
                self._generate_hmac_key()
        
        # Start maintenance tasks if auto-rotation is enabled
        if auto_rotate_keys:
            self._start_maintenance_tasks()
    
    def _generate_hmac_key(self):
        """Generate a new HMAC secret key."""
        key_id = str(uuid.uuid4())
        secret = os.urandom(32).hex()
        
        with self.lock:
            self.keys[key_id] = {
                'secret': secret,
                'created_at': time.time(),
                'type': 'hmac'
            }
            self.current_key_id = key_id
        
        logger.info(f"Generated new HMAC key with ID: {key_id}")
    
    def _set_hmac_key(self, secret: str):
        """Set an existing HMAC secret key."""
        key_id = str(uuid.uuid4())
        
        with self.lock:
            self.keys[key_id] = {
                'secret': secret,
                'created_at': time.time(),
                'type': 'hmac'
            }
            self.current_key_id = key_id
        
        logger.info(f"Set HMAC key with ID: {key_id}")
    
    def _generate_rsa_key(self):
        """Generate a new RSA key pair."""
        key_id = str(uuid.uuid4())
        
        # Generate private key
        private_key = rsa.generate_private_key(
            public_exponent=65537,
            key_size=2048,
            backend=default_backend()
        )
        
        # Get public key
        public_key = private_key.public_key()
        
        # Serialize keys
        private_pem = private_key.private_bytes(
            encoding=serialization.Encoding.PEM,
            format=serialization.PrivateFormat.PKCS8,
            encryption_algorithm=serialization.NoEncryption()
        ).decode('utf-8')
        
        public_pem = public_key.public_bytes(
            encoding=serialization.Encoding.PEM,
            format=serialization.PublicFormat.SubjectPublicKeyInfo
        ).decode('utf-8')
        
        with self.lock:
            self.keys[key_id] = {
                'private_key': private_pem,
                'public_key': public_pem,
                'created_at': time.time(),
                'type': 'rsa'
            }
            self.current_key_id = key_id
        
        logger.info(f"Generated new RSA key pair with ID: {key_id}")
    
    def _get_current_key(self) -> Tuple[str, Dict[str, Any]]:
        """Get the current signing key."""
        with self.lock:
            return self.current_key_id, self.keys[self.current_key_id]
    
    def _get_key_by_id(self, key_id: str) -> Optional[Dict[str, Any]]:
        """Get a key by its ID."""
        with self.lock:
            return self.keys.get(key_id)
    
    def _start_maintenance_tasks(self):
        """Start background maintenance tasks."""
        # Key rotation thread
        rotation_thread = threading.Thread(
            target=self._key_rotation_task,
            daemon=True,
            name="token-key-rotation"
        )
        rotation_thread.start()
        
        # Token cleanup thread
        cleanup_thread = threading.Thread(
            target=self._token_cleanup_task,
            daemon=True,
            name="token-cleanup"
        )
        cleanup_thread.start()
    
    def _key_rotation_task(self):
        """Background task for key rotation."""
        while True:
            try:
                time.sleep(60)  # Check every minute
                
                current_key_id, key = self._get_current_key()
                key_age = time.time() - key['created_at']
                
                if key_age >= self.key_rotation_interval:
                    logger.info(f"Rotating key {current_key_id} after {key_age:.2f} seconds")
                    
                    # Generate new key based on algorithm
                    if self.use_rsa:
                        self._generate_rsa_key()
                    else:
                        self._generate_hmac_key()
                    
                    # Keep old keys for a while to validate existing tokens
                    # Remove keys older than 2x the rotation interval
                    self._clean_old_keys()
            except Exception as e:
                logger.error(f"Error in key rotation task: {str(e)}")
    
    def _token_cleanup_task(self):
        """Background task for cleaning up expired revoked tokens."""
        while True:
            try:
                time.sleep(60)  # Check every minute
                
                now = time.time()
                if now - getattr(self, '_last_cleanup', 0) >= self.clean_interval:
                    self._clean_revoked_tokens()
                    self._last_cleanup = now
            except Exception as e:
                logger.error(f"Error in token cleanup task: {str(e)}")
    
    def _clean_old_keys(self):
        """Clean up old keys."""
        now = time.time()
        keys_to_remove = []
        
        with self.lock:
            # Keep current key
            current_key_id = self.current_key_id
            
            # Find keys older than 2x rotation interval
            for key_id, key in self.keys.items():
                if (key_id != current_key_id and 
                    now - key['created_at'] > 2 * self.key_rotation_interval):
                    keys_to_remove.append(key_id)
            
            # Remove old keys
            for key_id in keys_to_remove:
                del self.keys[key_id]
                logger.info(f"Removed old key: {key_id}")
    
    def _clean_revoked_tokens(self):
        """Clean up expired revoked tokens."""
        now = time.time()
        tokens_to_remove = []
        
        with self.lock:
            for token_id, expiry in self.revoked_tokens.items():
                if now >= expiry:
                    tokens_to_remove.append(token_id)
            
            for token_id in tokens_to_remove:
                del self.revoked_tokens[token_id]
            
            if tokens_to_remove:
                logger.info(f"Cleaned up {len(tokens_to_remove)} expired revoked tokens")
    
    def create_token(
        self,
        subject: str,
        scopes: List[str] = None,
        expiration: Optional[int] = None,
        custom_claims: Dict[str, Any] = None
    ) -> str:
        """
        Create a new JWT token.
        
        Args:
            subject: Token subject (service name, user ID, etc.)
            scopes: List of permission scopes
            expiration: Custom expiration time in seconds (overrides default)
            custom_claims: Additional custom claims to include
            
        Returns:
            JWT token string
        """
        # Get current key
        key_id, key = self._get_current_key()
        
        # Prepare token payload
        now = int(time.time())
        exp = now + (expiration or self.token_expiration)
        
        payload = {
            'iss': self.issuer,
            'sub': subject,
            'iat': now,
            'exp': exp,
            'jti': str(uuid.uuid4()),
            'kid': key_id
        }
        
        # Add scopes if provided
        if scopes:
            payload['scopes'] = scopes
        
        # Add custom claims if provided
        if custom_claims:
            for claim, value in custom_claims.items():
                if claim not in payload:
                    payload[claim] = value
        
        # Sign token
        if self.use_rsa:
            token = jwt.encode(
                payload,
                key['private_key'],
                algorithm=self.algorithm
            )
        else:
            token = jwt.encode(
                payload,
                key['secret'],
                algorithm=self.algorithm
            )
        
        # Convert bytes to string if needed
        if isinstance(token, bytes):
            token = token.decode('utf-8')
        
        return token
    
    def validate_token(self, token: str) -> Dict[str, Any]:
        """
        Validate a JWT token.
        
        Args:
            token: JWT token string
            
        Returns:
            Dict with token claims
            
        Raises:
            jwt.InvalidTokenError: If the token is invalid
        """
        # Decode token header to get key ID
        try:
            header = jwt.get_unverified_header(token)
            key_id = header.get('kid')
            
            if not key_id:
                raise jwt.InvalidTokenError("Token missing key ID")
            
            # Get key by ID
            key = self._get_key_by_id(key_id)
            if not key:
                raise jwt.InvalidTokenError("Unknown key ID")
            
            # Check algorithm
            token_alg = header.get('alg')
            if token_alg != self.algorithm:
                raise jwt.InvalidTokenError(f"Invalid algorithm: {token_alg}")
            
            # Check if token is revoked
            # Get the token ID without verifying the token
            unverified_claims = jwt.decode(
                token,
                options={"verify_signature": False}
            )
            token_id = unverified_claims.get('jti')
            
            with self.lock:
                if token_id in self.revoked_tokens:
                    raise jwt.InvalidTokenError("Token has been revoked")
            
            # Verify token
            if key['type'] == 'rsa':
                claims = jwt.decode(
                    token,
                    key['public_key'],
                    algorithms=[self.algorithm],
                    issuer=self.issuer
                )
            else:
                claims = jwt.decode(
                    token,
                    key['secret'],
                    algorithms=[self.algorithm],
                    issuer=self.issuer
                )
            
            return claims
        
        except jwt.PyJWTError as e:
            logger.warning(f"Token validation failed: {str(e)}")
            raise
    
    def revoke_token(self, token: Union[str, Dict[str, Any]]):
        """
        Revoke a JWT token.
        
        Args:
            token: JWT token string or decoded claims
            
        Raises:
            ValueError: If the token is invalid or missing required claims
        """
        try:
            # Get token claims if a token string is provided
            claims = token if isinstance(token, dict) else self.validate_token(token)
            
            # Get token ID and expiration
            token_id = claims.get('jti')
            expiration = claims.get('exp')
            
            if not token_id or not expiration:
                raise ValueError("Token missing required claims (jti, exp)")
            
            # Add to revoked tokens
            with self.lock:
                self.revoked_tokens[token_id] = expiration
            
            logger.info(f"Revoked token with ID: {token_id}")
        
        except (jwt.PyJWTError, ValueError) as e:
            logger.warning(f"Token revocation failed: {str(e)}")
            raise
    
    def get_public_jwks(self) -> Dict[str, Any]:
        """
        Get public keys in JWKS format for clients to verify tokens.
        Only applicable for RSA keys.
        
        Returns:
            Dict with JWKS keys
        """
        if not self.use_rsa:
            raise ValueError("JWKS is only available for RSA keys")
        
        jwks = {
            'keys': []
        }
        
        with self.lock:
            for key_id, key_data in self.keys.items():
                if key_data['type'] == 'rsa':
                    # Extract components for JWK
                    # This is a simplified version - real implementation should
                    # extract n and e values from the public key
                    jwks['keys'].append({
                        'kid': key_id,
                        'kty': 'RSA',
                        'alg': self.algorithm,
                        'use': 'sig',
                        # Provide public key in PEM format
                        'x5c': [key_data['public_key']]
                    })
        
        return jwks
    
    def get_stats(self) -> Dict[str, Any]:
        """
        Get statistics about the token manager.
        
        Returns:
            Dict with statistics
        """
        with self.lock:
            return {
                'active_keys': len(self.keys),
                'revoked_tokens': len(self.revoked_tokens),
                'algorithm': self.algorithm,
                'token_expiration': self.token_expiration,
                'auto_rotate_keys': self.auto_rotate_keys,
                'key_rotation_interval': self.key_rotation_interval
            }

# Global token manager instance
_token_manager_instance = None

def get_token_manager(
    secret_key: Optional[str] = None,
    token_expiration: Optional[int] = None,
    algorithm: Optional[str] = None
) -> TokenManager:
    """
    Get the global token manager instance.
    
    Args:
        secret_key: Optional secret key for token signing
        token_expiration: Optional token expiration time in seconds
        algorithm: Optional JWT signing algorithm
        
    Returns:
        TokenManager instance
    """
    global _token_manager_instance
    
    if _token_manager_instance is None:
        # Get configuration from environment variables
        env_secret = os.environ.get('JWT_SECRET_KEY')
        env_expiration = os.environ.get('JWT_EXPIRATION')
        env_algorithm = os.environ.get('JWT_ALGORITHM', 'HS256')
        env_issuer = os.environ.get('JWT_ISSUER', 'cryoprotect-auth')
        
        # Use parameters or environment variables
        secret = secret_key or env_secret
        expiration = token_expiration or (int(env_expiration) if env_expiration else 3600)
        alg = algorithm or env_algorithm
        
        # Create token manager
        _token_manager_instance = TokenManager(
            secret_key=secret,
            token_expiration=expiration,
            algorithm=alg,
            issuer=env_issuer
        )
    
    return _token_manager_instance

# Convenience functions
def create_service_token(
    service_name: str,
    scopes: List[str] = None,
    expiration: Optional[int] = None
) -> str:
    """
    Create a service-to-service authentication token.
    
    Args:
        service_name: Name of the service
        scopes: List of permission scopes
        expiration: Custom expiration time in seconds
        
    Returns:
        JWT token string
    """
    token_manager = get_token_manager()
    
    # Add service role claim
    custom_claims = {
        'role': 'service_role',
        'client': service_name
    }
    
    return token_manager.create_token(
        subject=service_name,
        scopes=scopes,
        expiration=expiration,
        custom_claims=custom_claims
    )

def validate_service_token(token: str) -> Dict[str, Any]:
    """
    Validate a service-to-service authentication token.
    
    Args:
        token: JWT token string
        
    Returns:
        Dict with token claims
        
    Raises:
        jwt.InvalidTokenError: If the token is invalid
        ValueError: If the token is not a service token
    """
    token_manager = get_token_manager()
    
    # Validate token
    claims = token_manager.validate_token(token)
    
    # Check service role
    if claims.get('role') != 'service_role':
        raise ValueError("Token is not a service token")
    
    return claims

def revoke_service_token(token: Union[str, Dict[str, Any]]):
    """
    Revoke a service-to-service authentication token.
    
    Args:
        token: JWT token string or decoded claims
        
    Raises:
        ValueError: If the token is invalid or missing required claims
    """
    token_manager = get_token_manager()
    token_manager.revoke_token(token)

# Main function for testing
def main():
    """Main function for command-line usage."""
    import argparse
    
    parser = argparse.ArgumentParser(description='JWT Token Manager')
    parser.add_argument('--create', action='store_true', help='Create a new token')
    parser.add_argument('--validate', help='Validate a token')
    parser.add_argument('--revoke', help='Revoke a token')
    parser.add_argument('--service', default='test-service', help='Service name for token creation')
    parser.add_argument('--scopes', help='Comma-separated list of scopes')
    parser.add_argument('--expiration', type=int, help='Token expiration in seconds')
    parser.add_argument('--secret', help='Secret key for token signing')
    parser.add_argument('--algorithm', default='HS256', help='JWT algorithm (HS256, RS256, etc.)')
    parser.add_argument('--stats', action='store_true', help='Show token manager stats')
    
    args = parser.parse_args()
    
    # Create token manager
    token_manager = get_token_manager(
        secret_key=args.secret,
        token_expiration=args.expiration,
        algorithm=args.algorithm
    )
    
    if args.create:
        scopes = args.scopes.split(',') if args.scopes else None
        token = create_service_token(
            service_name=args.service,
            scopes=scopes,
            expiration=args.expiration
        )
        print(f"Token: {token}")
    
    elif args.validate:
        try:
            claims = validate_service_token(args.validate)
            print(f"Token is valid: {json.dumps(claims, indent=2)}")
        except Exception as e:
            print(f"Token validation failed: {str(e)}")
    
    elif args.revoke:
        try:
            revoke_service_token(args.revoke)
            print(f"Token revoked successfully")
        except Exception as e:
            print(f"Token revocation failed: {str(e)}")
    
    elif args.stats:
        stats = token_manager.get_stats()
        print(f"Token Manager Stats: {json.dumps(stats, indent=2)}")
    
    else:
        parser.print_help()

if __name__ == "__main__":
    main()