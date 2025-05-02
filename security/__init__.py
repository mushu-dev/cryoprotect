"""
Security package for CryoProtect v2.

This package provides security-related functionality including:
- Encryption at rest for sensitive data
- Key management
- Secure hashing
"""

from security.encryption import (
    EncryptionService,
    get_encryption_service,
    encrypt_fields,
    decrypt_fields
)

from security.encryption_config import (
    get_encrypted_fields,
    get_hashed_fields,
    is_sensitive_table
)

__all__ = [
    'EncryptionService',
    'get_encryption_service',
    'encrypt_fields',
    'decrypt_fields',
    'get_encrypted_fields',
    'get_hashed_fields',
    'is_sensitive_table'
]