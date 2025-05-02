# Encryption at Rest Implementation Guide

This document provides comprehensive documentation for the encryption at rest implementation in CryoProtect v2. It covers the encryption architecture, key management, usage patterns, and security considerations.

## Table of Contents

1. [Overview](#overview)
2. [Architecture](#architecture)
3. [Key Management](#key-management)
4. [Usage Guide](#usage-guide)
5. [Field-Level Encryption](#field-level-encryption)
6. [File Encryption](#file-encryption)
7. [Backup Encryption](#backup-encryption)
8. [Key Rotation](#key-rotation)
9. [Security Considerations](#security-considerations)
10. [Troubleshooting](#troubleshooting)

## Overview

The encryption at rest implementation in CryoProtect v2 provides protection for sensitive data stored in the database, files, and backups. It uses the Fernet symmetric encryption scheme from the Python cryptography library, which provides authenticated encryption with AES-128 in CBC mode and PKCS7 padding, with HMAC using SHA256 for authentication.

Key features:
- Field-level encryption for sensitive database fields
- File encryption for exports and backups
- Key management with support for key rotation
- Secure hashing for comparison of sensitive data
- Integration with existing backup systems

## Architecture

The encryption implementation consists of the following components:

1. **EncryptionService**: The core service that provides encryption, decryption, and key management functionality.
2. **Encryption Configuration**: Defines which fields in each model should be encrypted.
3. **Field Encryption Decorators**: Decorators for model methods to automatically encrypt/decrypt fields.
4. **Backup Integration**: Integration with the existing backup system to encrypt backup files.
5. **Migration Utility**: A utility to encrypt existing data in the database.

The encryption service is implemented as a singleton to ensure consistent key usage throughout the application.

## Key Management

Encryption keys are stored in the `config/keys` directory by default. The following files are used:

- `primary.key`: The current primary encryption key
- `key_metadata.json`: Metadata about all keys, including creation date, status, and purpose

Keys are automatically generated if they don't exist when the encryption service is initialized. The service supports multiple keys to enable key rotation without losing the ability to decrypt data encrypted with older keys.

### Key Security

The encryption keys are critical security assets and must be protected:

1. Keys should be stored in a secure location with restricted access
2. In production, consider using a key management service (KMS) like AWS KMS, HashiCorp Vault, or Azure Key Vault
3. Keys should be backed up securely to prevent data loss
4. Key files should have restricted permissions (0600)

## Usage Guide

### Basic Usage

```python
from security.encryption import get_encryption_service

# Get the encryption service
encryption_service = get_encryption_service()

# Encrypt data
sensitive_data = "This is sensitive information"
encrypted_data = encryption_service.encrypt(sensitive_data)

# Decrypt data
decrypted_data = encryption_service.decrypt(encrypted_data)
```

### Encrypting Dictionary Values

```python
# Encrypt specific fields in a dictionary
user_data = {
    "id": "user123",
    "name": "John Doe",
    "email": "john@example.com",
    "role": "admin"
}

fields_to_encrypt = ["name", "email"]
encrypted_data = encryption_service.encrypt_dict_values(user_data, fields_to_encrypt)

# Decrypt specific fields
decrypted_data = encryption_service.decrypt_dict_values(encrypted_data, fields_to_encrypt)
```

## Field-Level Encryption

The system provides decorators to automatically encrypt and decrypt fields in model methods:

```python
from security.encryption import encrypt_fields, decrypt_fields
from security.encryption_config import get_encrypted_fields

class UserProfile:
    @encrypt_fields(['email', 'name'])
    def create(self, data):
        # The email and name fields will be automatically encrypted
        return supabase.table('user_profile').insert(data).execute()
    
    @decrypt_fields(['email', 'name'])
    def get(self, user_id):
        # The email and name fields will be automatically decrypted
        return supabase.table('user_profile').select('*').eq('user_id', user_id).execute()
```

You can also use the configuration-based approach:

```python
@encrypt_fields(get_encrypted_fields('user_profile'))
def create_user_profile(data):
    # Fields defined in encryption_config.py will be encrypted
    return supabase.table('user_profile').insert(data).execute()
```

## File Encryption

The encryption service provides methods to encrypt and decrypt files:

```python
# Encrypt a file
encrypted_file_path = encryption_service.encrypt_file('path/to/sensitive/file.txt')

# Decrypt a file
decrypted_file_path = encryption_service.decrypt_file('path/to/encrypted/file.txt.enc')
```

## Backup Encryption

The encryption service is integrated with the existing backup system. Backups are automatically encrypted if encryption is enabled in the backup configuration:

```python
# In backup_config.json
{
    "encryption": {
        "enabled": true
    }
}
```

The `database_backup.py` script uses the encryption service to encrypt backup files before storing them.

## Key Rotation

Key rotation is an important security practice to limit the impact of key compromise. The encryption service supports key rotation:

```python
# Rotate the primary encryption key
rotation_info = encryption_service.rotate_key()

# The rotation_info contains information about the old and new keys
# This can be used to re-encrypt data with the new key
```

After rotating the key, you should re-encrypt sensitive data with the new key. The `apply_encryption.py` script can be used for this purpose:

```bash
python security/apply_encryption.py
```

## Security Considerations

1. **Key Protection**: Protect encryption keys as you would protect passwords or other credentials.
2. **Regular Key Rotation**: Rotate encryption keys regularly (e.g., every 90 days).
3. **Secure Configuration**: Ensure that the encryption configuration correctly identifies all sensitive fields.
4. **Backup Security**: Ensure that encrypted backups are stored securely and that keys are backed up separately.
5. **Access Control**: Limit access to decryption functionality to authorized users and services.
6. **Logging**: Monitor and log access to sensitive data and encryption operations.

## Troubleshooting

### Common Issues

1. **InvalidToken Error**: This occurs when trying to decrypt data with the wrong key or when the data is corrupted. Ensure you're using the correct key and that the data hasn't been modified.

2. **Missing Keys**: If encryption keys are missing, the service will generate new ones by default. This can lead to data loss if the original keys are needed to decrypt existing data. Always back up keys securely.

3. **Performance Impact**: Field-level encryption can impact query performance, especially for fields used in WHERE clauses. Consider using database-level encryption for large datasets.

### Debugging

The encryption service logs operations to the application logs. You can increase the log level for more detailed information:

```python
import logging
logging.getLogger('security.encryption').setLevel(logging.DEBUG)
```

### Testing

A comprehensive test suite is provided in `tests/test_encryption.py`. Run these tests to verify that the encryption system is working correctly:

```bash
python tests/test_encryption.py
```

---

This guide provides an overview of the encryption at rest implementation in CryoProtect v2. For more detailed information, refer to the source code and comments in the `security` package.