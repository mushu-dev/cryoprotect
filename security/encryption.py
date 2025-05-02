"""
Encryption Service for CryoProtect v2

This module provides encryption services for sensitive data at rest, including:
1. Field-level encryption for database fields
2. File encryption for backups and exports
3. Key management (loading, generating, and rotating keys)

The implementation uses Fernet symmetric encryption from the cryptography library,
which provides authenticated encryption with AES-128 in CBC mode and PKCS7 padding,
with HMAC using SHA256 for authentication.
"""

import os
import base64
import json
import logging
import hashlib
from typing import Any, Dict, Optional, Union, List, Tuple
from datetime import datetime, timedelta
from pathlib import Path

from cryptography.fernet import Fernet, InvalidToken
from cryptography.hazmat.primitives import hashes
from cryptography.hazmat.primitives.kdf.pbkdf2 import PBKDF2HMAC
from cryptography.hazmat.backends import default_backend

# Configure logging
logger = logging.getLogger(__name__)
handler = logging.StreamHandler()
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
handler.setFormatter(formatter)
logger.addHandler(handler)
logger.setLevel(logging.INFO)

# Default paths
DEFAULT_KEY_DIR = os.environ.get('ENCRYPTION_KEY_DIR', 'config/keys')
DEFAULT_PRIMARY_KEY_FILE = os.path.join(DEFAULT_KEY_DIR, 'primary.key')
DEFAULT_KEY_METADATA_FILE = os.path.join(DEFAULT_KEY_DIR, 'key_metadata.json')

class EncryptionService:
    """Service for handling data encryption/decryption and key management."""
    
    def __init__(self, 
                 primary_key_file: Optional[str] = None,
                 key_metadata_file: Optional[str] = None,
                 auto_create_keys: bool = True):
        """
        Initialize the encryption service.
        
        Args:
            primary_key_file: Path to the primary encryption key file
            key_metadata_file: Path to the key metadata file
            auto_create_keys: Whether to automatically create keys if they don't exist
        """
        self.primary_key_file = primary_key_file or DEFAULT_PRIMARY_KEY_FILE
        self.key_metadata_file = key_metadata_file or DEFAULT_KEY_METADATA_FILE
        self.key_dir = os.path.dirname(self.primary_key_file)
        
        # Ensure key directory exists
        os.makedirs(self.key_dir, exist_ok=True)
        
        # Load or generate keys
        if auto_create_keys:
            self.primary_key, self.key_id = self._load_or_generate_primary_key()
            self.fernet = Fernet(self.primary_key)
            self.all_keys = self._load_all_keys()
        else:
            if not os.path.exists(self.primary_key_file):
                raise FileNotFoundError(f"Primary key file not found: {self.primary_key_file}")
            self.primary_key, self.key_id = self._load_primary_key()
            self.fernet = Fernet(self.primary_key)
            self.all_keys = self._load_all_keys()
    
    def _load_or_generate_primary_key(self) -> Tuple[bytes, str]:
        """Load existing primary key or generate a new one."""
        if os.path.exists(self.primary_key_file):
            return self._load_primary_key()
        else:
            return self._generate_new_primary_key()
    
    def _load_primary_key(self) -> Tuple[bytes, str]:
        """Load the primary key from file."""
        try:
            with open(self.primary_key_file, 'rb') as f:
                key_data = f.read()
            
            # Extract key ID from filename or metadata
            key_id = self._get_key_id_from_metadata(os.path.basename(self.primary_key_file))
            
            logger.info(f"Loaded primary encryption key: {self.primary_key_file}")
            return key_data, key_id
        except Exception as e:
            logger.error(f"Error loading primary key: {str(e)}")
            raise
    
    def _generate_new_primary_key(self) -> Tuple[bytes, str]:
        """Generate a new primary key and save it to file."""
        try:
            key = Fernet.generate_key()
            key_id = self._generate_key_id()
            
            # Save the key
            with open(self.primary_key_file, 'wb') as f:
                f.write(key)
            
            # Update metadata
            self._update_key_metadata(key_id, self.primary_key_file, is_primary=True)
            
            logger.info(f"Generated new primary encryption key: {self.primary_key_file}")
            return key, key_id
        except Exception as e:
            logger.error(f"Error generating primary key: {str(e)}")
            raise
    
    def _load_all_keys(self) -> Dict[str, bytes]:
        """Load all encryption keys from the key directory."""
        keys = {}
        
        # Always add the primary key
        keys[self.key_id] = self.primary_key
        
        # Load metadata to find other keys
        if os.path.exists(self.key_metadata_file):
            try:
                with open(self.key_metadata_file, 'r') as f:
                    metadata = json.load(f)
                
                # Load all active keys
                for key_info in metadata.get('keys', []):
                    if key_info.get('active', False) and not key_info.get('is_primary', False):
                        key_path = key_info.get('path')
                        key_id = key_info.get('id')
                        if key_path and key_id and os.path.exists(key_path):
                            with open(key_path, 'rb') as f:
                                keys[key_id] = f.read()
            except Exception as e:
                logger.error(f"Error loading key metadata: {str(e)}")
        
        logger.info(f"Loaded {len(keys)} encryption keys")
        return keys
    
    def _generate_key_id(self) -> str:
        """Generate a unique key ID."""
        timestamp = datetime.utcnow().strftime('%Y%m%d%H%M%S')
        random_suffix = os.urandom(4).hex()
        return f"key_{timestamp}_{random_suffix}"
    
    def _get_key_id_from_metadata(self, key_filename: str) -> str:
        """Get the key ID from metadata or generate a new one."""
        if os.path.exists(self.key_metadata_file):
            try:
                with open(self.key_metadata_file, 'r') as f:
                    metadata = json.load(f)
                
                for key_info in metadata.get('keys', []):
                    if os.path.basename(key_info.get('path', '')) == key_filename:
                        return key_info.get('id')
            except Exception:
                pass
        
        # If not found, generate a new ID
        return self._generate_key_id()
    
    def _update_key_metadata(self, key_id: str, key_path: str, is_primary: bool = False) -> None:
        """Update the key metadata file."""
        metadata = {'keys': []}
        
        # Load existing metadata if available
        if os.path.exists(self.key_metadata_file):
            try:
                with open(self.key_metadata_file, 'r') as f:
                    metadata = json.load(f)
            except Exception:
                pass
        
        # Update primary key status if needed
        if is_primary:
            for key in metadata.get('keys', []):
                if key.get('is_primary', False):
                    key['is_primary'] = False
        
        # Add or update key info
        key_exists = False
        for key in metadata.get('keys', []):
            if key.get('id') == key_id:
                key['path'] = key_path
                key['is_primary'] = is_primary
                key['active'] = True
                key['created_at'] = key.get('created_at', datetime.utcnow().isoformat())
                key['updated_at'] = datetime.utcnow().isoformat()
                key_exists = True
                break
        
        if not key_exists:
            metadata.setdefault('keys', []).append({
                'id': key_id,
                'path': key_path,
                'is_primary': is_primary,
                'active': True,
                'created_at': datetime.utcnow().isoformat(),
                'updated_at': datetime.utcnow().isoformat()
            })
        
        # Save metadata
        with open(self.key_metadata_file, 'w') as f:
            json.dump(metadata, f, indent=2)
    
    def encrypt(self, data: Union[str, bytes], key_id: Optional[str] = None) -> bytes:
        """
        Encrypt data using the specified key or the primary key.
        
        Args:
            data: The data to encrypt (string or bytes)
            key_id: Optional key ID to use for encryption (defaults to primary key)
            
        Returns:
            Encrypted data as bytes with key ID prefix
        """
        if data is None:
            return None
        
        # Convert string to bytes if needed
        if isinstance(data, str):
            data = data.encode('utf-8')
        
        # Use specified key or primary key
        if key_id and key_id in self.all_keys:
            fernet = Fernet(self.all_keys[key_id])
            using_key_id = key_id
        else:
            fernet = self.fernet
            using_key_id = self.key_id
        
        # Encrypt the data
        encrypted_data = fernet.encrypt(data)
        
        # Prefix with key ID for later decryption
        prefixed_data = f"{using_key_id}:".encode('utf-8') + encrypted_data
        
        return prefixed_data
    
    def decrypt(self, encrypted_data: Union[str, bytes]) -> bytes:
        """
        Decrypt data that was encrypted with any of the available keys.
        
        Args:
            encrypted_data: The encrypted data to decrypt
            
        Returns:
            Decrypted data as bytes
            
        Raises:
            ValueError: If the data format is invalid
            InvalidToken: If decryption fails
        """
        if encrypted_data is None:
            return None
        
        # Convert string to bytes if needed
        if isinstance(encrypted_data, str):
            encrypted_data = encrypted_data.encode('utf-8')
        
        # Extract key ID and actual encrypted data
        try:
            parts = encrypted_data.split(b':', 1)
            if len(parts) != 2:
                raise ValueError("Invalid encrypted data format: missing key ID prefix")
            
            key_id = parts[0].decode('utf-8')
            actual_encrypted_data = parts[1]
        except Exception as e:
            logger.error(f"Error parsing encrypted data: {str(e)}")
            raise ValueError(f"Invalid encrypted data format: {str(e)}")
        
        # Find the appropriate key
        if key_id in self.all_keys:
            fernet = Fernet(self.all_keys[key_id])
        else:
            raise ValueError(f"Encryption key not found: {key_id}")
        
        # Decrypt the data
        try:
            decrypted_data = fernet.decrypt(actual_encrypted_data)
            return decrypted_data
        except InvalidToken as e:
            logger.error(f"Invalid token or corrupted data: {str(e)}")
            raise
        except Exception as e:
            logger.error(f"Error decrypting data: {str(e)}")
            raise
    
    def encrypt_dict_values(self, data: Dict[str, Any], fields_to_encrypt: List[str]) -> Dict[str, Any]:
        """
        Encrypt specific fields in a dictionary.
        
        Args:
            data: Dictionary containing data to encrypt
            fields_to_encrypt: List of field names to encrypt
            
        Returns:
            Dictionary with encrypted values
        """
        if not data or not fields_to_encrypt:
            return data
        
        result = data.copy()
        
        for field in fields_to_encrypt:
            if field in result and result[field] is not None:
                # Convert to string if not already
                if not isinstance(result[field], (str, bytes)):
                    value_str = json.dumps(result[field])
                else:
                    value_str = result[field]
                
                # Encrypt and encode as base64 for storage
                encrypted_value = self.encrypt(value_str)
                result[field] = base64.b64encode(encrypted_value).decode('utf-8')
        
        return result
    
    def decrypt_dict_values(self, data: Dict[str, Any], fields_to_decrypt: List[str]) -> Dict[str, Any]:
        """
        Decrypt specific fields in a dictionary.
        
        Args:
            data: Dictionary containing encrypted data
            fields_to_decrypt: List of field names to decrypt
            
        Returns:
            Dictionary with decrypted values
        """
        if not data or not fields_to_decrypt:
            return data
        
        result = data.copy()
        
        for field in fields_to_decrypt:
            if field in result and result[field] is not None:
                try:
                    # Decode from base64 and decrypt
                    encrypted_value = base64.b64decode(result[field])
                    decrypted_value = self.decrypt(encrypted_value)
                    
                    # Try to parse as JSON if it looks like JSON
                    try:
                        if decrypted_value.startswith(b'{') or decrypted_value.startswith(b'['):
                            result[field] = json.loads(decrypted_value)
                        else:
                            result[field] = decrypted_value.decode('utf-8')
                    except (json.JSONDecodeError, UnicodeDecodeError):
                        # If not valid JSON or UTF-8, return as bytes
                        result[field] = decrypted_value
                except Exception as e:
                    logger.error(f"Error decrypting field {field}: {str(e)}")
                    # Keep the original value on error
        
        return result
    
    def rotate_key(self) -> Dict[str, Any]:
        """
        Generate a new primary key and return information needed for re-encryption.
        
        Returns:
            Dictionary with old and new key information
        """
        # Store old key info
        old_key_id = self.key_id
        old_key = self.primary_key
        old_fernet = self.fernet
        
        # Generate new key
        new_key = Fernet.generate_key()
        new_key_id = self._generate_key_id()
        new_key_file = os.path.join(self.key_dir, f"{new_key_id}.key")
        
        # Save new key
        with open(new_key_file, 'wb') as f:
            f.write(new_key)
        
        # Update primary key file to point to new key
        with open(self.primary_key_file, 'wb') as f:
            f.write(new_key)
        
        # Update metadata
        self._update_key_metadata(new_key_id, new_key_file, is_primary=True)
        
        # Update instance variables
        self.primary_key = new_key
        self.key_id = new_key_id
        self.fernet = Fernet(new_key)
        self.all_keys[new_key_id] = new_key
        
        logger.info(f"Key rotation complete: {old_key_id} -> {new_key_id}")
        
        return {
            'old_key_id': old_key_id,
            'new_key_id': new_key_id,
            'old_key': old_key,
            'new_key': new_key
        }
    
    def encrypt_file(self, file_path: str) -> str:
        """
        Encrypt a file.
        
        Args:
            file_path: Path to the file to encrypt
            
        Returns:
            Path to the encrypted file
        """
        try:
            with open(file_path, 'rb') as f:
                data = f.read()
            
            encrypted_data = self.encrypt(data)
            encrypted_file = f"{file_path}.enc"
            
            with open(encrypted_file, 'wb') as f:
                f.write(encrypted_data)
            
            logger.info(f"Encrypted file: {file_path} -> {encrypted_file}")
            return encrypted_file
        except Exception as e:
            logger.error(f"Error encrypting file {file_path}: {str(e)}")
            raise
    
    def decrypt_file(self, encrypted_file_path: str, output_path: Optional[str] = None) -> str:
        """
        Decrypt an encrypted file.
        
        Args:
            encrypted_file_path: Path to the encrypted file
            output_path: Optional path for the decrypted file
            
        Returns:
            Path to the decrypted file
        """
        try:
            with open(encrypted_file_path, 'rb') as f:
                encrypted_data = f.read()
            
            decrypted_data = self.decrypt(encrypted_data)
            
            if not output_path:
                if encrypted_file_path.endswith('.enc'):
                    output_path = encrypted_file_path[:-4]
                else:
                    output_path = f"{encrypted_file_path}.dec"
            
            with open(output_path, 'wb') as f:
                f.write(decrypted_data)
            
            logger.info(f"Decrypted file: {encrypted_file_path} -> {output_path}")
            return output_path
        except Exception as e:
            logger.error(f"Error decrypting file {encrypted_file_path}: {str(e)}")
            raise
    
    def derive_key_from_password(self, password: str, salt: Optional[bytes] = None) -> Tuple[bytes, bytes]:
        """
        Derive an encryption key from a password using PBKDF2.
        
        Args:
            password: The password to derive the key from
            salt: Optional salt (generated if not provided)
            
        Returns:
            Tuple of (key, salt)
        """
        if salt is None:
            salt = os.urandom(16)
        
        kdf = PBKDF2HMAC(
            algorithm=hashes.SHA256(),
            length=32,
            salt=salt,
            iterations=100000,
            backend=default_backend()
        )
        
        key = base64.urlsafe_b64encode(kdf.derive(password.encode()))
        return key, salt
    
    def hash_sensitive_data(self, data: str) -> str:
        """
        Create a secure hash of sensitive data for comparison purposes.
        
        Args:
            data: The data to hash
            
        Returns:
            Secure hash of the data
        """
        if not data:
            return None
        
        # Use SHA-256 for hashing
        hash_obj = hashlib.sha256(data.encode())
        return hash_obj.hexdigest()


# Singleton instance for app-wide use
_encryption_service_instance = None

def get_encryption_service() -> EncryptionService:
    """Get or create the singleton encryption service instance."""
    global _encryption_service_instance
    if _encryption_service_instance is None:
        _encryption_service_instance = EncryptionService()
    return _encryption_service_instance


# Field encryption decorators for models
def encrypt_fields(fields_to_encrypt: List[str]):
    """
    Decorator for model create/update methods to encrypt specified fields.
    
    Args:
        fields_to_encrypt: List of field names to encrypt
    """
    def decorator(func):
        def wrapper(*args, **kwargs):
            # Get the data argument (usually first or second argument)
            data = None
            if len(args) > 1 and isinstance(args[1], dict):
                data = args[1]
            elif 'data' in kwargs:
                data = kwargs['data']
            
            if data:
                # Get encryption service
                encryption_service = get_encryption_service()
                
                # Encrypt fields
                encrypted_data = encryption_service.encrypt_dict_values(data, fields_to_encrypt)
                
                # Replace data in args or kwargs
                if len(args) > 1 and isinstance(args[1], dict):
                    args_list = list(args)
                    args_list[1] = encrypted_data
                    args = tuple(args_list)
                elif 'data' in kwargs:
                    kwargs['data'] = encrypted_data
            
            # Call the original function
            return func(*args, **kwargs)
        return wrapper
    return decorator


def decrypt_fields(fields_to_decrypt: List[str]):
    """
    Decorator for model get methods to decrypt specified fields.
    
    Args:
        fields_to_decrypt: List of field names to decrypt
    """
    def decorator(func):
        def wrapper(*args, **kwargs):
            # Call the original function
            result = func(*args, **kwargs)
            
            if result:
                # Get encryption service
                encryption_service = get_encryption_service()
                
                # Handle different return types
                if isinstance(result, dict):
                    return encryption_service.decrypt_dict_values(result, fields_to_decrypt)
                elif isinstance(result, list):
                    return [encryption_service.decrypt_dict_values(item, fields_to_decrypt) 
                            if isinstance(item, dict) else item 
                            for item in result]
            
            return result
        return wrapper
    return decorator