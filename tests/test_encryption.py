#!/usr/bin/env python3
"""
Test Encryption Service

This script tests the encryption service to ensure it works correctly.
It verifies key generation, field-level encryption, file encryption,
and integration with database models.

Usage:
    python tests/test_encryption.py
"""

import os
import sys
import json
import tempfile
import unittest
from datetime import datetime
from pathlib import Path

# Add parent directory to path to allow imports
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from security.encryption import EncryptionService, encrypt_fields, decrypt_fields
from security.encryption_config import get_encrypted_fields, get_hashed_fields

class TestEncryptionService(unittest.TestCase):
    """Test cases for the encryption service."""
    
    def setUp(self):
        """Set up test environment."""
        # Create a temporary directory for test keys
        self.temp_dir = tempfile.TemporaryDirectory()
        self.key_dir = Path(self.temp_dir.name) / "keys"
        self.key_dir.mkdir(exist_ok=True)
        
        # Create a test encryption service
        self.primary_key_file = self.key_dir / "primary.key"
        self.key_metadata_file = self.key_dir / "key_metadata.json"
        self.encryption_service = EncryptionService(
            primary_key_file=str(self.primary_key_file),
            key_metadata_file=str(self.key_metadata_file),
            auto_create_keys=True
        )
    
    def tearDown(self):
        """Clean up after tests."""
        self.temp_dir.cleanup()
    
    def test_key_generation(self):
        """Test key generation and loading."""
        # Verify that the primary key was created
        self.assertTrue(self.primary_key_file.exists())
        self.assertTrue(self.key_metadata_file.exists())
        
        # Verify that the key metadata contains the primary key
        with open(self.key_metadata_file, 'r') as f:
            metadata = json.load(f)
        
        self.assertIn('keys', metadata)
        self.assertTrue(len(metadata['keys']) > 0)
        
        # Verify that at least one key is marked as primary
        primary_keys = [k for k in metadata['keys'] if k.get('is_primary', False)]
        self.assertTrue(len(primary_keys) > 0)
    
    def test_encryption_decryption(self):
        """Test basic encryption and decryption."""
        # Test with string data
        original_text = "This is sensitive data that needs encryption"
        encrypted_data = self.encryption_service.encrypt(original_text)
        decrypted_data = self.encryption_service.decrypt(encrypted_data)
        
        self.assertNotEqual(original_text, encrypted_data)
        self.assertEqual(original_text, decrypted_data.decode('utf-8'))
        
        # Test with dictionary data
        original_dict = {"name": "John Doe", "email": "john@example.com", "age": 30}
        fields_to_encrypt = ["name", "email"]
        
        encrypted_dict = self.encryption_service.encrypt_dict_values(original_dict, fields_to_encrypt)
        
        # Verify that specified fields are encrypted
        self.assertNotEqual(original_dict["name"], encrypted_dict["name"])
        self.assertNotEqual(original_dict["email"], encrypted_dict["email"])
        
        # Verify that non-specified fields are unchanged
        self.assertEqual(original_dict["age"], encrypted_dict["age"])
        
        # Decrypt and verify
        decrypted_dict = self.encryption_service.decrypt_dict_values(encrypted_dict, fields_to_encrypt)
        self.assertEqual(original_dict["name"], decrypted_dict["name"])
        self.assertEqual(original_dict["email"], decrypted_dict["email"])
        self.assertEqual(original_dict["age"], decrypted_dict["age"])
    
    def test_key_rotation(self):
        """Test key rotation."""
        # Encrypt data with the original key
        original_text = "Data encrypted with the original key"
        encrypted_data = self.encryption_service.encrypt(original_text)
        
        # Rotate the key
        rotation_info = self.encryption_service.rotate_key()
        
        # Verify that we can still decrypt data encrypted with the old key
        decrypted_data = self.encryption_service.decrypt(encrypted_data)
        self.assertEqual(original_text, decrypted_data.decode('utf-8'))
        
        # Encrypt new data with the new key
        new_text = "Data encrypted with the new key"
        new_encrypted_data = self.encryption_service.encrypt(new_text)
        
        # Verify that we can decrypt data encrypted with the new key
        new_decrypted_data = self.encryption_service.decrypt(new_encrypted_data)
        self.assertEqual(new_text, new_decrypted_data.decode('utf-8'))
        
        # Verify that the key metadata was updated
        with open(self.key_metadata_file, 'r') as f:
            metadata = json.load(f)
        
        # There should be at least two keys now
        self.assertTrue(len(metadata['keys']) >= 2)
        
        # Only one key should be marked as primary
        primary_keys = [k for k in metadata['keys'] if k.get('is_primary', False)]
        self.assertEqual(len(primary_keys), 1)
        
        # The primary key should be the new key
        self.assertEqual(primary_keys[0]['id'], rotation_info['new_key_id'])
    
    def test_file_encryption(self):
        """Test file encryption and decryption."""
        # Create a test file
        test_file = Path(self.temp_dir.name) / "test_file.txt"
        test_content = "This is a test file with sensitive content"
        with open(test_file, 'w') as f:
            f.write(test_content)
        
        # Encrypt the file
        encrypted_file = self.encryption_service.encrypt_file(str(test_file))
        
        # Verify that the encrypted file exists and is different from the original
        self.assertTrue(os.path.exists(encrypted_file))
        with open(encrypted_file, 'rb') as f:
            encrypted_content = f.read()
        self.assertNotIn(test_content.encode(), encrypted_content)
        
        # Decrypt the file
        decrypted_file = self.encryption_service.decrypt_file(encrypted_file)
        
        # Verify that the decrypted content matches the original
        with open(decrypted_file, 'r') as f:
            decrypted_content = f.read()
        self.assertEqual(test_content, decrypted_content)
    
    def test_password_derived_key(self):
        """Test key derivation from password."""
        password = "secure-password-123"
        
        # Derive a key from the password
        key, salt = self.encryption_service.derive_key_from_password(password)
        
        # Derive the same key again with the same salt
        key2, _ = self.encryption_service.derive_key_from_password(password, salt)
        
        # Verify that the keys match
        self.assertEqual(key, key2)
        
        # Derive a key with a different password
        different_password = "different-password-456"
        key3, _ = self.encryption_service.derive_key_from_password(different_password, salt)
        
        # Verify that the keys are different
        self.assertNotEqual(key, key3)
    
    def test_hash_sensitive_data(self):
        """Test secure hashing of sensitive data."""
        data = "sensitive-data-to-hash"
        
        # Hash the data
        hash1 = self.encryption_service.hash_sensitive_data(data)
        
        # Hash the same data again
        hash2 = self.encryption_service.hash_sensitive_data(data)
        
        # Verify that the hashes match
        self.assertEqual(hash1, hash2)
        
        # Hash different data
        different_data = "different-data-to-hash"
        hash3 = self.encryption_service.hash_sensitive_data(different_data)
        
        # Verify that the hashes are different
        self.assertNotEqual(hash1, hash3)
    
    def test_field_decorators(self):
        """Test the field encryption decorators."""
        # Define a mock model class with encryption decorators
        class MockModel:
            @encrypt_fields(['name', 'email'])
            def create(self, data):
                return data
            
            @decrypt_fields(['name', 'email'])
            def get(self, id):
                # This would normally fetch from the database
                # For testing, we'll just return encrypted data
                return self.encrypted_data
        
        # Create an instance and test
        mock_model = MockModel()
        
        # Test create with encryption
        original_data = {"name": "Jane Doe", "email": "jane@example.com", "role": "admin"}
        encrypted_data = mock_model.create(original_data)
        
        # Verify that specified fields are encrypted
        self.assertNotEqual(original_data["name"], encrypted_data["name"])
        self.assertNotEqual(original_data["email"], encrypted_data["email"])
        
        # Verify that non-specified fields are unchanged
        self.assertEqual(original_data["role"], encrypted_data["role"])
        
        # Store the encrypted data for the get test
        mock_model.encrypted_data = encrypted_data
        
        # Test get with decryption
        decrypted_data = mock_model.get(1)  # ID doesn't matter for this test
        
        # Verify that the data was decrypted correctly
        self.assertEqual(original_data["name"], decrypted_data["name"])
        self.assertEqual(original_data["email"], decrypted_data["email"])
        self.assertEqual(original_data["role"], decrypted_data["role"])


class TestEncryptionConfig(unittest.TestCase):
    """Test cases for the encryption configuration."""
    
    def test_get_encrypted_fields(self):
        """Test getting encrypted fields for models."""
        # Test a model with encrypted fields
        user_profile_fields = get_encrypted_fields('user_profile')
        self.assertTrue(len(user_profile_fields) > 0)
        self.assertIn('email', user_profile_fields)
        
        # Test a model without encrypted fields
        nonexistent_model_fields = get_encrypted_fields('nonexistent_model')
        self.assertEqual(len(nonexistent_model_fields), 0)
    
    def test_get_hashed_fields(self):
        """Test getting hashed fields for models."""
        # Test a model with hashed fields
        user_profile_fields = get_hashed_fields('user_profile')
        self.assertTrue(len(user_profile_fields) > 0)
        self.assertIn('password_reset_token', user_profile_fields)
        
        # Test a model without hashed fields
        nonexistent_model_fields = get_hashed_fields('nonexistent_model')
        self.assertEqual(len(nonexistent_model_fields), 0)


def run_tests():
    """Run the encryption tests."""
    unittest.main()

if __name__ == "__main__":
    run_tests()