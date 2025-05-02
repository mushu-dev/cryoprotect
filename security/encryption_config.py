"""
Encryption Configuration for CryoProtect v2

This module defines which fields in each model should be encrypted.
It provides a centralized configuration for field-level encryption.
"""

# Dictionary mapping model names to lists of fields that should be encrypted
ENCRYPTED_FIELDS = {
    # User-related models
    'user_profile': [
        'email',
        'name'
    ],
    
    # Experiment-related models
    'experiments': [
        'experimental_conditions',
        'text_value'  # May contain sensitive experimental data
    ],
    
    # Toxicity-related models
    'toxicity_data': [
        'value',
        'notes'
    ],
    
    # Lab verification models
    'lab_verifications': [
        'equipment_used',
        'comments'
    ],
    
    # Project-related models
    'projects': [
        'description'  # May contain sensitive project information
    ],
    
    # Shared resources
    'shared_resources': [
        'access_key'  # If implemented for secure sharing
    ]
}

# Dictionary mapping model names to fields that should be hashed for comparison
# rather than stored in encrypted form (e.g., passwords, API keys)
HASHED_FIELDS = {
    'user_profile': [
        'password_reset_token'
    ],
    'api_keys': [
        'key_value'
    ]
}

# Dictionary defining which models/tables need encryption for backups
SENSITIVE_TABLES = [
    'user_profile',
    'experiments',
    'toxicity_data',
    'lab_verifications',
    'projects',
    'shared_resources'
]

def get_encrypted_fields(model_name):
    """
    Get the list of fields that should be encrypted for a given model.
    
    Args:
        model_name: The name of the model/table
        
    Returns:
        List of field names to encrypt
    """
    return ENCRYPTED_FIELDS.get(model_name, [])

def get_hashed_fields(model_name):
    """
    Get the list of fields that should be hashed for a given model.
    
    Args:
        model_name: The name of the model/table
        
    Returns:
        List of field names to hash
    """
    return HASHED_FIELDS.get(model_name, [])

def is_sensitive_table(table_name):
    """
    Check if a table contains sensitive data that should be encrypted in backups.
    
    Args:
        table_name: The name of the table
        
    Returns:
        True if the table contains sensitive data, False otherwise
    """
    return table_name in SENSITIVE_TABLES