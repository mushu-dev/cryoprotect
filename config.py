"""
CryoProtect v2 - Hierarchical Environment Configuration System

This module implements a robust, hierarchical configuration system for CryoProtect v2, supporting:
- A BaseConfig with environment-specific subclasses
- Type validation, required/optional enforcement, and runtime error handling
- Fallback and override precedence (env vars, .env, Docker/CI/CD secrets, hardcoded defaults)
- Per-environment overrides and extensibility
- Compatibility with Docker, CI/CD, and security best practices
"""

import os
import json
import sys
import typing
from typing import Any, Dict, List, Optional, Type, TypeVar, Union, get_type_hints
from pathlib import Path
from urllib.parse import urlparse
from dotenv import load_dotenv

# Load environment variables from .env file
load_dotenv()

# Type variable for config class types
ConfigType = TypeVar('ConfigType', bound='BaseConfig')


class ConfigurationError(Exception):
    """Exception raised for configuration errors."""
    pass


class BaseConfig:
    """
    Base configuration class that defines common configuration variables,
    their types, required/optional status, and default values.
    
    All configuration variables should be defined as class variables with
    type annotations. Required variables should not have default values.
    """
    
    # Flask configuration
    SECRET_KEY: str = os.getenv('SECRET_KEY', 'dev-key-please-change-in-production')
    DEBUG: bool = False
    TESTING: bool = False
    
    # API settings
    API_TITLE: str = 'CryoProtect Analyzer API'
    API_VERSION: str = 'v1'
    OPENAPI_VERSION: str = '3.0.2'
    OPENAPI_URL_PREFIX: str = '/'
    OPENAPI_SWAGGER_UI_PATH: str = '/swagger'
    OPENAPI_SWAGGER_UI_URL: str = 'https://cdn.jsdelivr.net/npm/swagger-ui-dist/'
    
    # Supabase connection settings
    SUPABASE_URL: str  # Required, no default
    SUPABASE_KEY: str  # Required, no default
    SUPABASE_SERVICE_KEY: Optional[str] = None
    
    # Authentication
    SUPABASE_USER: Optional[str] = None
    SUPABASE_PASSWORD: Optional[str] = None
    
    # Connection pool settings
    SUPABASE_MIN_CONNECTIONS: int = 2
    SUPABASE_MAX_CONNECTIONS: int = 10
    SUPABASE_CONNECTION_TIMEOUT: int = 30
    SUPABASE_CONNECTION_LIFETIME: int = 3600
    SUPABASE_IDLE_TIMEOUT: int = 300
    
    # ChEMBL API configuration
    CHEMBL_API_URL: Optional[str] = "https://www.ebi.ac.uk/chembl/api/data"
    CHEMBL_API_KEY: Optional[str] = None  # API key for ChEMBL authentication
    CHEMBL_API_KEY_REQUIRED: bool = False  # Set to True if API key is required
    CHEMBL_API_DELAY: float = 0.3  # Seconds between API requests
    CHEMBL_ID_FILE: Optional[str] = None  # Path to file containing ChEMBL IDs
    
    # ChEMBL client configuration
    CHEMBL_CACHE_DIR: str = "cache/chembl"
    CHEMBL_REQUESTS_PER_SECOND: float = 5.0  # Maximum requests per second (spec: ≤5 requests/sec)
    CHEMBL_MAX_RETRIES: int = 5
    CHEMBL_FAILURE_THRESHOLD: int = 3
    CHEMBL_RECOVERY_TIMEOUT: int = 60
    CHEMBL_CACHE_TTL: int = 86400 * 30  # 30 days
    CHEMBL_MEMORY_CACHE_SIZE: int = 1000
    CHEMBL_MEMORY_THRESHOLD: float = 80.0  # Memory usage threshold (percentage) to trigger rate limiting
    
    # Batch processing and checkpointing
    CHEMBL_BATCH_SIZE: int = 100  # Default batch size (spec: 100-500 molecules per batch)
    CHEMBL_MEMORY_CHECK_FREQUENCY: int = 10  # Check memory usage every N molecules within a batch
    CHECKPOINT_DIR: str = "checkpoints"
    
    # Logging configuration
    LOG_LEVEL: str = "INFO"
    LOG_FILE: str = "cryoprotectant_analysis.log"
    LOG_TO_FILE: bool = True
    
    # Security configuration
    CORS_ORIGINS: List[str] = ["http://localhost:5000"]
    
    # Feature flags
    ENABLE_EXPERIMENTAL_FEATURES: bool = False
    
    # Caching configuration
    CACHE_TYPE: str = "memory"
    CACHE_REDIS_URL: Optional[str] = None
    CACHE_DEFAULT_TIMEOUT: int = 300
    
    # Rate limiting configuration
    RATE_LIMIT_ENABLED: bool = True
    RATE_LIMIT_STORAGE_URL: str = "memory://"
    RATE_LIMIT_STRATEGY: str = "fixed-window"
    RATE_LIMIT_BY: str = "hybrid"
    RATE_LIMIT_HEADERS_ENABLED: bool = True
    RATE_LIMIT_RETRY_AFTER: int = 60
    RATE_LIMIT_ROLES: Dict[str, List[str]] = {
        "admin": ["5000 per day", "500 per hour", "100 per minute"],
        "premium": ["2000 per day", "200 per hour", "40 per minute"],
        "basic": ["1000 per day", "100 per hour", "20 per minute"]
    }
    
    def __init__(self):
        """Initialize the configuration with environment variables."""
        self._load_from_env()
        self.validate()
    
    def _load_from_env(self):
        """Load configuration from environment variables."""
        type_hints = get_type_hints(self.__class__)
        
        for name, type_hint in type_hints.items():
            # Skip private attributes
            if name.startswith('_'):
                continue
            
            # Check if this is an environment-specific variable
            env_prefix = self._get_env_prefix()
            env_var_name = f"{env_prefix}{name}" if env_prefix else name
            
            # Get the value from environment variables with fallback logic
            value = self._get_env_var(
                name=name,
                env_var_name=env_var_name,
                type_hint=type_hint
            )
            
            if value is not None:
                setattr(self, name, value)
    
    def _get_env_prefix(self) -> str:
        """
        Get the environment prefix for this configuration class.
        Override in subclasses to provide environment-specific prefixes.
        """
        return ""
    
    def _get_env_var(self, name: str, env_var_name: str, type_hint: Type) -> Any:
        """
        Get an environment variable with the given name and convert it to the specified type.
        
        Args:
            name: The name of the configuration variable
            env_var_name: The name of the environment variable
            type_hint: The expected type of the variable
            
        Returns:
            The value of the environment variable converted to the specified type,
            or the default value if the environment variable is not set.
        """
        # Check if the variable has a default value
        default_value = getattr(self.__class__, name, None)
        required = not hasattr(self.__class__, name)
        
        # Check for Docker secrets
        docker_secret_path = f"/run/secrets/{env_var_name}"
        if os.path.exists(docker_secret_path):
            with open(docker_secret_path, 'r') as f:
                value = f.read().strip()
            return self._convert_value(name, value, type_hint)
        
        # Check environment variables
        value = os.getenv(env_var_name)
        if value is not None:
            return self._convert_value(name, value, type_hint)
        
        # Check environment variables without prefix
        if env_var_name != name:
            value = os.getenv(name)
            if value is not None:
                return self._convert_value(name, value, type_hint)
        
        # Return default value or None
        return default_value
    
    def _convert_value(self, name: str, value: str, type_hint: Type) -> Any:
        """
        Convert a string value to the specified type.
        
        Args:
            name: The name of the configuration variable
            value: The string value to convert
            type_hint: The expected type of the variable
            
        Returns:
            The converted value
            
        Raises:
            ConfigurationError: If the value cannot be converted to the specified type
        """
        # Handle Optional types
        if hasattr(type_hint, "__origin__") and type_hint.__origin__ is Union:
            if type(None) in type_hint.__args__:
                # This is an Optional type
                if value == "":
                    return None
                # Get the actual type (excluding None)
                actual_type = next(arg for arg in type_hint.__args__ if arg is not type(None))
                return self._convert_value(name, value, actual_type)
        
        try:
            # Handle basic types
            if type_hint is str:
                return value
            elif type_hint is bool:
                return value.lower() in ('true', 'yes', '1', 'y', 'on')
            elif type_hint is int:
                return int(value)
            elif type_hint is float:
                return float(value)
            # Handle container types
            elif hasattr(type_hint, "__origin__"):
                if type_hint.__origin__ is list or type_hint.__origin__ is List:
                    return json.loads(value)
                elif type_hint.__origin__ is dict or type_hint.__origin__ is Dict:
                    return json.loads(value)
            
            # Default: try to use the type as a constructor
            return type_hint(value)
        except (ValueError, TypeError, json.JSONDecodeError) as e:
            source = self._get_value_source(name) if hasattr(self, '_get_value_source') else "Unknown source"
            raise ConfigurationError(
                f"Invalid value for {name}: '{value}'. Expected type: {type_hint}. Source: {source}. Error: {str(e)}\n\n"
                f"Please check your environment variables, .env file, or Docker secrets."
            )
    
    def validate(self):
        """
        Validate the configuration.
        
        Performs comprehensive validation of all configuration variables:
        - Checks if required variables are present
        - Validates types for all variables
        - Performs format validation for specific variables (URLs, emails, etc.)
        - Validates environment-specific requirements
        
        Raises:
            ConfigurationError: If any validation fails, with detailed error messages
        """
        errors = []
        type_hints = get_type_hints(self.__class__)
        env_name = self.__class__.__name__.replace('Config', '')
        
        for name, type_hint in type_hints.items():
            # Skip private attributes
            if name.startswith('_'):
                continue
            
            # Check if the variable is required
            required = not hasattr(self.__class__, name)
            value = getattr(self, name, None)
            source = self._get_value_source(name)
            
            # Check if required variables are present
            if required and value is None:
                errors.append(f"Required configuration variable '{name}' is missing in {env_name} environment")
                continue
            
            # Skip type checking for None values in optional fields
            if value is None and not required:
                continue
            
            # Check types
            if not self._check_type(value, type_hint):
                errors.append(
                    f"Invalid type for '{name}' in {env_name} environment. Expected {type_hint}, got {type(value)}. "
                    f"Source: {source}"
                )
                continue
            
            # Format validation for specific variable types
            self._validate_format(name, value, errors)
            
            # Environment-specific validation
            self._validate_environment_specific(name, value, errors)
        
        # Separate critical errors from warnings
        critical_errors = []
        warnings = []
        
        for error in errors:
            if error.startswith("Warning:"):
                warnings.append(error)
            else:
                critical_errors.append(error)
        
        # Log warnings but don't fail validation for them in development/testing/base config
        # Also include custom config classes that inherit directly from BaseConfig
        is_dev_or_test = (isinstance(self, DevelopmentConfig) or
                         isinstance(self, TestingConfig) or
                         self.__class__ is BaseConfig or
                         self.__class__.__bases__[0] is BaseConfig)
        
        if warnings and is_dev_or_test:
            print("\nConfiguration Warnings:")
            for warning in warnings:
                print(f"- {warning}")
            print()
        
        # Only raise an error for critical errors
        if critical_errors:
            error_message = f"Configuration validation failed for {env_name} environment:\n" + "\n".join(f"- {error}" for error in critical_errors)
            if warnings:
                error_message += "\n\nWarnings:\n" + "\n".join(f"- {warning}" for warning in warnings)
            error_message += "\n\nPlease check your environment variables, .env file, or Docker secrets."
            raise ConfigurationError(error_message)
        elif warnings and not is_dev_or_test:
            # For non-development environments, warnings are treated as errors
            error_message = f"Configuration validation failed for {env_name} environment (warnings treated as errors):\n" + "\n".join(f"- {warning}" for warning in warnings)
            error_message += "\n\nPlease check your environment variables, .env file, or Docker secrets."
            raise ConfigurationError(error_message)
    
    def _check_type(self, value: Any, expected_type: Type) -> bool:
        """
        Check if a value matches the expected type.
        
        Args:
            value: The value to check
            expected_type: The expected type
            
        Returns:
            True if the value matches the expected type, False otherwise
        """
        # Handle Optional types
        if hasattr(expected_type, "__origin__") and expected_type.__origin__ is Union:
            return any(self._check_type(value, arg) for arg in expected_type.__args__)
        
        # Handle container types
        if hasattr(expected_type, "__origin__"):
            if expected_type.__origin__ is list or expected_type.__origin__ is List:
                return isinstance(value, list)
            elif expected_type.__origin__ is dict or expected_type.__origin__ is Dict:
                return isinstance(value, dict)
        
        # Handle basic types
        return isinstance(value, expected_type)
    
    def _validate_format(self, name: str, value: Any, errors: List[str]):
        """
        Validate the format of specific configuration variables.
        
        Args:
            name: The name of the configuration variable
            value: The value to validate
            errors: List to append error messages to
        """
        if value is None:
            return
            
        # URL validation
        if name.endswith('_URL') and isinstance(value, str):
            # Special case for memory:// URLs used in caching and rate limiting
            if value in ['memory://', 'simple://']:
                return
                
            try:
                result = urlparse(value)
                if not all([result.scheme, result.netloc]):
                    errors.append(f"Invalid URL format for {name}: '{value}'. URLs must include scheme and host.")
            except Exception as e:
                errors.append(f"Invalid URL format for {name}: '{value}'. Error: {str(e)}")
        
        # Email validation
        elif name.endswith('_EMAIL') and isinstance(value, str):
            if '@' not in value or '.' not in value.split('@')[1]:
                errors.append(f"Invalid email format for {name}: '{value}'. Must be a valid email address.")
        
        # Path validation
        elif name.endswith('_PATH') and isinstance(value, str):
            # Skip validation for API paths and Swagger UI paths
            if name in ['OPENAPI_SWAGGER_UI_PATH', 'OPENAPI_URL_PREFIX']:
                return
                
            if not os.path.exists(value) and not Path(value).is_absolute():
                errors.append(f"Warning: Path for {name} does not exist: '{value}'")
        
        # JSON validation
        elif name.endswith('_JSON') and isinstance(value, str):
            try:
                json.loads(value)
            except json.JSONDecodeError as e:
                errors.append(f"Invalid JSON format for {name}: '{value}'. Error: {str(e)}")
        
        # Security validation
        if name == 'SECRET_KEY' and value == 'dev-key-please-change-in-production':
            if isinstance(self, ProductionConfig):
                errors.append(f"Insecure default SECRET_KEY detected in production environment. Please set a strong, unique SECRET_KEY.")
            else:
                errors.append(f"Warning: Using default SECRET_KEY in {self.__class__.__name__}. This is fine for development but must be changed in production.")
    
    def _validate_environment_specific(self, name: str, value: Any, errors: List[str]):
        """
        Perform environment-specific validation.
        
        Args:
            name: The name of the configuration variable
            value: The value to validate
            errors: List to append error messages to
        """
        # Production-specific validations
        if isinstance(self, ProductionConfig):
            if name == 'DEBUG' and value:
                errors.append("DEBUG should be False in production environment")
            
            if name == 'TESTING' and value:
                errors.append("TESTING should be False in production environment")
            
            if name == 'ENABLE_EXPERIMENTAL_FEATURES' and value:
                errors.append("Warning: ENABLE_EXPERIMENTAL_FEATURES is True in production environment")
    
    def _get_value_source(self, name: str) -> str:
        """
        Determine the source of a configuration value.
        
        Args:
            name: The name of the configuration variable
            
        Returns:
            A string indicating the source of the value (env var, Docker secret, default)
        """
        env_prefix = self._get_env_prefix()
        env_var_name = f"{env_prefix}{name}" if env_prefix else name
        
        # Check Docker secrets
        docker_secret_path = f"/run/secrets/{env_var_name}"
        if os.path.exists(docker_secret_path):
            return f"Docker secret ({docker_secret_path})"
        
        # Check environment variables with prefix
        if os.getenv(env_var_name) is not None:
            return f"Environment variable ({env_var_name})"
        
        # Check environment variables without prefix
        if os.getenv(name) is not None:
            return f"Environment variable ({name})"
        
        # Must be a default value
        return "Default value"
    
    @classmethod
    def from_env(cls: Type[ConfigType]) -> ConfigType:
        """
        Create a configuration instance based on the current environment.
        
        Returns:
            An instance of the appropriate configuration class for the current environment
        """
        # Check for FLASK_ENV or APP_ENV environment variables
        flask_env = os.getenv('FLASK_ENV', os.getenv('APP_ENV', 'development')).lower()
        config_class = config_classes.get(flask_env, DevelopmentConfig)
        return config_class()
    
    def as_dict(self) -> Dict[str, Any]:
        """
        Convert the configuration to a dictionary.
        
        Returns:
            A dictionary containing all configuration variables
        """
        result = {}
        for key in dir(self):
            # Skip private attributes and methods
            if not key.startswith('_') and not callable(getattr(self, key)):
                result[key] = getattr(self, key)
        return result
    
    def print_config(self, include_secrets: bool = False):
        """
        Print the configuration to the console.
        
        Args:
            include_secrets: Whether to include secret values in the output
        """
        print(f"\n{'=' * 50}")
        print(f"Configuration: {self.__class__.__name__}")
        print(f"{'=' * 50}")
        
        for key, value in sorted(self.as_dict().items()):
            # Skip methods and private attributes
            if callable(value) or key.startswith('_'):
                continue
            
            # Mask secret values
            if not include_secrets and any(secret in key.lower() for secret in ['key', 'password', 'secret', 'token']):
                if value:
                    value = f"{value[:4]}...{value[-4:]}" if len(str(value)) > 8 else "****"
            
            print(f"{key}: {value}")
        print(f"{'=' * 50}\n")


class DevelopmentConfig(BaseConfig):
    """Development environment configuration."""
    
    DEBUG: bool = True
    TESTING: bool = False
    LOG_LEVEL: str = "DEBUG"
    ENABLE_EXPERIMENTAL_FEATURES: bool = True
    
    def _get_env_prefix(self) -> str:
        return "DEVELOPMENT_"


class TestingConfig(BaseConfig):
    """Testing environment configuration."""
    
    DEBUG: bool = True
    TESTING: bool = True
    LOG_LEVEL: str = "DEBUG"
    
    def _get_env_prefix(self) -> str:
        return "TESTING_"


class StagingConfig(BaseConfig):
    """Staging environment configuration."""
    
    DEBUG: bool = False
    TESTING: bool = False
    LOG_LEVEL: str = "INFO"
    CORS_ORIGINS: List[str] = ["https://staging.cryoprotect-analyzer.com"]
    ENABLE_EXPERIMENTAL_FEATURES: bool = True
    
    def _get_env_prefix(self) -> str:
        return "STAGING_"


class ProductionConfig(BaseConfig):
    """Production environment configuration."""
    
    DEBUG: bool = False
    TESTING: bool = False
    LOG_LEVEL: str = "WARNING"
    CORS_ORIGINS: List[str] = ["https://cryoprotect-analyzer.com"]
    ENABLE_EXPERIMENTAL_FEATURES: bool = False
    CACHE_TYPE: str = "redis"
    
    def _get_env_prefix(self) -> str:
        return "PRODUCTION_"


# Dictionary with different configuration environments
config_classes = {
    'development': DevelopmentConfig,
    'testing': TestingConfig,
    'staging': StagingConfig,
    'production': ProductionConfig,
}

# Create the active configuration instance
active_config = BaseConfig.from_env()

# For module-level imports in other files, provide these variables
# These will be set based on the active configuration
SUPABASE_URL = active_config.SUPABASE_URL
SUPABASE_KEY = active_config.SUPABASE_KEY
SUPABASE_SERVICE_KEY = active_config.SUPABASE_SERVICE_KEY
SUPABASE_USER = active_config.SUPABASE_USER
SUPABASE_PASSWORD = active_config.SUPABASE_PASSWORD
SUPABASE_MIN_CONNECTIONS = active_config.SUPABASE_MIN_CONNECTIONS
SUPABASE_MAX_CONNECTIONS = active_config.SUPABASE_MAX_CONNECTIONS
SUPABASE_CONNECTION_TIMEOUT = active_config.SUPABASE_CONNECTION_TIMEOUT
SUPABASE_CONNECTION_LIFETIME = active_config.SUPABASE_CONNECTION_LIFETIME
SUPABASE_IDLE_TIMEOUT = active_config.SUPABASE_IDLE_TIMEOUT

# Additional module-level variables for convenience
DEBUG = active_config.DEBUG
TESTING = active_config.TESTING
SECRET_KEY = active_config.SECRET_KEY
API_TITLE = active_config.API_TITLE
API_VERSION = active_config.API_VERSION
LOG_LEVEL = active_config.LOG_LEVEL
CORS_ORIGINS = active_config.CORS_ORIGINS
ENABLE_EXPERIMENTAL_FEATURES = active_config.ENABLE_EXPERIMENTAL_FEATURES
CACHE_TYPE = active_config.CACHE_TYPE
CACHE_REDIS_URL = active_config.CACHE_REDIS_URL
CACHE_DEFAULT_TIMEOUT = active_config.CACHE_DEFAULT_TIMEOUT

# ChEMBL API configuration
CHEMBL_API_URL = active_config.CHEMBL_API_URL
CHEMBL_API_KEY = active_config.CHEMBL_API_KEY
CHEMBL_API_KEY_REQUIRED = active_config.CHEMBL_API_KEY_REQUIRED
CHEMBL_API_DELAY = active_config.CHEMBL_API_DELAY
CHEMBL_CACHE_DIR = active_config.CHEMBL_CACHE_DIR
CHEMBL_REQUESTS_PER_SECOND = active_config.CHEMBL_REQUESTS_PER_SECOND
CHEMBL_MAX_RETRIES = active_config.CHEMBL_MAX_RETRIES
CHEMBL_FAILURE_THRESHOLD = active_config.CHEMBL_FAILURE_THRESHOLD
CHEMBL_RECOVERY_TIMEOUT = active_config.CHEMBL_RECOVERY_TIMEOUT
CHEMBL_CACHE_TTL = active_config.CHEMBL_CACHE_TTL
CHEMBL_MEMORY_CACHE_SIZE = active_config.CHEMBL_MEMORY_CACHE_SIZE
CHEMBL_MEMORY_THRESHOLD = active_config.CHEMBL_MEMORY_THRESHOLD
CHEMBL_BATCH_SIZE = active_config.CHEMBL_BATCH_SIZE
CHEMBL_MEMORY_CHECK_FREQUENCY = active_config.CHEMBL_MEMORY_CHECK_FREQUENCY
CHECKPOINT_DIR = active_config.CHECKPOINT_DIR


def load_environment_variables() -> None:
    """
    Load and normalize environment variables.
    Ensures variables are consistently available regardless of naming convention.
    """
    # DB_* and SUPABASE_DB_* variables normalization
    db_vars = {
        'HOST': os.getenv('DB_HOST') or os.getenv('SUPABASE_DB_HOST'),
        'PORT': os.getenv('DB_PORT') or os.getenv('SUPABASE_DB_PORT') or '5432',
        'NAME': os.getenv('DB_NAME') or os.getenv('SUPABASE_DB_NAME') or 'postgres',
        'USER': os.getenv('DB_USER') or os.getenv('SUPABASE_DB_USER'),
        'PASSWORD': os.getenv('DB_PASSWORD') or os.getenv('SUPABASE_DB_PASSWORD'),
        'IP_ADDRESS': os.getenv('DB_IP_ADDRESS') or os.getenv('SUPABASE_DB_IP_ADDRESS'),
        'MIN_CONNECTIONS': os.getenv('DB_MIN_CONNECTIONS') or os.getenv('SUPABASE_DB_MIN_CONNECTIONS') or '1',
        'MAX_CONNECTIONS': os.getenv('DB_MAX_CONNECTIONS') or os.getenv('SUPABASE_DB_MAX_CONNECTIONS') or '10'
    }
    
    # Set both DB_* and SUPABASE_DB_* variables
    for key, value in db_vars.items():
        if value:
            if not os.getenv(f'DB_{key}'):
                os.environ[f'DB_{key}'] = value
            if not os.getenv(f'SUPABASE_DB_{key}'):
                os.environ[f'SUPABASE_DB_{key}'] = value

def get_db_config() -> Dict[str, Any]:
    """
    Get database configuration based on connection mode.
    
    Returns:
        Dict containing database configuration parameters
    """
    # Ensure environment variables are loaded and normalized
    load_environment_variables()
    
    connection_mode = os.getenv('DB_CONNECTION_MODE', 'auto').lower()
    
    if connection_mode == 'local' or connection_mode == 'auto':
        return {
            'local': {
                'host': os.getenv('LOCAL_DB_HOST', 'localhost'),
                'port': os.getenv('LOCAL_DB_PORT', '5432'),
                'database': os.getenv('LOCAL_DB_NAME', 'cryoprotect'),
                'user': os.getenv('LOCAL_DB_USER', 'postgres'),
                'password': os.getenv('LOCAL_DB_PASSWORD', ''),
                'min_connections': int(os.getenv('LOCAL_DB_MIN_CONNECTIONS', '1')),
                'max_connections': int(os.getenv('LOCAL_DB_MAX_CONNECTIONS', '5'))
            }
        }
    
    if connection_mode == 'supabase' or connection_mode == 'auto':
        # Ensure we have the IP address resolved if available
        ip_address = os.getenv('SUPABASE_DB_IP_ADDRESS')
        
        return {
            'supabase': {
                'host': os.getenv('SUPABASE_DB_HOST'),
                'port': os.getenv('SUPABASE_DB_PORT', '5432'),
                'database': os.getenv('SUPABASE_DB_NAME', 'postgres'),
                'user': os.getenv('SUPABASE_DB_USER'),
                'password': os.getenv('SUPABASE_DB_PASSWORD'),
                'ip_address': ip_address,
                'min_connections': int(os.getenv('SUPABASE_DB_MIN_CONNECTIONS', '1')),
                'max_connections': int(os.getenv('SUPABASE_DB_MAX_CONNECTIONS', '10'))
            }
        }
    
    if connection_mode == 'mcp' or connection_mode == 'auto':
        return {
            'mcp': {
                'project_id': os.getenv('SUPABASE_PROJECT_ID')
            }
        }
    
    # Default to empty config if no valid connection mode
    return {}


def validate_config():
    """
    Validate the configuration and exit if invalid.
    This function should be called at application startup.
    
    Returns:
        True if the configuration is valid
        
    Exits with status code 1 if the configuration is invalid
    """
    try:
        # Perform standard configuration validation
        active_config.validate()
        
        # Additional validation for ChEMBL API credentials
        # Check if ChEMBL API key is required but missing
        # This is application-specific logic that goes beyond the basic type/presence validation
        if getattr(active_config, 'CHEMBL_API_KEY_REQUIRED', False) and not active_config.CHEMBL_API_KEY:
            raise ConfigurationError(
                "ChEMBL API Key is required but not provided. "
                "Please set the CHEMBL_API_KEY environment variable."
            )
            
        print(f"Configuration validated successfully for {active_config.__class__.__name__}")
        return True
    except ConfigurationError as e:
        print(f"\n{'!' * 80}", file=sys.stderr)
        print(f"CONFIGURATION ERROR", file=sys.stderr)
        print(f"{'!' * 80}\n", file=sys.stderr)
        print(f"{str(e)}", file=sys.stderr)
        print(f"\nApplication startup aborted due to configuration errors.", file=sys.stderr)
        print(f"Please fix the above issues and restart the application.", file=sys.stderr)
        print(f"\n{'!' * 80}\n", file=sys.stderr)
        sys.exit(1)


# If this module is run directly, print the configuration
if __name__ == "__main__":
    active_config.print_config()