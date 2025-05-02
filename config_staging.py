"""
CryoProtect v2 - Staging Environment Configuration

This module defines the staging environment configuration for CryoProtect v2.
It extends the BaseConfig class from the main config module.
"""

from config import BaseConfig

class StagingConfig(BaseConfig):
    """
    Configuration for the staging environment.
    
    This class overrides the BaseConfig with staging-specific settings.
    """
    # Flask configuration
    DEBUG: bool = False
    TESTING: bool = False
    
    # API configuration
    API_TITLE: str = "CryoProtect Analyzer API (Staging)"
    API_VERSION: str = "v1"
    OPENAPI_VERSION: str = "3.0.2"
    
    # Logging configuration
    LOG_LEVEL: str = "INFO"
    
    # Security configuration
    CORS_ORIGINS: list = ["https://staging.cryoprotect-analyzer.com"]
    
    # Feature flags
    ENABLE_EXPERIMENTAL_FEATURES: bool = True
    
    def _get_env_prefix(self) -> str:
        """
        Get the environment prefix for this configuration class.
        
        Returns:
            The environment prefix for staging-specific variables
        """
        return "STAGING_"