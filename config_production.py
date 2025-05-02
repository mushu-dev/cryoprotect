"""
CryoProtect v2 - Production Environment Configuration

This module defines the production environment configuration for CryoProtect v2.
It extends the BaseConfig class from the main config module.
"""

from config import BaseConfig

class ProductionConfig(BaseConfig):
    """
    Configuration for the production environment.
    
    This class overrides the BaseConfig with production-specific settings.
    """
    # Flask configuration
    DEBUG: bool = False
    TESTING: bool = False
    
    # API configuration
    API_TITLE: str = "CryoProtect Analyzer API"
    API_VERSION: str = "v1"
    OPENAPI_VERSION: str = "3.0.2"
    
    # Logging configuration
    LOG_LEVEL: str = "WARNING"
    
    # Security configuration
    CORS_ORIGINS: list = ["https://cryoprotect-analyzer.com"]
    
    # Feature flags
    ENABLE_EXPERIMENTAL_FEATURES: bool = False
    
    # Performance optimization
    CACHE_TYPE: str = "redis"
    
    def _get_env_prefix(self) -> str:
        """
        Get the environment prefix for this configuration class.
        
        Returns:
            The environment prefix for production-specific variables
        """
        return "PRODUCTION_"