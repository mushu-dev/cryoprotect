"""
Heroku-specific configuration for CryoProtect.
This allows the app to utilize Heroku's environment variables and services.
"""

import os
from urllib.parse import urlparse
from config import BaseConfig, config_classes

class HerokuConfig(BaseConfig):
    """Heroku configuration settings with automatic environment variable parsing."""
    
    DEBUG = False
    TESTING = False
    
    def __init__(self):
        """Initialize Heroku config with environment variables from Heroku."""
        super().__init__()
        
        # Load configuration from Heroku environment variables
        self.HEROKU = True
        self.DEBUG = os.environ.get('DEBUG', 'False').lower() == 'true'
        self.TESTING = os.environ.get('TESTING', 'False').lower() == 'true'
        
        # Parse DATABASE_URL for PostgreSQL connection
        database_url = os.environ.get('DATABASE_URL')
        if database_url:
            # Heroku provides PostgreSQL URLs in the format:
            # postgres://username:password@host:port/database
            # Parse this URL to extract connection parameters for Supabase
            parsed_url = urlparse(database_url)
            self.DB_HOST = parsed_url.hostname
            self.DB_PORT = parsed_url.port or 5432
            self.DB_NAME = parsed_url.path.lstrip('/')
            self.DB_USER = parsed_url.username
            self.DB_PASSWORD = parsed_url.password
        
        # Set Supabase configuration from environment variables
        self.SUPABASE_URL = os.environ.get('SUPABASE_URL')
        self.SUPABASE_KEY = os.environ.get('SUPABASE_KEY')
        self.SUPABASE_SERVICE_KEY = os.environ.get('SUPABASE_SERVICE_KEY')
        
        # Redis configuration for rate limiting and caching
        redis_url = os.environ.get('REDIS_URL')
        if redis_url:
            self.REDIS_URL = redis_url
            self.RATE_LIMIT_STORAGE_URL = redis_url
        
        # Enable connection pooling with reasonable defaults for Heroku
        self.SUPABASE_CONNECTION_POOL_ENABLED = True
        self.SUPABASE_MIN_CONNECTIONS = int(os.environ.get('SUPABASE_MIN_CONNECTIONS', 2))
        self.SUPABASE_MAX_CONNECTIONS = int(os.environ.get('SUPABASE_MAX_CONNECTIONS', 10))
        self.SUPABASE_CONNECTION_TIMEOUT = int(os.environ.get('SUPABASE_CONNECTION_TIMEOUT', 30))
        
        # API configuration
        self.API_TITLE = os.environ.get('API_TITLE', 'CryoProtect API')
        self.API_VERSION = os.environ.get('API_VERSION', '1.0.0')
        self.OPENAPI_VERSION = os.environ.get('OPENAPI_VERSION', '3.0.2')
        
        # RDKit configuration - defaults to disabled on Heroku
        self.RDKIT_ENABLED = os.environ.get('RDKIT_ENABLED', 'False').lower() == 'true'
        self.RDKIT_FALLBACK = os.environ.get('RDKIT_FALLBACK', 'True').lower() == 'true'
        self.RDKIT_SERVICE_URL = os.environ.get('RDKIT_SERVICE_URL', None)
    
    def _get_env_prefix(self) -> str:
        return ""

# Add the Heroku config to the config classes
config_classes['heroku'] = HerokuConfig

# If running on Heroku, automatically use the Heroku config
if os.environ.get('DYNO') and not os.environ.get('FLASK_ENV'):
    os.environ['FLASK_ENV'] = 'heroku'

# Make the Heroku config available for direct import
active_heroku_config = HerokuConfig()
