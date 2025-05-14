"""
Minimal configuration file for Heroku deployment.
This enables the app to use the DATABASE_URL and other environment variables.
"""
import os
from config import BaseConfig, config_classes

class HerokuConfig(BaseConfig):
    """Heroku configuration that automatically picks up environment variables."""
    
    DEBUG = False
    TESTING = False
    
    def _get_env_prefix(self) -> str:
        return ""

# Add the Heroku config to the config classes
config_classes['heroku'] = HerokuConfig

# If FLASK_ENV is not set but we're on Heroku, use HerokuConfig
if os.environ.get('DYNO') and not os.environ.get('FLASK_ENV'):
    os.environ['FLASK_ENV'] = 'heroku'
