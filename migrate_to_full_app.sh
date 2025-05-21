#!/bin/bash
# Script to migrate from simplified Heroku deployment to the full application

# Color definitions
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
NC='\033[0m' # No Color

echo -e "${YELLOW}Starting migration from simplified to full CryoProtect app on Heroku...${NC}"

# Check if heroku CLI is installed
if ! command -v heroku &> /dev/null; then
    echo -e "${RED}Error: Heroku CLI is not installed. Please install it first.${NC}"
    exit 1
fi

# Check if git is installed
if ! command -v git &> /dev/null; then
    echo -e "${RED}Error: Git is not installed. Please install it first.${NC}"
    exit 1
fi

# Check if we're logged in to Heroku
if ! heroku whoami &> /dev/null; then
    echo -e "${RED}Error: Not logged in to Heroku. Please run 'heroku login' first.${NC}"
    exit 1
fi

# Verify app name
APP_NAME="cryoprotect"
echo -e "${YELLOW}Verifying app ${APP_NAME} exists...${NC}"
if ! heroku apps:info --app $APP_NAME &> /dev/null; then
    echo -e "${RED}Error: App '$APP_NAME' not found. Please provide the correct app name.${NC}"
    exit 1
fi

# Create a dedicated branch for the full app deployment
BRANCH_NAME="heroku-full-app"
echo -e "${YELLOW}Creating branch $BRANCH_NAME for full app deployment...${NC}"
git checkout -b $BRANCH_NAME

# Check if we have uncommitted changes
if [[ -n $(git status -s) ]]; then
    echo -e "${YELLOW}You have uncommitted changes. Stashing them...${NC}"
    git stash
fi

# Update Procfile to use the full app
echo -e "${YELLOW}Updating Procfile to use full app.py...${NC}"
echo "web: gunicorn app:app --log-file -" > Procfile
git add Procfile

# Ensure we have essential dependencies in requirements.txt
echo -e "${YELLOW}Checking and updating dependencies in requirements.txt...${NC}"
if ! grep -q "numpy" requirements.txt; then
    echo -e "${YELLOW}Adding scientific packages to requirements.txt...${NC}"
    cat >> requirements.txt << 'EOF'

# Scientific packages
numpy==1.26.3
scipy==1.13.0  
scikit-learn==1.5.0
pandas==2.2.0
matplotlib==3.8.0
seaborn==0.13.0
pillow==10.2.0
openpyxl==3.1.2
xlsxwriter==3.1.9
reportlab==4.0.9
EOF
    git add requirements.txt
fi

# Create enhanced config file for Heroku environment
echo -e "${YELLOW}Creating enhanced config_heroku.py...${NC}"
cat > config_heroku.py << 'EOF'
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
EOF

git add config_heroku.py

# Update app.py to import the Heroku config
echo -e "${YELLOW}Updating app.py to import Heroku config...${NC}"
if ! grep -q "import config_heroku" app.py; then
    # Add import at the top, right after the existing imports
    sed -i '1s/^/# Import Heroku config\nimport config_heroku\n\n/' app.py
    git add app.py
fi

# Create enhanced database initialization script for Heroku
echo -e "${YELLOW}Creating enhanced setup_database.py for database initialization...${NC}"
cat > setup_database.py << 'EOF'
#!/usr/bin/env python3
"""
Database initialization script for Heroku deployment.
This script is run during the release phase to ensure the database is properly set up.
"""

import os
import sys
import logging
import time
from urllib.parse import urlparse
import subprocess

try:
    import psycopg2
    from psycopg2 import sql
except ImportError:
    # Try to install psycopg2 if it's not available
    subprocess.check_call([sys.executable, "-m", "pip", "install", "psycopg2-binary"])
    import psycopg2
    from psycopg2 import sql

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[logging.StreamHandler(sys.stdout)]
)
logger = logging.getLogger('setup_database')

def wait_for_database(connection_params, max_retries=5, retry_delay=3):
    """Wait for the database to be available."""
    retry_count = 0
    
    while retry_count < max_retries:
        try:
            logger.info(f"Attempting to connect to database (attempt {retry_count + 1}/{max_retries})...")
            
            # Connect to the database
            conn = psycopg2.connect(**connection_params)
            conn.close()
            
            logger.info("Database connection successful!")
            return True
        except psycopg2.OperationalError as e:
            logger.warning(f"Database connection failed: {e}")
            retry_count += 1
            if retry_count < max_retries:
                logger.info(f"Retrying in {retry_delay} seconds...")
                time.sleep(retry_delay)
    
    logger.error("Failed to connect to the database after maximum retries.")
    return False

def set_environment_variables(connection_params):
    """Set environment variables for database utilities."""
    os.environ['SUPABASE_DB_HOST'] = connection_params['host']
    os.environ['SUPABASE_DB_PORT'] = str(connection_params['port'])
    os.environ['SUPABASE_DB_NAME'] = connection_params['database']
    os.environ['SUPABASE_DB_USER'] = connection_params['user']
    os.environ['SUPABASE_DB_PASSWORD'] = connection_params['password']
    
    logger.info(f"Database connection info: {connection_params['host']}:{connection_params['port']}/{connection_params['database']}")

def check_database_structure(connection_params):
    """Check if the database has the required tables."""
    required_tables = [
        'molecules',
        'property_types',
        'molecular_properties',
        'mixtures',
        'mixture_components'
    ]
    
    try:
        # Connect to the database
        conn = psycopg2.connect(**connection_params)
        cursor = conn.cursor()
        
        # Query for existing tables
        cursor.execute("SELECT table_name FROM information_schema.tables WHERE table_schema = 'public';")
        existing_tables = [row[0] for row in cursor.fetchall()]
        
        # Check if all required tables exist
        missing_tables = [table for table in required_tables if table not in existing_tables]
        
        cursor.close()
        conn.close()
        
        if missing_tables:
            logger.warning(f"Missing tables: {', '.join(missing_tables)}")
            return False
        else:
            logger.info("All required tables exist.")
            return True
    except Exception as e:
        logger.error(f"Error checking database structure: {e}")
        return False

def run_migrations():
    """Run database migrations if they exist."""
    if os.path.exists('migrations'):
        logger.info("Applying database migrations...")
        try:
            # Run migrations using apply_migrations.py if it exists
            if os.path.exists('apply_migrations.py'):
                subprocess.run([sys.executable, 'apply_migrations.py'], check=True)
                logger.info("Migrations applied successfully")
            # If migrations directory exists but no apply script, try to apply SQL files
            else:
                migrate_dir_path = os.path.join(os.getcwd(), 'migrations')
                migration_files = sorted([f for f in os.listdir(migrate_dir_path) if f.endswith('.sql')])
                
                if migration_files:
                    logger.info(f"Found {len(migration_files)} migration SQL files.")
                    for sql_file in migration_files:
                        file_path = os.path.join(migrate_dir_path, sql_file)
                        logger.info(f"Applying migration file: {sql_file}")
                        try:
                            # Use apply_file_sql.py if it exists
                            if os.path.exists('apply_file_sql.py'):
                                subprocess.run([sys.executable, 'apply_file_sql.py', file_path], check=True)
                            # Otherwise use psql directly if DATABASE_URL is available
                            elif os.environ.get('DATABASE_URL'):
                                cmd = f'psql {os.environ.get("DATABASE_URL")} -f {file_path}'
                                subprocess.run(cmd, shell=True, check=True)
                        except subprocess.CalledProcessError as e:
                            logger.error(f"Error applying migration file {sql_file}: {e}")
                else:
                    logger.warning("No SQL migration files found in migrations directory.")
        except subprocess.CalledProcessError as e:
            logger.error(f"Error applying migrations: {e}")
            # Don't exit on error, try to continue with other setup
    else:
        logger.info("No migrations directory found. Skipping migrations.")

def setup_database():
    """Set up the database with required tables and initial data."""
    logger.info("Starting database setup for Heroku deployment...")
    
    # Get database URL from environment
    database_url = os.environ.get('DATABASE_URL')
    if not database_url:
        logger.error("DATABASE_URL environment variable not set.")
        return False
    
    # Parse the URL to get connection parameters
    parsed_url = urlparse(database_url)
    connection_params = {
        'host': parsed_url.hostname,
        'port': parsed_url.port or 5432,
        'database': parsed_url.path.lstrip('/'),
        'user': parsed_url.username,
        'password': parsed_url.password
    }
    
    # Set environment variables for database utilities
    set_environment_variables(connection_params)
    
    # Wait for the database to be available
    if not wait_for_database(connection_params):
        return False
    
    # Run database migrations
    run_migrations()
    
    # Check if required tables exist
    if check_database_structure(connection_params):
        logger.info("Database structure is already set up.")
    else:
        # If migrations didn't create tables, try running population scripts
        if os.path.exists('populate_database.py'):
            logger.info("Running database population script...")
            try:
                subprocess.run([sys.executable, 'populate_database.py'], check=True)
                logger.info("Database population completed successfully")
            except subprocess.CalledProcessError as e:
                logger.error(f"Error populating database: {e}")
                return False
    
    # Run Supabase RLS scripts if they exist
    if os.path.exists('apply_rls_simple.py'):
        logger.info("Applying RLS policies...")
        try:
            subprocess.run([sys.executable, 'apply_rls_simple.py'], check=True)
            logger.info("RLS policies applied successfully")
        except subprocess.CalledProcessError as e:
            logger.error(f"Error applying RLS policies: {e}")
            # Don't exit on error, try to continue
    
    logger.info("Database setup completed successfully.")
    return True

if __name__ == '__main__':
    success = setup_database()
    if not success:
        logger.error("Database setup failed.")
        sys.exit(1)
    sys.exit(0)
EOF

git add setup_database.py

# Update the release command in Procfile to run the setup script
echo -e "${YELLOW}Adding release command to Procfile...${NC}"
sed -i '1s/^/release: python setup_database.py\n/' Procfile
git add Procfile

# Commit the changes
echo -e "${YELLOW}Committing changes...${NC}"
git commit -m "Prepare full app for Heroku deployment"

# Push to Heroku
echo -e "${YELLOW}Pushing to Heroku...${NC}"
git push heroku $BRANCH_NAME:master --force

# Scale dynos
echo -e "${YELLOW}Scaling web dyno...${NC}"
heroku ps:scale web=1 --app $APP_NAME

# Check the logs for any errors
echo -e "${YELLOW}Checking logs for errors...${NC}"
heroku logs --tail --app $APP_NAME

echo -e "${GREEN}Migration to full app completed!${NC}"
echo -e "${YELLOW}You can verify it's working by visiting:${NC}"
echo -e "${GREEN}https://$APP_NAME.herokuapp.com/health${NC}"
echo -e "${YELLOW}If there are issues, check the logs with:${NC}"
echo -e "${GREEN}heroku logs --tail --app $APP_NAME${NC}"