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
    
    # Check if required tables exist
    if check_database_structure(connection_params):
        logger.info("Database structure is already set up.")
    else:
        # Create tables using the Heroku-specific script if it exists
        if os.path.exists('create_heroku_tables.py'):
            logger.info("Creating database tables using Heroku-specific script...")
            try:
                subprocess.run([sys.executable, 'create_heroku_tables.py'], check=True)
                logger.info("Database tables created successfully.")
                return True
            except subprocess.CalledProcessError as e:
                logger.error(f"Error creating database tables: {e}")
                return False
        # Otherwise try existing database population scripts
        elif os.path.exists('populate_database.py'):
            logger.info("Running database population script...")
            try:
                subprocess.run([sys.executable, 'populate_database.py'], check=True)
                logger.info("Database population completed successfully")
            except subprocess.CalledProcessError as e:
                logger.error(f"Error populating database: {e}")
                return False
        else:
            # No table creation scripts available
            logger.error("No database creation scripts found. Create create_heroku_tables.py or populate_database.py")
            return False
    
    logger.info("Database setup completed successfully.")
    return True

if __name__ == '__main__':
    success = setup_database()
    if not success:
        logger.error("Database setup failed.")
        sys.exit(1)
    sys.exit(0)
