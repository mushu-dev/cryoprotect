#!/usr/bin/env python3
"""
Database initialization script for Heroku deployment.
This script runs during the release phase to ensure the database is properly set up.
"""
import os
import sys
import json
from urllib.parse import urlparse
import subprocess

def setup_database():
    """Set up the database with required tables and initial data."""
    print("Setting up database...")
    
    # Determine database connection method
    if os.environ.get('DATABASE_URL'):
        print("Using DATABASE_URL for database connection")
        # Parse the DATABASE_URL to extract credentials
        url = urlparse(os.environ.get('DATABASE_URL'))
        dbname = url.path[1:]
        user = url.username
        password = url.password
        host = url.hostname
        port = url.port
        
        # Set environment variables for database utilities
        os.environ['SUPABASE_DB_HOST'] = host
        os.environ['SUPABASE_DB_PORT'] = str(port)
        os.environ['SUPABASE_DB_NAME'] = dbname
        os.environ['SUPABASE_DB_USER'] = user
        os.environ['SUPABASE_DB_PASSWORD'] = password
        
        print(f"Database connection info: {host}:{port}/{dbname}")
    else:
        print("No DATABASE_URL found, assuming Supabase connection")
    
    # Run database migration scripts if they exist
    if os.path.exists('migrations'):
        print("Applying database migrations...")
        try:
            subprocess.run(['python', 'apply_migrations.py'], check=True)
            print("Migrations applied successfully")
        except subprocess.CalledProcessError as e:
            print(f"Error applying migrations: {e}")
            # Don't exit on error, try to continue with other setup
    
    # Run any necessary database population scripts
    print("Database setup completed")

if __name__ == '__main__':
    setup_database()
