"""
Database configuration module for CryoProtect.
This module provides a centralized configuration for all database connections.
"""

import os
from urllib.parse import urlparse
from dotenv import load_dotenv

# Load environment variables from .env file
load_dotenv()

# Supabase project ID
SUPABASE_PROJECT_ID = "tsdlmynydfuypiugmkev"

# Connection options
USE_SSL = True  # Required for Supabase
CONNECTION_TIMEOUT = 30  # seconds
MAX_CONNECTIONS = 10
APPLICATION_NAME = "CryoProtect"

# Production Supabase URL and key
SUPABASE_URL = os.environ.get("SUPABASE_URL")
SUPABASE_KEY = os.environ.get("SUPABASE_KEY")
SUPABASE_SERVICE_KEY = os.environ.get("SUPABASE_SERVICE_KEY")

# Direct Postgres connection details
DB_HOST = os.environ.get("DB_HOST")
DB_PORT = os.environ.get("DB_PORT", "5432")
DB_NAME = os.environ.get("DB_NAME")
DB_USER = os.environ.get("DB_USER")
DB_PASSWORD = os.environ.get("DB_PASSWORD")
DB_URL = os.environ.get("DATABASE_URL")

def get_db_connection_params():
    """
    Get the database connection parameters based on available environment variables.
    Returns a dictionary suitable for use with psycopg2.connect()
    """
    # If full connection URL is provided, parse it
    if DB_URL:
        url = urlparse(DB_URL)
        return {
            "host": url.hostname,
            "port": url.port or "5432",
            "dbname": url.path[1:],
            "user": url.username,
            "password": url.password,
            "application_name": APPLICATION_NAME,
            "sslmode": "require" if USE_SSL else "prefer"
        }
    # Otherwise use individual environment variables
    elif DB_HOST and DB_NAME and DB_USER and DB_PASSWORD:
        return {
            "host": DB_HOST,
            "port": DB_PORT,
            "dbname": DB_NAME,
            "user": DB_USER,
            "password": DB_PASSWORD,
            "application_name": APPLICATION_NAME,
            "sslmode": "require" if USE_SSL else "prefer"
        }
    # For local development with default PostgreSQL config
    else:
        return {
            "host": "localhost",
            "port": "5432",
            "dbname": "postgres",
            "user": "postgres",
            "password": "", 
            "application_name": APPLICATION_NAME,
            "sslmode": "disable"  # Local connections typically don't need SSL
        }

def get_connection_pool_config():
    """Get configuration for connection pooling."""
    return {
        "minconn": 1,
        "maxconn": MAX_CONNECTIONS,
        "dsn": None,  # Will be set at runtime
        "connection_timeout": CONNECTION_TIMEOUT
    }

def is_supabase_configured():
    """Check if Supabase is properly configured."""
    return bool(SUPABASE_URL and SUPABASE_KEY)

def is_direct_db_configured():
    """Check if direct database connection is properly configured."""
    return bool(DB_URL or (DB_HOST and DB_NAME and DB_USER and DB_PASSWORD))

def get_supabase_config():
    """Get Supabase configuration."""
    return {
        "url": SUPABASE_URL,
        "key": SUPABASE_KEY,
        "service_key": SUPABASE_SERVICE_KEY,
        "project_id": SUPABASE_PROJECT_ID
    }

# Export the main connection string for use in scripts
CONNECTION_STRING = DB_URL or f"postgresql://{DB_USER}:{DB_PASSWORD}@{DB_HOST}:{DB_PORT}/{DB_NAME}"