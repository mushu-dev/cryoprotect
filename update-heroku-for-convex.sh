#!/bin/bash
# Update Heroku app to use Convex
# This script patches the Heroku app to use Convex instead of Supabase

set -e

# Configuration variables with defaults
HEROKU_APP_NAME=${HEROKU_APP_NAME:-cryoprotect}

echo "==== CryoProtect Convex Database Integration ===="
echo "Target Heroku app: $HEROKU_APP_NAME"
echo

# Verify Heroku CLI is available
if ! command -v heroku &> /dev/null; then
    echo "❌ Heroku CLI is not installed or not available in PATH"
    echo "Please install the Heroku CLI: https://devcenter.heroku.com/articles/heroku-cli"
    exit 1
fi

# Verify logged in to Heroku
echo "Verifying Heroku authentication..."
heroku auth:whoami || {
    echo "❌ Not logged in to Heroku"
    echo "Please run 'heroku login' to authenticate"
    exit 1
}

# Verify app exists
echo "Verifying app exists..."
heroku apps:info --app $HEROKU_APP_NAME || {
    echo "❌ App $HEROKU_APP_NAME does not exist or you don't have access"
    exit 1
}

# Check if we need to set environment variables
echo "Checking for required environment variables..."
heroku config:get CONVEX_URL --app $HEROKU_APP_NAME > /dev/null || {
    echo "Setting CONVEX_URL environment variable..."
    heroku config:set CONVEX_URL="https://dynamic-mink-63.convex.cloud" --app $HEROKU_APP_NAME
}

heroku config:get USE_CONVEX --app $HEROKU_APP_NAME > /dev/null || {
    echo "Setting USE_CONVEX environment variable..."
    heroku config:set USE_CONVEX="true" --app $HEROKU_APP_NAME
}

heroku config:get CONVEX_DEPLOYMENT_KEY --app $HEROKU_APP_NAME > /dev/null || {
    echo "Setting CONVEX_DEPLOYMENT_KEY environment variable..."
    heroku config:set CONVEX_DEPLOYMENT_KEY="prod_DsfZ1kX1FcUx2i5bI9ovW9GsCsSGGzPt" --app $HEROKU_APP_NAME
}

# Create a temporary deployment with Convex integration
echo "Creating temporary deployment with Convex integration..."
TEMP_DIR=$(mktemp -d)
mkdir -p "$TEMP_DIR/api"
mkdir -p "$TEMP_DIR/database"

# Copy CORS configuration to temporary directory
cp api/cors_config.py "$TEMP_DIR/api/"

# Create Convex adapter in temporary directory
cat > "$TEMP_DIR/database/convex_adapter.py" << 'EOL'
"""
Convex database adapter for CryoProtect.

This module provides a compatibility layer between the existing Supabase-based
API and the new Convex database backend.
"""

import os
import json
import requests
import logging
from urllib.parse import urljoin

logger = logging.getLogger(__name__)

class ConvexAdapter:
    """
    Adapter class that provides Supabase-like interface for Convex.
    This allows for gradual migration of the codebase.
    """
    
    def __init__(self, url=None, key=None):
        """
        Initialize the Convex adapter.
        
        Args:
            url (str): The Convex URL. If not provided, it's read from CONVEX_URL env var.
            key (str): The Convex deployment key. If not provided, it's read from CONVEX_DEPLOYMENT_KEY env var.
        """
        self.url = url or os.environ.get('CONVEX_URL', 'https://dynamic-mink-63.convex.cloud')
        self.key = key or os.environ.get('CONVEX_DEPLOYMENT_KEY', '')
        
        # Ensure the URL ends with a slash for proper joining
        if not self.url.endswith('/'):
            self.url = self.url + '/'
            
        self.headers = {
            'Content-Type': 'application/json',
            'Authorization': f'Bearer {self.key}' if self.key else None
        }
        
        logger.info("Initialized Convex adapter with URL: %s", self.url)
    
    def table(self, table_name):
        """
        Get a reference to a Convex table (collection).
        
        Args:
            table_name (str): The name of the table.
            
        Returns:
            TableAdapter: An adapter for the specified table.
        """
        return TableAdapter(self, table_name)
    
    def auth(self):
        """
        Get a reference to the Convex auth methods.
        
        Returns:
            AuthAdapter: An adapter for authentication methods.
        """
        return AuthAdapter(self)
    
    def execute_query(self, action, path, params=None):
        """
        Execute a query against Convex.
        
        Args:
            action (str): The HTTP method (GET, POST, PUT, DELETE).
            path (str): The API path.
            params (dict): The parameters to send with the request.
            
        Returns:
            dict: The response from Convex.
        """
        url = urljoin(self.url, path)
        
        try:
            if action == 'GET':
                response = requests.get(url, params=params, headers=self.headers)
            elif action == 'POST':
                response = requests.post(url, json=params, headers=self.headers)
            elif action == 'PUT':
                response = requests.put(url, json=params, headers=self.headers)
            elif action == 'DELETE':
                response = requests.delete(url, json=params, headers=self.headers)
            else:
                raise ValueError(f"Unsupported action: {action}")
            
            response.raise_for_status()
            return response.json()
        except Exception as e:
            logger.error("Error executing Convex query: %s", str(e))
            raise

class TableAdapter:
    """Adapter for Convex tables to provide Supabase-like interface."""
    
    def __init__(self, adapter, table_name):
        """
        Initialize the table adapter.
        
        Args:
            adapter (ConvexAdapter): The parent Convex adapter.
            table_name (str): The name of the table.
        """
        self.adapter = adapter
        self.table_name = table_name
        self.filters = {}
        self.order_clauses = []
        self.limit_value = None
        self.offset_value = None
        
    def select(self, columns="*"):
        """
        Select columns from the table.
        
        Args:
            columns (str): Comma-separated column names or "*" for all columns.
            
        Returns:
            TableAdapter: Self, for method chaining.
        """
        self.columns = columns
        return self
    
    def eq(self, column, value):
        """
        Add an equality filter.
        
        Args:
            column (str): The column name.
            value: The value to compare against.
            
        Returns:
            TableAdapter: Self, for method chaining.
        """
        self.filters[column] = value
        return self
    
    def neq(self, column, value):
        """
        Add a not-equal filter.
        
        Args:
            column (str): The column name.
            value: The value to compare against.
            
        Returns:
            TableAdapter: Self, for method chaining.
        """
        self.filters[f"{column}_neq"] = value
        return self
    
    def gt(self, column, value):
        """
        Add a greater-than filter.
        
        Args:
            column (str): The column name.
            value: The value to compare against.
            
        Returns:
            TableAdapter: Self, for method chaining.
        """
        self.filters[f"{column}_gt"] = value
        return self
    
    def gte(self, column, value):
        """
        Add a greater-than-or-equal filter.
        
        Args:
            column (str): The column name.
            value: The value to compare against.
            
        Returns:
            TableAdapter: Self, for method chaining.
        """
        self.filters[f"{column}_gte"] = value
        return self
    
    def lt(self, column, value):
        """
        Add a less-than filter.
        
        Args:
            column (str): The column name.
            value: The value to compare against.
            
        Returns:
            TableAdapter: Self, for method chaining.
        """
        self.filters[f"{column}_lt"] = value
        return self
    
    def lte(self, column, value):
        """
        Add a less-than-or-equal filter.
        
        Args:
            column (str): The column name.
            value: The value to compare against.
            
        Returns:
            TableAdapter: Self, for method chaining.
        """
        self.filters[f"{column}_lte"] = value
        return self
    
    def order(self, column, order="asc"):
        """
        Add an ordering clause.
        
        Args:
            column (str): The column name.
            order (str): Either "asc" or "desc".
            
        Returns:
            TableAdapter: Self, for method chaining.
        """
        self.order_clauses.append((column, order))
        return self
    
    def limit(self, limit):
        """
        Set a limit on the number of results.
        
        Args:
            limit (int): The maximum number of results to return.
            
        Returns:
            TableAdapter: Self, for method chaining.
        """
        self.limit_value = limit
        return self
    
    def offset(self, offset):
        """
        Set an offset for the results.
        
        Args:
            offset (int): The offset for results.
            
        Returns:
            TableAdapter: Self, for method chaining.
        """
        self.offset_value = offset
        return self
    
    def execute(self):
        """
        Execute the query.
        
        Returns:
            Response: A Supabase-like response object.
        """
        # Construct query parameters
        params = {
            "table": self.table_name,
            "filters": self.filters,
        }
        
        if hasattr(self, 'columns'):
            params["columns"] = self.columns
            
        if self.order_clauses:
            params["order"] = self.order_clauses
            
        if self.limit_value is not None:
            params["limit"] = self.limit_value
            
        if self.offset_value is not None:
            params["offset"] = self.offset_value
        
        # Execute the query and convert to Supabase-like response
        try:
            result = self.adapter.execute_query('POST', 'api/query', params)
            return Response(result.get('data', []), result.get('error'))
        except Exception as e:
            logger.error("Error executing Convex query: %s", str(e))
            return Response([], str(e))
    
    def insert(self, data):
        """
        Insert data into the table.
        
        Args:
            data (dict or list): The data to insert.
            
        Returns:
            Response: A Supabase-like response object.
        """
        params = {
            "table": self.table_name,
            "data": data
        }
        
        try:
            result = self.adapter.execute_query('POST', 'api/insert', params)
            return Response(result.get('data', []), result.get('error'))
        except Exception as e:
            logger.error("Error inserting into Convex: %s", str(e))
            return Response([], str(e))
    
    def update(self, data):
        """
        Update data in the table.
        
        Args:
            data (dict): The data to update.
            
        Returns:
            Response: A Supabase-like response object.
        """
        params = {
            "table": self.table_name,
            "data": data,
            "filters": self.filters
        }
        
        try:
            result = self.adapter.execute_query('POST', 'api/update', params)
            return Response(result.get('data', []), result.get('error'))
        except Exception as e:
            logger.error("Error updating in Convex: %s", str(e))
            return Response([], str(e))
    
    def delete(self):
        """
        Delete data from the table.
        
        Returns:
            Response: A Supabase-like response object.
        """
        params = {
            "table": self.table_name,
            "filters": self.filters
        }
        
        try:
            result = self.adapter.execute_query('POST', 'api/delete', params)
            return Response(result.get('data', []), result.get('error'))
        except Exception as e:
            logger.error("Error deleting from Convex: %s", str(e))
            return Response([], str(e))

class AuthAdapter:
    """Adapter for Convex auth to provide Supabase-like interface."""
    
    def __init__(self, adapter):
        """
        Initialize the auth adapter.
        
        Args:
            adapter (ConvexAdapter): The parent Convex adapter.
        """
        self.adapter = adapter
    
    def sign_in_with_password(self, credentials):
        """
        Sign in with email and password.
        
        Args:
            credentials (dict): Dict with email and password.
            
        Returns:
            Response: A Supabase-like response object.
        """
        params = {
            "email": credentials.get('email'),
            "password": credentials.get('password')
        }
        
        try:
            result = self.adapter.execute_query('POST', 'api/auth/signin', params)
            return Response(result.get('data', {}), result.get('error'))
        except Exception as e:
            logger.error("Error signing in with Convex: %s", str(e))
            return Response({}, str(e))
    
    def sign_up(self, credentials):
        """
        Sign up with email and password.
        
        Args:
            credentials (dict): Dict with email and password.
            
        Returns:
            Response: A Supabase-like response object.
        """
        params = {
            "email": credentials.get('email'),
            "password": credentials.get('password')
        }
        
        try:
            result = self.adapter.execute_query('POST', 'api/auth/signup', params)
            return Response(result.get('data', {}), result.get('error'))
        except Exception as e:
            logger.error("Error signing up with Convex: %s", str(e))
            return Response({}, str(e))
    
    def sign_out(self):
        """
        Sign out the current user.
        
        Returns:
            Response: A Supabase-like response object.
        """
        try:
            result = self.adapter.execute_query('POST', 'api/auth/signout', {})
            return Response(result.get('data', {}), result.get('error'))
        except Exception as e:
            logger.error("Error signing out with Convex: %s", str(e))
            return Response({}, str(e))

class Response:
    """A Supabase-like response object."""
    
    def __init__(self, data, error=None):
        """
        Initialize the response.
        
        Args:
            data: The response data.
            error: Any error that occurred.
        """
        self.data = data
        self.error = error

def create_client(url=None, key=None, use_convex=None):
    """
    Create a Convex client based on environment configuration.
    
    Args:
        url (str): The Convex URL.
        key (str): The Convex deployment key.
        use_convex (bool): Force use of Convex if True, Supabase if False.
        
    Returns:
        ConvexAdapter or supabase.Client: The appropriate client.
    """
    # Determine whether to use Convex or fall back to Supabase
    use_convex_env = os.environ.get('USE_CONVEX', '').lower() in ('true', 'yes', '1')
    should_use_convex = use_convex if use_convex is not None else use_convex_env
    
    if should_use_convex:
        logger.info("Using Convex adapter for database operations")
        return ConvexAdapter(url, key)
    else:
        # Fall back to Supabase
        logger.info("Using Supabase for database operations")
        try:
            from supabase import create_client as create_supabase_client
            
            supabase_url = os.environ.get('SUPABASE_URL', '')
            supabase_key = os.environ.get('SUPABASE_KEY', '')
            
            return create_supabase_client(supabase_url, supabase_key)
        except ImportError:
            logger.error("Supabase client library not available")
            return None
EOL

# Create database factory in temporary directory
cat > "$TEMP_DIR/database/db_factory.py" << 'EOL'
"""
Database client factory for CryoProtect.

This module provides a factory function for creating the appropriate database
client (Supabase or Convex) based on configuration.
"""

import os
import logging

logger = logging.getLogger(__name__)

def get_db_client(force_convex=None):
    """
    Get the appropriate database client based on configuration.
    
    Args:
        force_convex (bool): If provided, forces use of Convex (True) or Supabase (False).
        
    Returns:
        Object: Either a Convex adapter or Supabase client.
    """
    # Determine whether to use Convex or Supabase
    use_convex_env = os.environ.get('USE_CONVEX', '').lower() in ('true', 'yes', '1')
    should_use_convex = force_convex if force_convex is not None else use_convex_env
    
    if should_use_convex:
        logger.info("Using Convex for database operations")
        from database.convex_adapter import create_client
        return create_client()
    else:
        logger.info("Using Supabase for database operations")
        try:
            from supabase import create_client
            
            supabase_url = os.environ.get('SUPABASE_URL', '')
            supabase_key = os.environ.get('SUPABASE_KEY', '')
            
            if not supabase_url or not supabase_key:
                raise ValueError("Supabase URL and key are required when not using Convex")
            
            return create_client(supabase_url, supabase_key)
        except ImportError:
            logger.error("Supabase client library not available")
            return None
EOL

# Create an empty __init__.py file for the database package
touch "$TEMP_DIR/database/__init__.py"

# Create modified simple_app.py that uses our database factory
cat > "$TEMP_DIR/simple_app.py" << 'EOL'
#!/usr/bin/env python3
"""
Simplified CryoProtect API for Heroku deployment without RDKit dependencies.
This version uses the standardized CORS configuration from cors_config.py
and supports both Supabase and Convex databases.
"""

import os
import logging
from flask import Flask, jsonify, request
import psycopg2
from urllib.parse import urlparse
from api.cors_config import configure_cors
from database.db_factory import get_db_client

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger('simple_app')

app = Flask(__name__)

# Apply standard CORS configuration from cors_config.py
app = configure_cors(app)

# Get database URL from environment (for direct Postgres access)
DATABASE_URL = os.environ.get('DATABASE_URL')
if DATABASE_URL:
    # Parse the URL to get connection parameters
    parsed_url = urlparse(DATABASE_URL)
    DB_HOST = parsed_url.hostname
    DB_PORT = parsed_url.port or 5432
    DB_NAME = parsed_url.path.lstrip('/')
    DB_USER = parsed_url.username
    DB_PASSWORD = parsed_url.password
    logger.info(f"Database connection info: {DB_HOST}:{DB_PORT}/{DB_NAME}")
else:
    logger.warning("DATABASE_URL environment variable not set.")

def get_db_connection():
    """Get a connection to the database."""
    if not DATABASE_URL:
        raise ValueError("DATABASE_URL not set")
    
    conn = psycopg2.connect(
        host=DB_HOST,
        port=DB_PORT,
        dbname=DB_NAME,
        user=DB_USER,
        password=DB_PASSWORD
    )
    return conn

@app.route('/')
def index():
    # Check if we're using Convex
    use_convex = os.environ.get('USE_CONVEX', '').lower() in ('true', 'yes', '1')
    db_type = "Convex" if use_convex else "Supabase"
    
    return jsonify({
        'status': 'ok',
        'message': f'CryoProtect API is running with enhanced CORS support and {db_type} integration',
        'version': '1.0.2'
    })

@app.route('/api/connect')
def api_connect():
    """API connectivity endpoint for testing frontend-backend communication."""
    origin = request.headers.get('Origin', 'Unknown')
    logger.info(f"Connectivity test from origin: {origin}")
    
    # Check if we're using Convex
    use_convex = os.environ.get('USE_CONVEX', '').lower() in ('true', 'yes', '1')
    db_type = "Convex" if use_convex else "Supabase"
    
    return jsonify({
        'status': 'success',
        'message': f'Connected successfully to CryoProtect API with {db_type} integration',
        'origin': origin,
        'cors': 'enabled',
        'environment': os.environ.get('FLASK_ENV', 'production'),
        'database': db_type,
        'server': 'Heroku'
    })

@app.route('/health')
def health_check():
    """Health check endpoint."""
    try:
        # Check if we're using Convex or Supabase
        use_convex = os.environ.get('USE_CONVEX', '').lower() in ('true', 'yes', '1')
        
        if use_convex:
            # For Convex, we'll just check if we can create a client
            from database.convex_adapter import create_client
            client = create_client()
            db_status = "Convex client initialized"
        else:
            # For direct Postgres, we'll check the connection
            conn = get_db_connection()
            cursor = conn.cursor()
            cursor.execute("SELECT 1")
            cursor.close()
            conn.close()
            db_status = "Postgres connection verified"
        
        return jsonify({
            'status': 'ok',
            'database': db_status,
            'timestamp': 'now'
        })
    except Exception as e:
        logger.error(f"Health check failed: {str(e)}")
        return jsonify({
            'status': 'error',
            'message': str(e)
        }), 500

@app.route('/api/database/status')
def database_status():
    """Get database status information."""
    try:
        # Check if we're using Convex or Supabase
        use_convex = os.environ.get('USE_CONVEX', '').lower() in ('true', 'yes', '1')
        
        if use_convex:
            db_type = "Convex"
            db_url = os.environ.get('CONVEX_URL', '')
            
            # Mask the deployment key for security
            deployment_key = os.environ.get('CONVEX_DEPLOYMENT_KEY', '')
            masked_key = "***" + deployment_key[-4:] if deployment_key else "Not set"
            
            details = {
                "url": db_url,
                "deployment_key": masked_key
            }
        else:
            db_type = "Supabase"
            
            # Mask the keys for security
            key = os.environ.get('SUPABASE_KEY', '')
            masked_key = "***" + key[-4:] if key else "Not set"
            
            service_key = os.environ.get('SUPABASE_SERVICE_KEY', '')
            masked_service_key = "***" + service_key[-4:] if service_key else "Not set"
            
            details = {
                "url": os.environ.get('SUPABASE_URL', ''),
                "key": masked_key,
                "service_key": masked_service_key,
                "db_host": os.environ.get('SUPABASE_DB_HOST', '')
            }
        
        return jsonify({
            'status': 'success',
            'database_type': db_type,
            'details': details
        })
    except Exception as e:
        logger.error(f"Database status check failed: {str(e)}")
        return jsonify({
            'status': 'error',
            'message': str(e)
        }), 500

@app.route('/api/molecules')
def list_molecules():
    """List molecules endpoint."""
    try:
        # Get query parameters
        limit = request.args.get('limit', 10, type=int)
        offset = request.args.get('offset', 0, type=int)
        
        # Get the appropriate database client
        db = get_db_client()
        
        if db is None:
            return jsonify({
                'status': 'error',
                'message': 'Database client not available'
            }), 500
        
        # Use the unified client interface to query molecules
        response = db.table("molecules").select("*").limit(limit).offset(offset).execute()
        
        if response.error:
            return jsonify({
                'status': 'error',
                'message': str(response.error)
            }), 500
        
        return jsonify({
            'status': 'success',
            'data': response.data,
            'pagination': {
                'limit': limit,
                'offset': offset
            }
        })
    except Exception as e:
        logger.error(f"Error listing molecules: {str(e)}")
        return jsonify({
            'status': 'error',
            'message': str(e)
        }), 500

@app.route('/api/molecules/<int:molecule_id>')
def get_molecule(molecule_id):
    """Get a specific molecule."""
    try:
        # Get the appropriate database client
        db = get_db_client()
        
        if db is None:
            return jsonify({
                'status': 'error',
                'message': 'Database client not available'
            }), 500
        
        # Use the unified client interface to get the molecule
        response = db.table("molecules").select("*").eq("id", molecule_id).execute()
        
        if response.error:
            return jsonify({
                'status': 'error',
                'message': str(response.error)
            }), 500
        
        if not response.data:
            return jsonify({
                'status': 'error',
                'message': f'Molecule with ID {molecule_id} not found'
            }), 404
        
        return jsonify({
            'status': 'success',
            'data': response.data[0]
        })
    except Exception as e:
        logger.error(f"Error getting molecule {molecule_id}: {str(e)}")
        return jsonify({
            'status': 'error',
            'message': str(e)
        }), 500

if __name__ == '__main__':
    port = int(os.environ.get('PORT', 5000))
    app.run(host='0.0.0.0', port=port)
EOL

# Create a requirements.txt for the deployment
cat > "$TEMP_DIR/requirements.txt" << 'EOL'
flask==2.3.3
flask-cors==4.0.0
psycopg2-binary==2.9.9
gunicorn==21.2.0
requests==2.31.0
python-dotenv==1.0.0
EOL

# Create a Procfile
cat > "$TEMP_DIR/Procfile" << 'EOL'
web: gunicorn simple_app:app
EOL

# Create a temporary git repository and push to Heroku
cd "$TEMP_DIR"
git init
git add .
git config --local user.email "deploy@cryoprotect.app"
git config --local user.name "CryoProtect Deployment"
git commit -m "Deploy Convex integration for backend"

echo "Deploying Convex integration to Heroku..."
git push -f "https://git.heroku.com/$HEROKU_APP_NAME.git" HEAD:main

# Clean up
cd -
rm -rf "$TEMP_DIR"

echo "✅ Convex integration deployed to Heroku app: $HEROKU_APP_NAME"
echo "Verify the integration with:"
echo "  curl https://cryoprotect-8030e4025428.herokuapp.com/api/database/status"