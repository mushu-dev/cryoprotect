#!/bin/bash
# Deploy CORS configuration to Heroku app
# This script updates the Heroku app with the CORS configuration

set -e

# Configuration variables with defaults
HEROKU_APP_NAME=${HEROKU_APP_NAME:-cryoprotect}
DEPLOY_BRANCH=${DEPLOY_BRANCH:-main}

echo "==== CryoProtect CORS Configuration Deployment ===="
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
heroku config:get ALLOWED_ORIGINS --app $HEROKU_APP_NAME > /dev/null || {
    echo "Setting ALLOWED_ORIGINS environment variable..."
    heroku config:set ALLOWED_ORIGINS="https://cryoprotect.netlify.app,https://rdkit.cryoprotect.app,https://dynamic-mink-63.convex.cloud" --app $HEROKU_APP_NAME
}

heroku config:get CONVEX_URL --app $HEROKU_APP_NAME > /dev/null || {
    echo "Setting CONVEX_URL environment variable..."
    heroku config:set CONVEX_URL="https://dynamic-mink-63.convex.cloud" --app $HEROKU_APP_NAME
}

heroku config:get FRONTEND_URL --app $HEROKU_APP_NAME > /dev/null || {
    echo "Setting FRONTEND_URL environment variable..."
    heroku config:set FRONTEND_URL="https://cryoprotect.netlify.app" --app $HEROKU_APP_NAME
}

heroku config:get RDKIT_SERVICE_URL --app $HEROKU_APP_NAME > /dev/null || {
    echo "Setting RDKIT_SERVICE_URL environment variable..."
    heroku config:set RDKIT_SERVICE_URL="https://rdkit.cryoprotect.app" --app $HEROKU_APP_NAME
}

# Create a simple deployment to update the Heroku app
echo "Detecting which app version is currently deployed..."
if [ -f "simple_app.py" ]; then
    echo "Found simple_app.py - preparing for deployment"
    
    # Create temporary files to patch the application
    TEMP_DIR=$(mktemp -d)
    mkdir -p "$TEMP_DIR/api"
    
    # Copy CORS configuration to temporary directory
    cp api/cors_config.py "$TEMP_DIR/api/"
    
    # Create a modified version of simple_app.py that uses cors_config.py
    cat > "$TEMP_DIR/simple_app.py" << 'EOL'
#!/usr/bin/env python3
"""
Simplified CryoProtect API for Heroku deployment without RDKit dependencies.
This version uses the standardized CORS configuration from cors_config.py
"""

import os
import logging
from flask import Flask, jsonify, request
import psycopg2
from urllib.parse import urlparse
from api.cors_config import configure_cors

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger('simple_app')

app = Flask(__name__)

# Apply standard CORS configuration from cors_config.py
app = configure_cors(app)

# Get database URL from environment
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
    logger.error("DATABASE_URL environment variable not set.")

def get_db_connection():
    """Get a connection to the database."""
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
    return jsonify({
        'status': 'ok',
        'message': 'CryoProtect API is running with enhanced CORS support',
        'version': '1.0.1'
    })

@app.route('/api/connect')
def api_connect():
    """API connectivity endpoint for testing frontend-backend communication."""
    origin = request.headers.get('Origin', 'Unknown')
    logger.info(f"Connectivity test from origin: {origin}")
    
    return jsonify({
        'status': 'success',
        'message': 'Connected successfully to CryoProtect API',
        'origin': origin,
        'cors': 'enabled',
        'environment': os.environ.get('FLASK_ENV', 'production'),
        'server': 'Heroku'
    })

@app.route('/health')
def health_check():
    """Health check endpoint."""
    try:
        # Check database connection
        conn = get_db_connection()
        cursor = conn.cursor()
        cursor.execute("SELECT 1")
        cursor.close()
        conn.close()
        
        return jsonify({
            'status': 'ok',
            'database': 'connected',
            'timestamp': 'now'
        })
    except Exception as e:
        logger.error(f"Health check failed: {str(e)}")
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
        
        # Connect to the database
        conn = get_db_connection()
        cursor = conn.cursor()
        
        # Get molecules with pagination
        cursor.execute(
            """
            SELECT id, name, smiles, pubchem_cid, molecular_formula, molecular_weight
            FROM molecules
            ORDER BY id
            LIMIT %s OFFSET %s
            """,
            (limit, offset)
        )
        rows = cursor.fetchall()
        
        # Get total count
        cursor.execute("SELECT COUNT(*) FROM molecules")
        total = cursor.fetchone()[0]
        
        # Close cursor and connection
        cursor.close()
        conn.close()
        
        # Format results
        molecules = []
        for row in rows:
            molecules.append({
                'id': row[0],
                'name': row[1],
                'smiles': row[2],
                'pubchem_cid': row[3],
                'molecular_formula': row[4],
                'molecular_weight': row[5]
            })
        
        return jsonify({
            'status': 'success',
            'data': molecules,
            'pagination': {
                'total': total,
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
        # Connect to the database
        conn = get_db_connection()
        cursor = conn.cursor()
        
        # Get the molecule
        cursor.execute(
            """
            SELECT id, name, smiles, pubchem_cid, molecular_formula, molecular_weight
            FROM molecules
            WHERE id = %s
            """,
            (molecule_id,)
        )
        row = cursor.fetchone()
        
        # Close cursor and connection
        cursor.close()
        conn.close()
        
        if not row:
            return jsonify({
                'status': 'error',
                'message': f'Molecule with ID {molecule_id} not found'
            }), 404
        
        # Format result
        molecule = {
            'id': row[0],
            'name': row[1],
            'smiles': row[2],
            'pubchem_cid': row[3],
            'molecular_formula': row[4],
            'molecular_weight': row[5]
        }
        
        return jsonify({
            'status': 'success',
            'data': molecule
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
    git commit -m "Deploy CORS configuration for backend integration"
    
    echo "Deploying CORS configuration to Heroku..."
    git push -f "https://git.heroku.com/$HEROKU_APP_NAME.git" HEAD:main
    
    # Clean up
    cd -
    rm -rf "$TEMP_DIR"
else
    echo "❌ simple_app.py not found. Please run this script from the project root directory."
    exit 1
fi

echo "✅ CORS configuration deployed to Heroku app: $HEROKU_APP_NAME"
echo "Run the test-connection.js script to verify the configuration."