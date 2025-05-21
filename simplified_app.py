"""
Simplified CryoProtect app with consolidated database connection.
This version uses the unified database/core module for all database operations.
"""

import os
import sys
import logging
import json
from datetime import datetime
from flask import Flask, jsonify, request
import requests
from dotenv import load_dotenv

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler(sys.stdout)
    ]
)

# Load environment variables
load_dotenv()
logger = logging.getLogger(__name__)

# Import our core database module
try:
    # First try the new Fedora-specific database.core module
    from database.core import (
        test_connection,
        rest_request,
        execute_sql,
        get_table_data,
        count_records,
        get_schema_info
    )

    # Test if the module is working
    db_test = test_connection()
    if db_test['overall']['status'] == 'success':
        logger.info("Database core module successfully loaded")
        HAS_DB_CORE = True
    else:
        logger.warning(f"Database core module test failed: {db_test['overall']['message']}")
        HAS_DB_CORE = False

except (ImportError, AttributeError) as e:
    logger.warning(f"Failed to load new database core module: {str(e)}")
    logger.info("Falling back to compatibility mode")

    try:
        # Legacy path - just set flag to false to ensure app can still run
        HAS_DB_CORE = False
        logger.info("Running in compatibility mode without database core functions")
    except Exception as e2:
        logger.error(f"Failed to initialize compatibility mode: {str(e2)}")
        HAS_DB_CORE = False

# Create Flask app
app = Flask(__name__)

# Get Supabase credentials from environment for backward compatibility
SUPABASE_URL = os.environ.get('SUPABASE_URL')
SUPABASE_KEY = os.environ.get('SUPABASE_KEY')
SUPABASE_SERVICE_KEY = os.environ.get('SUPABASE_SERVICE_KEY')

if not SUPABASE_URL or not SUPABASE_KEY:
    logger.error("Supabase URL or key not configured. Please set SUPABASE_URL and SUPABASE_KEY environment variables.")
    # We'll continue running but warn that Supabase features will not work

@app.route('/')
def index():
    """Root endpoint with basic information"""
    return jsonify({
        'name': 'CryoProtect Simplified App',
        'status': 'running',
        'environment': os.environ.get('FLASK_ENV', 'development'),
        'python_version': sys.version,
        'supabase_configured': bool(SUPABASE_URL and SUPABASE_KEY)
    })

@app.route('/health')
def health():
    """Health check endpoint"""
    return jsonify({
        'status': 'healthy',
        'uptime': 'unknown'  # Would calculate from start time in a full app
    })

@app.route('/api/v1/health/connectivity')
def api_connectivity_check():
    """
    Simple API connectivity test endpoint for verifying frontend-backend connection.
    Returns basic information about the API and connection status.
    """
    try:
        # Get Heroku app name
        heroku_app_name = os.environ.get('HEROKU_APP_NAME', 'cryoprotect-8030e4025428')
        # Get Vercel-specific environment variables or netlify URL
        frontend_url = os.environ.get('VERCEL_FRONTEND_URL', os.environ.get('NETLIFY_URL', 'https://cryoprotect.netlify.app'))
        
        # Return connectivity information
        return jsonify({
            'status': 'connected',
            'api_version': os.environ.get('API_VERSION', '1.0.0'),
            'environment': os.environ.get('FLASK_ENV', 'production'),
            'deployment': {
                'backend': f"https://{heroku_app_name}.herokuapp.com",
                'frontend': frontend_url
            },
            'timestamp': str(datetime.now().isoformat()),
            'cors_enabled': True
        }), 200
    except Exception as e:
        logger.error(f"API connectivity check failed: {str(e)}")
        return jsonify({
            'status': 'error',
            'message': str(e),
            'timestamp': str(datetime.now().isoformat())
        }), 500

@app.route('/env')
def environment():
    """Display environment variables (safe ones only)"""
    safe_vars = {
        'FLASK_ENV': os.environ.get('FLASK_ENV', 'not set'),
        'FLASK_APP': os.environ.get('FLASK_APP', 'not set'),
        'PYTHON_VERSION': sys.version,
        'PYTHON_PATH': sys.executable,
        'PWD': os.environ.get('PWD', 'not set'),
        'HOST': os.environ.get('HOST', 'not set'),
        'PORT': os.environ.get('PORT', 'not set')
    }
    return jsonify(safe_vars)

@app.route('/test-supabase')
def test_supabase():
    """Test Supabase connectivity"""
    if not SUPABASE_URL or not SUPABASE_KEY:
        return jsonify({
            'status': 'error',
            'message': 'Supabase URL or key not configured'
        }), 500
    
    try:
        # Test Supabase REST API
        headers = {
            'apikey': SUPABASE_KEY,
            'Authorization': f'Bearer {SUPABASE_KEY}'
        }
        
        response = requests.get(f"{SUPABASE_URL}/rest/v1/", headers=headers, timeout=5)
        
        return jsonify({
            'status': 'success' if response.status_code < 300 else 'error',
            'status_code': response.status_code,
            'message': 'Supabase connection successful' if response.status_code < 300 else 'Supabase connection failed',
            'response': response.text[:200]  # First 200 chars of response
        })
    except Exception as e:
        logger.error(f'Supabase connection error: {str(e)}')
        return jsonify({
            'status': 'error',
            'message': f'Exception: {str(e)}'
        }), 500

@app.route('/supabase/tables')
def list_tables():
    """List tables in Supabase"""
    if not SUPABASE_URL or not SUPABASE_KEY:
        return jsonify({
            'status': 'error',
            'message': 'Supabase URL or key not configured'
        }), 500
    
    try:
        # Get table information from Supabase
        headers = {
            'apikey': SUPABASE_KEY,
            'Authorization': f'Bearer {SUPABASE_KEY}'
        }
        
        # First, get schema information
        response = requests.get(
            f"{SUPABASE_URL}/rest/v1/",
            headers=headers,
            timeout=5
        )
        
        if response.status_code >= 300:
            return jsonify({
                'status': 'error',
                'message': f'Failed to get schema information: {response.status_code}',
                'response': response.text[:200]
            }), response.status_code
        
        # Parse OpenAPI schema to get table names
        try:
            schema = response.json()
            paths = schema.get('paths', {})
            
            tables = []
            for path in paths:
                # Extract table name from path
                # Typical path format: /rest/v1/table_name
                parts = path.strip('/').split('/')
                if len(parts) >= 3 and parts[0] == 'rest' and parts[1] == 'v1':
                    table_name = parts[2]
                    if table_name and '.' not in table_name and '{' not in table_name:
                        tables.append(table_name)
            
            return jsonify({
                'status': 'success',
                'tables': sorted(list(set(tables)))
            })
        except Exception as e:
            logger.error(f'Error parsing schema: {str(e)}')
            return jsonify({
                'status': 'error',
                'message': f'Error parsing schema: {str(e)}',
                'response': response.text[:200]
            }), 500
    except Exception as e:
        logger.error(f'Supabase error: {str(e)}')
        return jsonify({
            'status': 'error',
            'message': f'Exception: {str(e)}'
        }), 500

@app.route('/dependencies')
def check_dependencies():
    """Check if important dependencies are available"""
    dependencies = {}
    
    # Function to safely get module version
    def get_module_version(module):
        try:
            version = getattr(module, '__version__', None)
            if version:
                return version
            version = getattr(module, 'VERSION', None)
            if version:
                return version
            version = getattr(module, 'version', None)
            if version:
                return version
            version = getattr(module, 'Version', None)
            if version:
                return version
            return "Available (version unknown)"
        except Exception:
            return "Available (version unknown)"
    
    # Check Flask
    dependencies['flask'] = {
        'status': 'available',
        'version': get_module_version(Flask)
    }
    
    # Check requests
    dependencies['requests'] = {
        'status': 'available',
        'version': get_module_version(requests)
    }
    
    # Check other important modules
    important_modules = [
        'numpy', 'pandas', 'rdkit', 'matplotlib', 'seaborn', 
        'psutil', 'sqlalchemy', 'marshmallow', 'flask_restful', 
        'flask_cors', 'reportlab', 'xlsxwriter', 'jwt', 'psycopg2'
    ]
    
    for module_name in important_modules:
        try:
            module = __import__(module_name)
            dependencies[module_name] = {
                'status': 'available',
                'version': get_module_version(module)
            }
        except ImportError as e:
            dependencies[module_name] = {
                'status': 'missing',
                'error': str(e)
            }
    
    # Special case for rdkit.Chem
    try:
        from rdkit import Chem
        dependencies['rdkit.Chem'] = {
            'status': 'available'
        }
    except ImportError as e:
        dependencies['rdkit.Chem'] = {
            'status': 'missing',
            'error': str(e)
        }
    
    return jsonify(dependencies)

@app.route('/api/sample/cryoprotectants', methods=['GET'])
def sample_cryoprotectants():
    """Return a sample list of cryoprotectants for demonstration"""
    # Sample data for testing
    sample_data = [
        {
            "id": 1,
            "name": "Glycerol",
            "molecular_formula": "C3H8O3",
            "molecular_weight": 92.09,
            "smiles": "C(C(CO)O)O",
            "inchi": "InChI=1S/C3H8O3/c4-1-3(6)2-5/h3-6H,1-2H2",
            "common_uses": ["Cryopreservation of cells", "Protein stabilization"]
        },
        {
            "id": 2,
            "name": "Dimethyl sulfoxide",
            "molecular_formula": "C2H6OS",
            "molecular_weight": 78.13,
            "smiles": "CS(=O)C",
            "inchi": "InChI=1S/C2H6OS/c1-4(2)3/h1-2H3",
            "common_uses": ["Cell freezing", "Solvent for pharmaceuticals"]
        },
        {
            "id": 3,
            "name": "Ethylene glycol",
            "molecular_formula": "C2H6O2",
            "molecular_weight": 62.07,
            "smiles": "OCCO",
            "inchi": "InChI=1S/C2H6O2/c3-1-2-4/h3-4H,1-2H2",
            "common_uses": ["Vitrification", "Antifreeze"]
        }
    ]
    
    return jsonify({
        "cryoprotectants": sample_data,
        "count": len(sample_data),
        "source": "sample data (not from database)"
    })

if __name__ == '__main__':
    port = int(os.environ.get('PORT', 5000))
    host = os.environ.get('HOST', '0.0.0.0')
    debug = os.environ.get('FLASK_ENV', 'development') == 'development'
    
    print(f"Starting simplified CryoProtect app on {host}:{port} (debug={debug})")
    print(f"Supabase URL: {'Configured' if SUPABASE_URL else 'Not configured'}")
    print(f"Supabase Key: {'Configured' if SUPABASE_KEY else 'Not configured'}")
    
    app.run(host=host, port=port, debug=debug)
