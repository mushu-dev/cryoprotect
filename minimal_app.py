import os
import logging
import sys
from flask import Flask, jsonify, request
import requests
import importlib

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler(sys.stdout)
    ]
)

logger = logging.getLogger(__name__)

# Create Flask app
app = Flask(__name__)

@app.route('/')
def index():
    """Root endpoint with basic information"""
    return jsonify({
        'name': 'CryoProtect Minimal App',
        'status': 'running',
        'environment': os.environ.get('FLASK_ENV', 'development'),
        'python_version': sys.version
    })

@app.route('/health')
def health():
    """Health check endpoint"""
    return jsonify({
        'status': 'healthy',
        'uptime': 'unknown'  # Would calculate from start time in a full app
    })

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
    supabase_url = os.environ.get('SUPABASE_URL')
    supabase_key = os.environ.get('SUPABASE_KEY')
    
    if not supabase_url or not supabase_key:
        return jsonify({
            'status': 'error',
            'message': 'Supabase URL or key not configured'
        }), 500
    
    try:
        # Test Supabase REST API
        headers = {
            'apikey': supabase_key,
            'Authorization': f'Bearer {supabase_key}'
        }
        
        response = requests.get(f"{supabase_url}/rest/v1/", headers=headers, timeout=5)
        
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

def get_module_version(module_name):
    """Get module version safely"""
    try:
        module = importlib.import_module(module_name)
        try:
            return getattr(module, '__version__')
        except AttributeError:
            try:
                return getattr(module, 'VERSION', 'Unknown')
            except AttributeError:
                try:
                    return getattr(module, 'version', 'Unknown')
                except AttributeError:
                    try: 
                        return getattr(module, 'Version', 'Unknown')
                    except AttributeError:
                        return "Available (version unknown)"
    except Exception as e:
        return f"Error: {str(e)}"
        
@app.route('/dependencies')
def check_dependencies():
    """Check if important dependencies are available"""
    dependencies = {}
    
    # List of dependencies to check
    deps_to_check = [
        'psutil',
        'pandas',
        'numpy',
        'flask',
        'rdkit',
        'matplotlib',
        'seaborn', 
        'reportlab',
        'xlsxwriter',
        'requests',
        'sqlalchemy',
        'psycopg2',
        'jwt',
        'pyjwt'
    ]
    
    for dep in deps_to_check:
        try:
            module = importlib.import_module(dep)
            dependencies[dep] = {
                'status': 'available',
                'version': get_module_version(dep)
            }
        except ImportError as e:
            dependencies[dep] = {
                'status': 'missing',
                'error': str(e)
            }
    
    # Special cases for complex imports
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
        
    try:
        import flask_restful
        dependencies['flask_restful'] = {
            'status': 'available',
            'version': get_module_version('flask_restful')
        }
    except ImportError as e:
        dependencies['flask_restful'] = {
            'status': 'missing',
            'error': str(e)
        }
        
    try:
        import flask_cors
        dependencies['flask_cors'] = {
            'status': 'available',
            'version': get_module_version('flask_cors')
        }
    except ImportError as e:
        dependencies['flask_cors'] = {
            'status': 'missing',
            'error': str(e)
        }
    
    try:
        import marshmallow
        dependencies['marshmallow'] = {
            'status': 'available',
            'version': get_module_version('marshmallow')
        }
    except ImportError as e:
        dependencies['marshmallow'] = {
            'status': 'missing',
            'error': str(e)
        }
        
    return jsonify(dependencies)

if __name__ == '__main__':
    port = int(os.environ.get('PORT', 5000))
    host = os.environ.get('HOST', '0.0.0.0')
    debug = os.environ.get('FLASK_ENV', 'development') == 'development'
    
    print(f"Starting minimal Flask app on {host}:{port} (debug={debug})")
    app.run(host=host, port=port, debug=debug)
