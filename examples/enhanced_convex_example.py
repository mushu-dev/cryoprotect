#!/usr/bin/env python3
"""
Enhanced Convex Adapter Example

This script demonstrates how to use the enhanced Convex adapter for
integrating Convex with a Flask application.

Usage:
    python enhanced_convex_example.py

Requirements:
    - Flask
    - dotenv
    - ConvexAdapter
    - AuthBridge
"""

import os
import sys
import logging
from flask import Flask, jsonify, request, g
from dotenv import load_dotenv

# Add parent directory to path for imports
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from database.enhanced_convex_adapter import ConvexAdapter
from database.auth_bridge import init_auth_bridge
from database.flask_convex_integration import init_flask_convex

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Load environment variables
load_dotenv()

# Create Flask app
app = Flask(__name__)

# Initialize Flask-Convex integration
flask_convex = init_flask_convex(app, {
    'convex_url': os.environ.get('CONVEX_URL', ''),
    'convex_key': os.environ.get('CONVEX_DEPLOYMENT_KEY', ''),
    'jwt_secret': os.environ.get('JWT_SECRET', 'your-jwt-secret-key')
})

# Get Convex adapter
convex = flask_convex.get_convex_adapter()

# Example routes

@app.route('/api/molecules', methods=['GET'])
def get_molecules():
    """Get all molecules."""
    try:
        # Get molecules from Convex
        molecules = convex.execute_query('api.molecules.list', {
            'limit': int(request.args.get('limit', 100))
        })
        
        return jsonify({
            'success': True,
            'data': molecules
        })
        
    except Exception as e:
        logger.error(f"Error getting molecules: {str(e)}")
        return jsonify({
            'success': False,
            'error': str(e)
        }), 500

@app.route('/api/molecules/<molecule_id>', methods=['GET'])
def get_molecule(molecule_id):
    """Get a specific molecule."""
    try:
        # Get molecule from Convex
        molecule = convex.execute_query('api.molecules.get', {
            'id': molecule_id
        })
        
        if not molecule:
            return jsonify({
                'success': False,
                'error': 'Molecule not found'
            }), 404
        
        return jsonify({
            'success': True,
            'data': molecule
        })
        
    except Exception as e:
        logger.error(f"Error getting molecule {molecule_id}: {str(e)}")
        return jsonify({
            'success': False,
            'error': str(e)
        }), 500

@app.route('/api/molecules', methods=['POST'])
def create_molecule():
    """Create a new molecule."""
    try:
        # Get request data
        data = request.json
        
        # Validate request data
        if not data or not data.get('name') or not data.get('formula'):
            return jsonify({
                'success': False,
                'error': 'Invalid request data'
            }), 400
        
        # Demonstrate transaction usage
        with flask_convex.transaction() as tx:
            # Create molecule in Convex
            molecule_id = convex.execute_query('api.molecules.create', {
                'data': data
            })
        
        return jsonify({
            'success': True,
            'data': {
                'id': molecule_id
            }
        }), 201
        
    except Exception as e:
        logger.error(f"Error creating molecule: {str(e)}")
        return jsonify({
            'success': False,
            'error': str(e)
        }), 500

@app.route('/api/molecules/<molecule_id>', methods=['PUT'])
def update_molecule(molecule_id):
    """Update a molecule."""
    try:
        # Get request data
        data = request.json
        
        # Validate request data
        if not data:
            return jsonify({
                'success': False,
                'error': 'Invalid request data'
            }), 400
        
        # Get current molecule to ensure it exists
        current = convex.execute_query('api.molecules.get', {
            'id': molecule_id
        })
        
        if not current:
            return jsonify({
                'success': False,
                'error': 'Molecule not found'
            }), 404
        
        # Update molecule in Convex
        success = convex.execute_query('api.molecules.update', {
            'id': molecule_id,
            'data': data
        })
        
        return jsonify({
            'success': success,
            'data': {
                'id': molecule_id
            }
        })
        
    except Exception as e:
        logger.error(f"Error updating molecule {molecule_id}: {str(e)}")
        return jsonify({
            'success': False,
            'error': str(e)
        }), 500

@app.route('/api/molecules/<molecule_id>', methods=['DELETE'])
def delete_molecule(molecule_id):
    """Delete a molecule."""
    try:
        # Delete molecule in Convex
        success = convex.execute_query('api.molecules.delete', {
            'id': molecule_id
        })
        
        return jsonify({
            'success': success
        })
        
    except Exception as e:
        logger.error(f"Error deleting molecule {molecule_id}: {str(e)}")
        return jsonify({
            'success': False,
            'error': str(e)
        }), 500

@app.route('/api/health', methods=['GET'])
def health():
    """Health check."""
    try:
        # Get connection info
        connection_info = convex.get_connection_info()
        
        return jsonify({
            'success': True,
            'status': 'healthy',
            'convex': connection_info
        })
        
    except Exception as e:
        logger.error(f"Error in health check: {str(e)}")
        return jsonify({
            'success': False,
            'status': 'unhealthy',
            'error': str(e)
        }), 500

@app.route('/api/login', methods=['POST'])
def login():
    """Example login endpoint."""
    try:
        # Get credentials
        data = request.json
        email = data.get('email')
        password = data.get('password')
        
        if not email or not password:
            return jsonify({
                'success': False,
                'error': 'Invalid credentials'
            }), 400
        
        # In a real app, you would verify credentials against a database
        # For this example, we'll create a token for any valid-looking email
        if '@' not in email:
            return jsonify({
                'success': False,
                'error': 'Invalid email format'
            }), 400
        
        # Create a user ID (in a real app, this would come from your database)
        user_id = 'user_' + email.split('@')[0]
        
        # Create user data
        user_data = {
            'email': email,
            'name': email.split('@')[0].title(),
            'role': 'user'
        }
        
        # Generate token using AuthBridge
        auth_bridge = flask_convex.get_auth_bridge()
        token = auth_bridge.generate_token(user_id, user_data)
        
        # Create Convex identity for frontend use
        convex_identity = auth_bridge.create_convex_identity(user_id, user_data)
        
        return jsonify({
            'success': True,
            'token': token,
            'user': user_data,
            'convex_identity': convex_identity
        })
        
    except Exception as e:
        logger.error(f"Error in login: {str(e)}")
        return jsonify({
            'success': False,
            'error': str(e)
        }), 500

# Run the app
if __name__ == '__main__':
    # Check if Convex adapter is connected
    if not convex.connected and not convex.connect():
        logger.error("Failed to connect to Convex")
        sys.exit(1)
    
    logger.info("Connected to Convex")
    logger.info("Started Flask app")
    
    app.run(host='0.0.0.0', port=5000, debug=True)