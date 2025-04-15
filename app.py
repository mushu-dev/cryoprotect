#!/usr/bin/env python3
"""
CryoProtect Analyzer API

This is the main Flask application for the CryoProtect Analyzer API.
It provides endpoints for accessing and manipulating data in the Supabase database.
"""

import os
import logging
from flask import Flask, jsonify, g, render_template
from flask_cors import CORS
from apispec import APISpec
from apispec.ext.marshmallow import MarshmallowPlugin
from flask_apispec.extension import FlaskApiSpec

from config import active_config
from api import init_app
from api.utils import get_supabase_client, authenticate_user

def create_app(config_object=None):
    """
    Create and configure the Flask application.
    
    Args:
        config_object: Configuration object to use
        
    Returns:
        Flask: Configured Flask application
    """
    app = Flask(__name__)
    
    # Load configuration
    if config_object:
        app.config.from_object(config_object)
    else:
        app.config.from_object(active_config)
    
    # Enable CORS
    CORS(app)
    
    # Configure logging
    if not app.debug:
        handler = logging.StreamHandler()
        handler.setLevel(logging.INFO)
        app.logger.addHandler(handler)
    
    # Initialize API
    init_app(app)
    
    # Configure API documentation
    app.config.update({
        'APISPEC_SPEC': APISpec(
            title=app.config['API_TITLE'],
            version=app.config['API_VERSION'],
            openapi_version=app.config['OPENAPI_VERSION'],
            plugins=[MarshmallowPlugin()],
        ),
        'APISPEC_SWAGGER_URL': '/swagger/',
        'APISPEC_SWAGGER_UI_URL': '/swagger-ui/'
    })
    docs = FlaskApiSpec(app)
    
    # Register error handlers
    @app.errorhandler(404)
    def not_found(error):
        return jsonify({'message': 'Resource not found'}), 404
    
    @app.errorhandler(500)
    def server_error(error):
        return jsonify({'message': 'Internal server error'}), 500
    
    # Register before request handler
    @app.before_request
    def before_request():
        # Initialize Supabase client
        get_supabase_client()
    
    # Register teardown request handler
    @app.teardown_request
    def teardown_request(exception=None):
        # Clean up resources
        if hasattr(g, 'supabase'):
            del g.supabase
    
    # Add health check endpoint
    @app.route('/health')
    def health_check():
        return jsonify({'status': 'ok', 'version': app.config['API_VERSION']}), 200
    
    # Add authentication endpoint
    @app.route('/auth/login', methods=['POST'])
    def login():
        user = authenticate_user()
        if user:
            return jsonify({
                'message': 'Authentication successful',
                'user': {
                    'id': user.id,
                    'email': user.email
                }
            }), 200
        else:
            return jsonify({'message': 'Authentication failed'}), 401
    
    return app

# Create the app instance
app = create_app()

# Add routes for the web interface
@app.route('/')
def index():
    return render_template('index.html')

@app.route('/molecules')
def molecules():
    return render_template('molecules.html')

@app.route('/molecules/rdkit')
def molecules_rdkit():
    return render_template('molecules_rdkit.html')

@app.route('/mixtures')
def mixtures():
    return render_template('mixtures.html')

@app.route('/predictions')
def predictions():
    return render_template('predictions.html')

@app.route('/experiments')
def experiments():
    return render_template('experiments.html')

@app.route('/comparisons')
def comparisons():
    return render_template('comparisons.html')

@app.route('/login')
def login_page():
    return render_template('login.html')

if __name__ == '__main__':
    # Try to authenticate with Supabase
    with app.app_context():
        user = authenticate_user()
        if user:
            app.logger.info(f"Authenticated as {user.email}")
        else:
            app.logger.warning("No authentication. Some operations may fail due to Row Level Security (RLS) policies.")
    
    # Run the application
    port = int(os.environ.get('PORT', 5000))
    app.run(host='0.0.0.0', port=port)