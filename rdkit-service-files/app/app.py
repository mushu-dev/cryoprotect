#!/usr/bin/env python3
"""
RDKit Microservice API - CryoProtect

This microservice provides RDKit functionality via a RESTful API.
It's designed to be deployed on Fly.io and connect to the main
CryoProtect application hosted on Heroku.
"""

import os
import logging
from flask import Flask, jsonify, request
from flask_cors import CORS
import time

# Import RDKit utilities (from current directory)
from app.rdkit_utils import (
    parse_molecule,
    calculate_all_properties,
    generate_molecule_svg,
    perform_substructure_search,
    calculate_similarity
)

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger('rdkit-service')

app = Flask(__name__)
CORS(app)  # Enable CORS for all routes

# Get API key from environment (if configured)
API_KEY = os.environ.get('API_KEY', '')
REQUIRE_API_KEY = bool(API_KEY)

def api_key_required(f):
    """Decorator to validate API key if configured."""
    def decorated_function(*args, **kwargs):
        if not REQUIRE_API_KEY:
            return f(*args, **kwargs)
            
        api_key = request.headers.get('X-API-Key')
        if api_key and api_key == API_KEY:
            return f(*args, **kwargs)
        return jsonify({"error": "Invalid or missing API key"}), 401
    
    decorated_function.__name__ = f.__name__
    return decorated_function

@app.route('/')
def index():
    """Root endpoint with service information."""
    return jsonify({
        'status': 'ok',
        'service': 'CryoProtect RDKit Microservice',
        'version': '1.0.0'
    })

@app.route('/health')
def health_check():
    """Health check endpoint for Fly.io monitoring."""
    # For faster health checks during deployment, don't load RDKit here
    # Just return a success response
    return jsonify({
        'status': 'ok',
        'service': 'CryoProtect RDKit Microservice',
        'timestamp': time.time()
    })

@app.route('/rdkit-check')
def rdkit_version_check():
    """Deep health check with RDKit loaded."""
    try:
        # Simple RDKit import check
        import rdkit
        version = rdkit.__version__
        
        return jsonify({
            'status': 'ok',
            'rdkit_version': version,
            'timestamp': time.time()
        })
    except Exception as e:
        logger.error(f"RDKit check failed: {str(e)}")
        return jsonify({
            'status': 'error',
            'message': str(e)
        }), 500

@app.route('/api/calculate-properties', methods=['POST'])
@api_key_required
def calculate_properties():
    """
    Calculate molecular properties using RDKit.
    
    Expected JSON payload:
    {
        "molecule_data": "CCO",
        "input_format": "smiles"
    }
    """
    try:
        data = request.json
        if not data:
            return jsonify({"error": "No JSON data provided"}), 400
            
        molecule_data = data.get('molecule_data')
        input_format = data.get('input_format', 'smiles')
        
        if not molecule_data:
            return jsonify({"error": "No molecule_data provided"}), 400
            
        # Calculate properties
        properties = calculate_all_properties(molecule_data, input_format)
        
        if "error" in properties:
            return jsonify({"error": properties["error"]}), 400
            
        return jsonify({
            "status": "success",
            "data": properties
        })
    except Exception as e:
        logger.error(f"Error calculating properties: {str(e)}")
        return jsonify({
            "status": "error",
            "message": str(e)
        }), 500

@app.route('/api/visualization', methods=['POST'])
@api_key_required
def generate_visualization():
    """
    Generate SVG visualization of a molecule.
    
    Expected JSON payload:
    {
        "molecule_data": "CCO",
        "input_format": "smiles",
        "width": 400,
        "height": 300,
        "highlight_atoms": [0, 1]
    }
    """
    try:
        data = request.json
        if not data:
            return jsonify({"error": "No JSON data provided"}), 400
            
        molecule_data = data.get('molecule_data')
        input_format = data.get('input_format', 'smiles')
        width = data.get('width', 400)
        height = data.get('height', 300)
        highlight_atoms = data.get('highlight_atoms')
        
        if not molecule_data:
            return jsonify({"error": "No molecule_data provided"}), 400
            
        # Generate SVG
        svg = generate_molecule_svg(
            molecule_data,
            input_format,
            width,
            height,
            highlight_atoms
        )
        
        if not svg:
            return jsonify({"error": "Failed to generate visualization"}), 400
            
        return jsonify({
            "status": "success",
            "data": {
                "svg": svg,
                "width": width,
                "height": height
            }
        })
    except Exception as e:
        logger.error(f"Error generating visualization: {str(e)}")
        return jsonify({
            "status": "error",
            "message": str(e)
        }), 500

@app.route('/api/substructure-search', methods=['POST'])
@api_key_required
def substructure_search():
    """
    Perform substructure search.
    
    Expected JSON payload:
    {
        "query_mol_data": "[OH]",
        "target_mol_data": "CCO",
        "query_format": "smarts",
        "target_format": "smiles"
    }
    """
    try:
        data = request.json
        if not data:
            return jsonify({"error": "No JSON data provided"}), 400
            
        query_mol_data = data.get('query_mol_data')
        target_mol_data = data.get('target_mol_data')
        query_format = data.get('query_format', 'smarts')
        target_format = data.get('target_format', 'smiles')
        
        if not query_mol_data or not target_mol_data:
            return jsonify({"error": "Missing molecule data"}), 400
            
        # Perform search
        result = perform_substructure_search(
            query_mol_data,
            target_mol_data,
            query_format,
            target_format
        )
        
        if "error" in result:
            return jsonify({"error": result["error"]}), 400
            
        return jsonify({
            "status": "success",
            "data": result
        })
    except Exception as e:
        logger.error(f"Error in substructure search: {str(e)}")
        return jsonify({
            "status": "error",
            "message": str(e)
        }), 500

@app.route('/api/similarity', methods=['POST'])
@api_key_required
def calculate_mol_similarity():
    """
    Calculate similarity between two molecules.
    
    Expected JSON payload:
    {
        "mol1_data": "CCO",
        "mol2_data": "CC(=O)O",
        "mol1_format": "smiles",
        "mol2_format": "smiles",
        "fingerprint_type": "morgan"
    }
    """
    try:
        data = request.json
        if not data:
            return jsonify({"error": "No JSON data provided"}), 400
            
        mol1_data = data.get('mol1_data')
        mol2_data = data.get('mol2_data')
        mol1_format = data.get('mol1_format', 'smiles')
        mol2_format = data.get('mol2_format', 'smiles')
        fingerprint_type = data.get('fingerprint_type', 'morgan')
        
        if not mol1_data or not mol2_data:
            return jsonify({"error": "Missing molecule data"}), 400
            
        # Calculate similarity
        result = calculate_similarity(
            mol1_data,
            mol2_data,
            mol1_format,
            mol2_format,
            fingerprint_type
        )
        
        if "error" in result:
            return jsonify({"error": result["error"]}), 400
            
        return jsonify({
            "status": "success",
            "data": result
        })
    except Exception as e:
        logger.error(f"Error calculating similarity: {str(e)}")
        return jsonify({
            "status": "error",
            "message": str(e)
        }), 500

if __name__ == '__main__':
    port = int(os.environ.get('PORT', 8080))
    app.run(host='0.0.0.0', port=port)