#!/usr/bin/env python3
"""
Simplified CryoProtect API for Heroku deployment without RDKit dependencies.
"""

import os
import logging
from flask import Flask, jsonify, request
from flask_cors import CORS
import psycopg2
from urllib.parse import urlparse

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger('simple_app')

app = Flask(__name__)

# Enable CORS with more specific configuration
CORS(app, resources={
    r"/*": {"origins": [
        "http://localhost:3000",                        # Local development
        "https://frontend-cryoprotect.vercel.app",      # Vercel production
        "https://www.cryoprotect.app",                  # Custom domain
        os.environ.get("VERCEL_FRONTEND_URL", "*")      # Dynamic frontend URL from env
    ], "methods": ["GET", "POST", "PUT", "DELETE", "OPTIONS"], 
       "allow_headers": ["Content-Type", "Authorization"]}
})

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
        'message': 'CryoProtect API is running in simplified mode',
        'version': '1.0.0'
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