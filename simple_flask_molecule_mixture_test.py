#!/usr/bin/env python3
"""
Simple Flask application for testing the Molecule and Mixture API integration.

This script creates a minimal Flask application that registers only the
necessary routes and models for testing the molecule and mixture API integration.
"""

from flask import Flask, request, jsonify, Blueprint
from flask_cors import CORS
import logging
import uuid
from datetime import datetime
import json
import os

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)

logger = logging.getLogger(__name__)

# Create Flask application
app = Flask(__name__)
CORS(app)  # Enable CORS for all routes

# Simple in-memory database for testing
db = {
    'molecules': [
        {
            'id': str(uuid.uuid4()),
            'name': 'Glycerol',
            'smiles': 'C(C(CO)O)O',
            'pubchem_cid': '753',
            'formula': 'C3H8O3',
            'molecular_weight': 92.09,
            'is_cryoprotectant': True,
            'created_at': datetime.now().isoformat(),
            'updated_at': datetime.now().isoformat()
        },
        {
            'id': str(uuid.uuid4()),
            'name': 'Dimethyl sulfoxide',
            'smiles': 'CS(=O)C',
            'pubchem_cid': '679',
            'formula': 'C2H6OS',
            'molecular_weight': 78.13,
            'is_cryoprotectant': True,
            'created_at': datetime.now().isoformat(),
            'updated_at': datetime.now().isoformat()
        },
        {
            'id': str(uuid.uuid4()),
            'name': 'Ethylene glycol',
            'smiles': 'C(CO)O',
            'pubchem_cid': '174',
            'formula': 'C2H6O2',
            'molecular_weight': 62.07,
            'is_cryoprotectant': True,
            'created_at': datetime.now().isoformat(),
            'updated_at': datetime.now().isoformat()
        }
    ],
    'properties': [],
    'mixtures': [],
    'mixture_components': []
}

# Initialize blueprints
molecule_bp = Blueprint('molecules', __name__)
mixture_bp = Blueprint('mixtures', __name__)

# Helper function to generate molecule properties
def generate_molecule_properties():
    """Generate properties for test molecules"""
    properties = []
    
    # Define some property types
    property_types = [
        {'name': 'melting_point', 'unit': 'C'},
        {'name': 'boiling_point', 'unit': 'C'},
        {'name': 'density', 'unit': 'g/mL'},
        {'name': 'log_p', 'unit': ''},
        {'name': 'solubility', 'unit': 'mg/mL'}
    ]
    
    # Generate values for glycerol
    properties.append({
        'id': str(uuid.uuid4()),
        'molecule_id': db['molecules'][0]['id'],
        'property_name': 'melting_point',
        'property_value': 18.2,
        'property_unit': 'C',
        'calculation_method': 'experimental',
        'created_at': datetime.now().isoformat(),
        'updated_at': datetime.now().isoformat()
    })
    
    properties.append({
        'id': str(uuid.uuid4()),
        'molecule_id': db['molecules'][0]['id'],
        'property_name': 'boiling_point',
        'property_value': 290.0,
        'property_unit': 'C',
        'calculation_method': 'experimental',
        'created_at': datetime.now().isoformat(),
        'updated_at': datetime.now().isoformat()
    })
    
    # Generate values for DMSO
    properties.append({
        'id': str(uuid.uuid4()),
        'molecule_id': db['molecules'][1]['id'],
        'property_name': 'melting_point',
        'property_value': 19.0,
        'property_unit': 'C',
        'calculation_method': 'experimental',
        'created_at': datetime.now().isoformat(),
        'updated_at': datetime.now().isoformat()
    })
    
    properties.append({
        'id': str(uuid.uuid4()),
        'molecule_id': db['molecules'][1]['id'],
        'property_name': 'boiling_point',
        'property_value': 189.0,
        'property_unit': 'C',
        'calculation_method': 'experimental',
        'created_at': datetime.now().isoformat(),
        'updated_at': datetime.now().isoformat()
    })
    
    # Generate values for ethylene glycol
    properties.append({
        'id': str(uuid.uuid4()),
        'molecule_id': db['molecules'][2]['id'],
        'property_name': 'melting_point',
        'property_value': -12.9,
        'property_unit': 'C',
        'calculation_method': 'experimental',
        'created_at': datetime.now().isoformat(),
        'updated_at': datetime.now().isoformat()
    })
    
    properties.append({
        'id': str(uuid.uuid4()),
        'molecule_id': db['molecules'][2]['id'],
        'property_name': 'boiling_point',
        'property_value': 197.3,
        'property_unit': 'C',
        'calculation_method': 'experimental',
        'created_at': datetime.now().isoformat(),
        'updated_at': datetime.now().isoformat()
    })
    
    return properties

# Generate some test mixtures
def generate_test_mixtures():
    """Generate test mixtures"""
    mixtures = []
    components = []
    
    # Create a DMSO-Glycerol mixture
    mixture_id = str(uuid.uuid4())
    mixtures.append({
        'id': mixture_id,
        'name': 'DMSO-Glycerol Mixture',
        'description': 'A test mixture of DMSO and Glycerol',
        'is_cryoprotectant_mixture': True,
        'user_id': 'test-user',
        'created_at': datetime.now().isoformat(),
        'updated_at': datetime.now().isoformat(),
        'total_concentration': 15.0,
        'cryoprotection_score': 85.3
    })
    
    # Add DMSO component
    components.append({
        'id': str(uuid.uuid4()),
        'mixture_id': mixture_id,
        'molecule_id': db['molecules'][1]['id'],
        'concentration': 10.0,
        'concentration_unit': 'mM',
        'created_at': datetime.now().isoformat(),
        'updated_at': datetime.now().isoformat()
    })
    
    # Add glycerol component
    components.append({
        'id': str(uuid.uuid4()),
        'mixture_id': mixture_id,
        'molecule_id': db['molecules'][0]['id'],
        'concentration': 5.0,
        'concentration_unit': 'mM',
        'created_at': datetime.now().isoformat(),
        'updated_at': datetime.now().isoformat()
    })
    
    # Create another mixture
    mixture_id = str(uuid.uuid4())
    mixtures.append({
        'id': mixture_id,
        'name': 'Triple Cryoprotectant Mix',
        'description': 'A test mixture of all three cryoprotectants',
        'is_cryoprotectant_mixture': True,
        'user_id': 'test-user',
        'created_at': datetime.now().isoformat(),
        'updated_at': datetime.now().isoformat(),
        'total_concentration': 30.0,
        'cryoprotection_score': 92.7
    })
    
    # Add components for the triple mix
    for i, molecule in enumerate(db['molecules']):
        components.append({
            'id': str(uuid.uuid4()),
            'mixture_id': mixture_id,
            'molecule_id': molecule['id'],
            'concentration': 10.0,
            'concentration_unit': 'mM',
            'created_at': datetime.now().isoformat(),
            'updated_at': datetime.now().isoformat()
        })
    
    return mixtures, components

# Helper function to generate provenance info
def generate_provenance(request):
    """Generate provenance information for database records"""
    return {
        'created_by': getattr(request, 'user_id', 'system'),
        'created_at': datetime.now().isoformat(),
        'source': 'api',
        'method': request.method,
        'version': '1.0'
    }

# Initialize the database with test data
def initialize_database():
    """Initialize the database with test data"""
    db['properties'] = generate_molecule_properties()
    db['mixtures'], db['mixture_components'] = generate_test_mixtures()

# Routes for molecules
@molecule_bp.route('', methods=['GET'])
def list_molecules():
    """List all molecules with pagination and filtering"""
    # Get query parameters
    page = int(request.args.get('page', 1))
    limit = int(request.args.get('limit', 10))
    search = request.args.get('search', '')
    sort_by = request.args.get('sort_by', 'name')
    sort_order = request.args.get('sort_order', 'asc')
    is_cryoprotectant = request.args.get('is_cryoprotectant')
    
    # Filter molecules
    filtered_molecules = db['molecules']
    
    if search:
        search = search.lower()
        filtered_molecules = [m for m in filtered_molecules if search in m['name'].lower() or search in m.get('smiles', '').lower()]
    
    if is_cryoprotectant is not None:
        is_cryo = is_cryoprotectant.lower() == 'true'
        filtered_molecules = [m for m in filtered_molecules if m.get('is_cryoprotectant') == is_cryo]
    
    # Apply sorting
    reverse = sort_order.lower() == 'desc'
    filtered_molecules = sorted(filtered_molecules, key=lambda m: m.get(sort_by, ''), reverse=reverse)
    
    # Apply pagination
    start_idx = (page - 1) * limit
    end_idx = start_idx + limit
    paginated_molecules = filtered_molecules[start_idx:end_idx]
    
    return jsonify({
        'molecules': paginated_molecules,
        'total': len(filtered_molecules),
        'page': page,
        'limit': limit
    })

@molecule_bp.route('/<molecule_id>', methods=['GET'])
def get_molecule(molecule_id):
    """Get molecule details by ID"""
    # Find molecule by ID
    molecule = next((m for m in db['molecules'] if m['id'] == molecule_id), None)
    
    if not molecule:
        return jsonify({'message': 'Molecule not found'}), 404
    
    return jsonify(molecule)

@molecule_bp.route('/<molecule_id>/properties', methods=['GET'])
def get_molecule_properties(molecule_id):
    """Get properties for a molecule"""
    # Find molecule by ID
    molecule = next((m for m in db['molecules'] if m['id'] == molecule_id), None)
    
    if not molecule:
        return jsonify({'message': 'Molecule not found'}), 404
    
    # Get properties for this molecule
    properties = [p for p in db['properties'] if p['molecule_id'] == molecule_id]
    
    return jsonify({
        'properties': properties
    })

@molecule_bp.route('/search', methods=['GET'])
def search_molecules():
    """Search molecules by name or SMILES"""
    query = request.args.get('query', '').lower()
    limit = int(request.args.get('limit', 10))
    
    # Search by name or SMILES
    results = [
        m for m in db['molecules'] 
        if query in m['name'].lower() or query in m.get('smiles', '').lower()
    ][:limit]
    
    return jsonify({
        'results': results
    })

@molecule_bp.route('/import/pubchem', methods=['POST'])
def import_from_pubchem():
    """Import a molecule from PubChem by CID"""
    data = request.get_json()
    cid = data.get('cid')
    
    if not cid:
        return jsonify({'message': 'PubChem CID is required'}), 400
    
    # Check if we already have this molecule
    existing = next((m for m in db['molecules'] if m.get('pubchem_cid') == cid), None)
    if existing:
        return jsonify(existing)
    
    # Create a new molecule with mock data
    molecule = {
        'id': str(uuid.uuid4()),
        'name': f'Imported Molecule {cid}',
        'smiles': 'CC(=O)OC1=CC=CC=C1C(=O)O',  # Mock SMILES
        'pubchem_cid': cid,
        'formula': 'C9H8O4',  # Mock formula
        'molecular_weight': 180.16,  # Mock weight
        'is_cryoprotectant': False,
        'created_at': datetime.now().isoformat(),
        'updated_at': datetime.now().isoformat()
    }
    
    # Add to database
    db['molecules'].append(molecule)
    
    return jsonify(molecule)

# Routes for mixtures
@mixture_bp.route('', methods=['GET'])
def list_mixtures():
    """List all mixtures with pagination and filtering"""
    # Get query parameters
    page = int(request.args.get('page', 1))
    limit = int(request.args.get('limit', 10))
    search = request.args.get('search', '')
    sort_by = request.args.get('sort_by', 'name')
    sort_order = request.args.get('sort_order', 'asc')
    is_cryoprotectant = request.args.get('is_cryoprotectant_mixture')
    user_id = request.args.get('user_id')
    
    # Filter mixtures
    filtered_mixtures = db['mixtures']
    
    if search:
        search = search.lower()
        filtered_mixtures = [m for m in filtered_mixtures if search in m['name'].lower() or (m.get('description') and search in m['description'].lower())]
    
    if is_cryoprotectant is not None:
        is_cryo = is_cryoprotectant.lower() == 'true'
        filtered_mixtures = [m for m in filtered_mixtures if m.get('is_cryoprotectant_mixture') == is_cryo]
    
    if user_id:
        filtered_mixtures = [m for m in filtered_mixtures if m.get('user_id') == user_id]
    
    # Apply sorting
    reverse = sort_order.lower() == 'desc'
    filtered_mixtures = sorted(filtered_mixtures, key=lambda m: m.get(sort_by, ''), reverse=reverse)
    
    # Apply pagination
    start_idx = (page - 1) * limit
    end_idx = start_idx + limit
    paginated_mixtures = filtered_mixtures[start_idx:end_idx]
    
    return jsonify({
        'mixtures': paginated_mixtures,
        'total': len(filtered_mixtures),
        'page': page,
        'limit': limit
    })

@mixture_bp.route('/<mixture_id>', methods=['GET'])
def get_mixture(mixture_id):
    """Get mixture details by ID"""
    # Find mixture by ID
    mixture = next((m for m in db['mixtures'] if m['id'] == mixture_id), None)
    
    if not mixture:
        return jsonify({'message': 'Mixture not found'}), 404
    
    # Add components
    components = [c for c in db['mixture_components'] if c['mixture_id'] == mixture_id]
    
    # Enhance components with molecule info
    for component in components:
        molecule_id = component['molecule_id']
        molecule = next((m for m in db['molecules'] if m['id'] == molecule_id), None)
        if molecule:
            component['molecule'] = molecule
    
    mixture = dict(mixture)  # Create a copy
    mixture['components'] = components
    
    return jsonify(mixture)

@mixture_bp.route('', methods=['POST'])
def create_mixture():
    """Create a new mixture"""
    data = request.get_json()
    
    # Validate input
    if 'name' not in data:
        return jsonify({'message': 'Mixture name is required'}), 400
    
    if 'components' not in data or not isinstance(data['components'], list) or len(data['components']) == 0:
        return jsonify({'message': 'Mixture must have at least one component'}), 400
    
    # Create mixture
    mixture_id = str(uuid.uuid4())
    mixture = {
        'id': mixture_id,
        'name': data['name'],
        'description': data.get('description', ''),
        'is_cryoprotectant_mixture': True,  # Default to true for test
        'user_id': 'test-user',  # Would be real user ID in production
        'created_at': datetime.now().isoformat(),
        'updated_at': datetime.now().isoformat(),
        'total_concentration': sum(comp.get('concentration', 0) for comp in data['components']),
        'cryoprotection_score': 75.0  # Dummy score for testing
    }
    
    # Create components
    components = []
    for comp_data in data['components']:
        # Validate component
        if 'molecule_id' not in comp_data:
            return jsonify({'message': 'Each component must have a molecule_id'}), 400
        
        if 'concentration' not in comp_data:
            return jsonify({'message': 'Each component must have a concentration'}), 400
        
        if 'concentration_unit' not in comp_data:
            return jsonify({'message': 'Each component must have a concentration_unit'}), 400
        
        # Create component
        component = {
            'id': str(uuid.uuid4()),
            'mixture_id': mixture_id,
            'molecule_id': comp_data['molecule_id'],
            'concentration': comp_data['concentration'],
            'concentration_unit': comp_data['concentration_unit'],
            'created_at': datetime.now().isoformat(),
            'updated_at': datetime.now().isoformat()
        }
        
        components.append(component)
    
    # Add to database
    db['mixtures'].append(mixture)
    db['mixture_components'].extend(components)
    
    # Return full mixture with components
    return get_mixture(mixture_id)

@mixture_bp.route('/<mixture_id>', methods=['PATCH'])
def update_mixture(mixture_id):
    """Update a mixture"""
    data = request.get_json()
    
    # Find mixture by ID
    mixture = next((m for m in db['mixtures'] if m['id'] == mixture_id), None)
    
    if not mixture:
        return jsonify({'message': 'Mixture not found'}), 404
    
    # Update fields
    if 'name' in data:
        mixture['name'] = data['name']
    
    if 'description' in data:
        mixture['description'] = data['description']
    
    if 'is_cryoprotectant_mixture' in data:
        mixture['is_cryoprotectant_mixture'] = data['is_cryoprotectant_mixture']
    
    mixture['updated_at'] = datetime.now().isoformat()
    
    return jsonify(mixture)

@mixture_bp.route('/<mixture_id>', methods=['DELETE'])
def delete_mixture(mixture_id):
    """Delete a mixture"""
    # Find mixture by ID
    mixture_idx = next((i for i, m in enumerate(db['mixtures']) if m['id'] == mixture_id), None)
    
    if mixture_idx is None:
        return jsonify({'message': 'Mixture not found'}), 404
    
    # Remove from database
    db['mixtures'].pop(mixture_idx)
    
    # Remove components
    db['mixture_components'] = [c for c in db['mixture_components'] if c['mixture_id'] != mixture_id]
    
    return jsonify({'message': 'Mixture deleted successfully'}), 200

@mixture_bp.route('/<mixture_id>/components', methods=['POST'])
def add_component(mixture_id):
    """Add a component to a mixture"""
    data = request.get_json()
    
    # Find mixture by ID
    mixture = next((m for m in db['mixtures'] if m['id'] == mixture_id), None)
    
    if not mixture:
        return jsonify({'message': 'Mixture not found'}), 404
    
    # Validate input
    if 'molecule_id' not in data:
        return jsonify({'message': 'molecule_id is required'}), 400
    
    if 'concentration' not in data:
        return jsonify({'message': 'concentration is required'}), 400
    
    if 'concentration_unit' not in data:
        return jsonify({'message': 'concentration_unit is required'}), 400
    
    # Create component
    component = {
        'id': str(uuid.uuid4()),
        'mixture_id': mixture_id,
        'molecule_id': data['molecule_id'],
        'concentration': data['concentration'],
        'concentration_unit': data['concentration_unit'],
        'created_at': datetime.now().isoformat(),
        'updated_at': datetime.now().isoformat()
    }
    
    # Add to database
    db['mixture_components'].append(component)
    
    # Update mixture total concentration
    mixture['total_concentration'] = sum(c['concentration'] for c in db['mixture_components'] if c['mixture_id'] == mixture_id)
    mixture['updated_at'] = datetime.now().isoformat()
    
    # Add molecule info
    molecule = next((m for m in db['molecules'] if m['id'] == data['molecule_id']), None)
    if molecule:
        component['molecule'] = molecule
    
    return jsonify(component), 201

@mixture_bp.route('/<mixture_id>/components/<component_id>', methods=['PATCH'])
def update_component(mixture_id, component_id):
    """Update a mixture component"""
    data = request.get_json()
    
    # Find component by ID
    component = next((c for c in db['mixture_components'] if c['id'] == component_id and c['mixture_id'] == mixture_id), None)
    
    if not component:
        return jsonify({'message': 'Component not found'}), 404
    
    # Update fields
    if 'concentration' in data:
        component['concentration'] = data['concentration']
    
    if 'concentration_unit' in data:
        component['concentration_unit'] = data['concentration_unit']
    
    component['updated_at'] = datetime.now().isoformat()
    
    # Update mixture total concentration
    mixture = next((m for m in db['mixtures'] if m['id'] == mixture_id), None)
    if mixture:
        mixture['total_concentration'] = sum(c['concentration'] for c in db['mixture_components'] if c['mixture_id'] == mixture_id)
        mixture['updated_at'] = datetime.now().isoformat()
    
    # Add molecule info
    molecule = next((m for m in db['molecules'] if m['id'] == component['molecule_id']), None)
    if molecule:
        component['molecule'] = molecule
    
    return jsonify(component)

@mixture_bp.route('/<mixture_id>/components/<component_id>', methods=['DELETE'])
def remove_component(mixture_id, component_id):
    """Remove a component from a mixture"""
    # Find component by ID
    component_idx = next((i for i, c in enumerate(db['mixture_components']) 
                         if c['id'] == component_id and c['mixture_id'] == mixture_id), None)
    
    if component_idx is None:
        return jsonify({'message': 'Component not found'}), 404
    
    # Remove from database
    db['mixture_components'].pop(component_idx)
    
    # Update mixture total concentration
    mixture = next((m for m in db['mixtures'] if m['id'] == mixture_id), None)
    if mixture:
        mixture['total_concentration'] = sum(c['concentration'] for c in db['mixture_components'] if c['mixture_id'] == mixture_id)
        mixture['updated_at'] = datetime.now().isoformat()
    
    return jsonify({'message': 'Component removed successfully'}), 200

@mixture_bp.route('/<mixture_id>/cryoprotection-score', methods=['GET'])
def get_cryoprotection_score(mixture_id):
    """Get cryoprotection score for a mixture"""
    # Find mixture by ID
    mixture = next((m for m in db['mixtures'] if m['id'] == mixture_id), None)
    
    if not mixture:
        return jsonify({'message': 'Mixture not found'}), 404
    
    # Return a placeholder score
    if 'cryoprotection_score' not in mixture:
        mixture['cryoprotection_score'] = 80.0
    
    return jsonify({
        'score': mixture['cryoprotection_score'],
        'details': {
            'method': 'simulation',
            'components': [
                {
                    'molecule_name': next((m['name'] for m in db['molecules'] if m['id'] == c['molecule_id']), 'Unknown'),
                    'contribution': c['concentration'] / mixture['total_concentration'] * 100 if mixture['total_concentration'] > 0 else 0
                }
                for c in db['mixture_components'] if c['mixture_id'] == mixture_id
            ]
        }
    })

@mixture_bp.route('/search', methods=['GET'])
def search_mixtures():
    """Search mixtures by name"""
    query = request.args.get('query', '').lower()
    limit = int(request.args.get('limit', 10))
    
    # Search by name or description
    results = [
        m for m in db['mixtures'] 
        if query in m['name'].lower() or 
        (m.get('description') and query in m['description'].lower())
    ][:limit]
    
    return jsonify({
        'results': results
    })

# Register the blueprints with URL prefixes
app.register_blueprint(molecule_bp, url_prefix='/api/v1/molecules')
app.register_blueprint(mixture_bp, url_prefix='/api/v1/mixtures')

# Home route
@app.route('/')
def home():
    """Home route"""
    return jsonify({
        'message': 'Molecule and Mixture API Test Server',
        'endpoints': {
            'GET /api/v1/molecules': 'List all molecules',
            'GET /api/v1/molecules/<id>': 'Get molecule details',
            'GET /api/v1/molecules/<id>/properties': 'Get molecule properties',
            'GET /api/v1/molecules/search': 'Search molecules',
            'POST /api/v1/molecules/import/pubchem': 'Import a molecule from PubChem',
            'GET /api/v1/mixtures': 'List all mixtures',
            'GET /api/v1/mixtures/<id>': 'Get mixture details',
            'POST /api/v1/mixtures': 'Create a new mixture',
            'PATCH /api/v1/mixtures/<id>': 'Update a mixture',
            'DELETE /api/v1/mixtures/<id>': 'Delete a mixture',
            'POST /api/v1/mixtures/<id>/components': 'Add a component to a mixture',
            'PATCH /api/v1/mixtures/<id>/components/<id>': 'Update a mixture component',
            'DELETE /api/v1/mixtures/<id>/components/<id>': 'Remove a component from a mixture',
            'GET /api/v1/mixtures/<id>/cryoprotection-score': 'Get cryoprotection score for a mixture',
            'GET /api/v1/mixtures/search': 'Search mixtures'
        }
    })

# Health check endpoint
@app.route('/api/v1/health')
def health_check():
    """Health check endpoint"""
    return jsonify({
        'status': 'ok',
        'timestamp': datetime.now().isoformat()
    })

# Helper endpoint to get test IDs
@app.route('/api/v1/test-ids')
def test_ids():
    """Get test molecule and mixture IDs"""
    if not db['molecules'] or not db['mixtures']:
        initialize_database()
    
    return jsonify({
        'molecule_ids': [m['id'] for m in db['molecules']],
        'mixture_ids': [m['id'] for m in db['mixtures']]
    })

if __name__ == '__main__':
    # Initialize database
    initialize_database()
    
    # Create data directory if it doesn't exist
    os.makedirs('data', exist_ok=True)
    
    # Run the app
    app.run(debug=True, port=5001, host='0.0.0.0')