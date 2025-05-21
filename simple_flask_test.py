#!/usr/bin/env python3
"""
Simple Flask application for testing the Experiment API integration.

This script creates a minimal Flask application that registers only the
necessary routes and models for testing the experiment API integration.
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
    'experiments': [],
    'results': [],
    'time_series': [],
    'protocols': [
        {
            'id': str(uuid.uuid4()),
            'name': 'Test Protocol',
            'description': 'A test protocol for experimentation',
            'version': '1.0',
            'steps': [],
            'created_at': datetime.now().isoformat(),
            'updated_at': datetime.now().isoformat(),
            'created_by': 'test-user',
            'is_template': False,
            'provenance': {
                'created_by': 'system',
                'created_at': datetime.now().isoformat(),
                'source': 'system',
                'method': 'initialization',
                'version': '1.0'
            }
        }
    ],
    'tissue_types': [
        {
            'id': str(uuid.uuid4()),
            'name': 'Test Tissue',
            'description': 'A test tissue type for experimentation',
            'category': 'test',
            'source': 'test',
            'properties': {},
            'created_at': datetime.now().isoformat(),
            'updated_at': datetime.now().isoformat()
        }
    ]
}

# Initialize blueprint
experiment_bp = Blueprint('experiments', __name__)

# Helper function to get test protocol ID
def get_test_protocol_id():
    """Get the ID of the test protocol"""
    return db['protocols'][0]['id']

# Helper function to get test tissue type ID
def get_test_tissue_type_id():
    """Get the ID of the test tissue type"""
    return db['tissue_types'][0]['id']

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

# Routes for experiments
@experiment_bp.route('', methods=['GET'])
def list_experiments():
    """List all experiments with pagination and filtering"""
    # Get query parameters
    page = int(request.args.get('page', 1))
    per_page = int(request.args.get('per_page', 10))
    sort_by = request.args.get('sort_by', 'start_date')
    sort_order = request.args.get('sort_order', 'desc')
    status = request.args.get('status')
    experiment_type = request.args.get('experiment_type')
    researcher = request.args.get('researcher')
    tissue_type_id = request.args.get('tissue_type_id')
    
    # Filter experiments
    filtered_experiments = db['experiments']
    
    if status:
        filtered_experiments = [e for e in filtered_experiments if e['status'] == status]
    
    if experiment_type:
        filtered_experiments = [e for e in filtered_experiments if e['experiment_type'] == experiment_type]
    
    if researcher:
        filtered_experiments = [e for e in filtered_experiments if researcher.lower() in e['researcher'].lower()]
    
    if tissue_type_id:
        filtered_experiments = [e for e in filtered_experiments if e['tissue_type_id'] == tissue_type_id]
    
    # Apply sorting
    reverse = sort_order.lower() == 'desc'
    filtered_experiments = sorted(filtered_experiments, key=lambda e: e.get(sort_by, ''), reverse=reverse)
    
    # Apply pagination
    start_idx = (page - 1) * per_page
    end_idx = start_idx + per_page
    paginated_experiments = filtered_experiments[start_idx:end_idx]
    
    return jsonify({
        'data': paginated_experiments,
        'total': len(filtered_experiments),
        'page': page,
        'per_page': per_page,
        'total_pages': (len(filtered_experiments) + per_page - 1) // per_page
    })

@experiment_bp.route('', methods=['POST'])
def create_experiment():
    """Create a new experiment"""
    data = request.get_json()
    
    # Create new experiment with generated ID
    experiment_id = str(uuid.uuid4())
    experiment = {
        'id': experiment_id,
        'name': data.get('name'),
        'description': data.get('description'),
        'protocol_id': data.get('protocol_id'),
        'tissue_type_id': data.get('tissue_type_id'),
        'experiment_type': data.get('experiment_type'),
        'start_date': data.get('start_date'),
        'end_date': data.get('end_date'),
        'status': data.get('status'),
        'researcher': data.get('researcher'),
        'lab_id': data.get('lab_id'),
        'equipment': data.get('equipment', []),
        'environmental_conditions': data.get('environmental_conditions', {}),
        'notes': data.get('notes'),
        'tags': data.get('tags', []),
        'provenance': generate_provenance(request)
    }
    
    # Save to database
    db['experiments'].append(experiment)
    
    return jsonify(experiment), 201

@experiment_bp.route('/<experiment_id>', methods=['GET'])
def get_experiment(experiment_id):
    """Get experiment details by ID"""
    # Find experiment by ID
    experiment = next((e for e in db['experiments'] if e['id'] == experiment_id), None)
    
    if not experiment:
        return jsonify({'message': 'Experiment not found'}), 404
    
    # Add results and protocol information
    experiment = dict(experiment)  # Create a copy to avoid modifying the original
    experiment['results'] = [r for r in db['results'] if r['experiment_id'] == experiment_id]
    experiment['time_series'] = [ts for ts in db['time_series'] if ts['experiment_id'] == experiment_id]
    
    # Add protocol information
    protocol_id = experiment['protocol_id']
    protocol = next((p for p in db['protocols'] if p['id'] == protocol_id), None)
    if protocol:
        experiment['protocol'] = protocol
    
    return jsonify(experiment)

@experiment_bp.route('/<experiment_id>', methods=['PATCH'])
def update_experiment(experiment_id):
    """Update an existing experiment"""
    data = request.get_json()
    
    # Find experiment by ID
    experiment = next((e for e in db['experiments'] if e['id'] == experiment_id), None)
    
    if not experiment:
        return jsonify({'message': 'Experiment not found'}), 404
    
    # Update fields
    for key, value in data.items():
        if key != 'id' and key != 'provenance':  # Don't allow changing ID or provenance
            experiment[key] = value
    
    return jsonify(experiment)

@experiment_bp.route('/<experiment_id>', methods=['DELETE'])
def delete_experiment(experiment_id):
    """Delete an experiment"""
    # Find experiment by ID
    experiment_idx = next((i for i, e in enumerate(db['experiments']) if e['id'] == experiment_id), None)
    
    if experiment_idx is None:
        return jsonify({'message': 'Experiment not found'}), 404
    
    # Remove from database
    db['experiments'].pop(experiment_idx)
    
    # Remove related results and time series
    db['results'] = [r for r in db['results'] if r['experiment_id'] != experiment_id]
    db['time_series'] = [ts for ts in db['time_series'] if ts['experiment_id'] != experiment_id]
    
    return jsonify({'message': 'Experiment deleted successfully'}), 200

@experiment_bp.route('/<experiment_id>/results', methods=['GET'])
def get_experiment_results(experiment_id):
    """Get all results for an experiment"""
    # Find experiment by ID
    experiment = next((e for e in db['experiments'] if e['id'] == experiment_id), None)
    
    if not experiment:
        return jsonify({'message': 'Experiment not found'}), 404
    
    # Get results for this experiment
    results = [r for r in db['results'] if r['experiment_id'] == experiment_id]
    
    return jsonify(results)

@experiment_bp.route('/<experiment_id>/results', methods=['POST'])
def add_experiment_result(experiment_id):
    """Add a result to an experiment"""
    data = request.get_json()
    
    # Find experiment by ID
    experiment = next((e for e in db['experiments'] if e['id'] == experiment_id), None)
    
    if not experiment:
        return jsonify({'message': 'Experiment not found'}), 404
    
    # Create new result with generated ID
    result_id = str(uuid.uuid4())
    result = {
        'id': result_id,
        'experiment_id': experiment_id,
        'tissue_type_id': data.get('tissue_type_id'),
        'molecule_id': data.get('molecule_id'),
        'mixture_id': data.get('mixture_id'),
        'concentration': data.get('concentration'),
        'concentration_unit': data.get('concentration_unit'),
        'viability_percentage': data.get('viability_percentage'),
        'recovery_rate': data.get('recovery_rate'),
        'functionality_score': data.get('functionality_score'),
        'uncertainty': data.get('uncertainty', {}),
        'result_details': data.get('result_details', {}),
        'notes': data.get('notes'),
        'protocol_step_id': data.get('protocol_step_id'),
        'timestamp': datetime.now().isoformat(),
        'provenance': generate_provenance(request)
    }
    
    # Save to database
    db['results'].append(result)
    
    return jsonify(result), 201

@experiment_bp.route('/results/<result_id>', methods=['GET'])
def get_result(result_id):
    """Get a specific experiment result"""
    # Find result by ID
    result = next((r for r in db['results'] if r['id'] == result_id), None)
    
    if not result:
        return jsonify({'message': 'Result not found'}), 404
    
    return jsonify(result)

@experiment_bp.route('/results/<result_id>', methods=['PATCH'])
def update_result(result_id):
    """Update an experiment result"""
    data = request.get_json()
    
    # Find result by ID
    result = next((r for r in db['results'] if r['id'] == result_id), None)
    
    if not result:
        return jsonify({'message': 'Result not found'}), 404
    
    # Update fields
    for key, value in data.items():
        if key != 'id' and key != 'provenance' and key != 'experiment_id':  # Don't allow changing these
            result[key] = value
    
    return jsonify(result)

@experiment_bp.route('/results/<result_id>', methods=['DELETE'])
def delete_result(result_id):
    """Delete an experiment result"""
    # Find result by ID
    result_idx = next((i for i, r in enumerate(db['results']) if r['id'] == result_id), None)
    
    if result_idx is None:
        return jsonify({'message': 'Result not found'}), 404
    
    # Remove from database
    db['results'].pop(result_idx)
    
    return jsonify({'message': 'Result deleted successfully'}), 200

@experiment_bp.route('/<experiment_id>/timeseries', methods=['GET'])
def get_experiment_timeseries(experiment_id):
    """Get time series data for an experiment"""
    # Find experiment by ID
    experiment = next((e for e in db['experiments'] if e['id'] == experiment_id), None)
    
    if not experiment:
        return jsonify({'message': 'Experiment not found'}), 404
    
    # Get parameter from query
    parameter = request.args.get('parameter')
    
    # Get time series for this experiment
    time_series = [ts for ts in db['time_series'] if ts['experiment_id'] == experiment_id]
    
    # Filter by parameter if specified
    if parameter:
        time_series = [ts for ts in time_series if ts['parameter'] == parameter]
    
    return jsonify(time_series)

@experiment_bp.route('/analyze', methods=['POST'])
def analyze_experiments():
    """Analyze experimental data"""
    data = request.get_json()
    experiment_ids = data.get('experiment_ids', [])
    
    # Find experiments by IDs
    experiments = [e for e in db['experiments'] if e['id'] in experiment_ids]
    
    if not experiments:
        return jsonify({'message': 'No experiments found for the provided IDs'}), 404
    
    # Get results for these experiments
    results = [r for r in db['results'] if r['experiment_id'] in experiment_ids]
    
    # Mock analysis for demonstration
    successful_count = sum(1 for r in results if r.get('viability_percentage', 0) > 50)
    total_count = len(results)
    
    # Calculate simple statistics
    viability_values = [r.get('viability_percentage', 0) for r in results if r.get('viability_percentage') is not None]
    recovery_values = [r.get('recovery_rate', 0) for r in results if r.get('recovery_rate') is not None]
    
    if not viability_values:
        viability_values = [0]
    if not recovery_values:
        recovery_values = [0]
    
    viability_values.sort()
    recovery_values.sort()
    
    # Simple calculation of statistics
    viability_stats = {
        'mean': sum(viability_values) / len(viability_values),
        'median': viability_values[len(viability_values) // 2],
        'std_dev': 0,  # Would calculate standard deviation in a real implementation
        'min': min(viability_values),
        'max': max(viability_values),
        'quartiles': [
            viability_values[len(viability_values) // 4],
            viability_values[len(viability_values) // 2],
            viability_values[3 * len(viability_values) // 4] if len(viability_values) > 2 else viability_values[-1]
        ]
    }
    
    recovery_stats = {
        'mean': sum(recovery_values) / len(recovery_values),
        'median': recovery_values[len(recovery_values) // 2],
        'std_dev': 0,  # Would calculate standard deviation in a real implementation
        'min': min(recovery_values),
        'max': max(recovery_values),
        'quartiles': [
            recovery_values[len(recovery_values) // 4],
            recovery_values[len(recovery_values) // 2],
            recovery_values[3 * len(recovery_values) // 4] if len(recovery_values) > 2 else recovery_values[-1]
        ]
    }
    
    # Mock trend data
    experiment_dates = sorted([datetime.fromisoformat(e['start_date'].replace('Z', '+00:00')) if isinstance(e['start_date'], str) else e['start_date'] for e in experiments])
    experiment_dates_str = [d.isoformat() for d in experiment_dates]
    
    viability_trend = {
        'parameter': 'viability',
        'timestamps': experiment_dates_str,
        'values': [r.get('viability_percentage', 0) for r in results[:len(experiment_dates)]],
        'trend_line': [v + (i * 0.5) for i, v in enumerate([50, 55, 60, 65, 70][:len(experiment_dates)])]
    }
    
    analysis_result = {
        'summary': {
            'total': total_count,
            'successful': successful_count,
            'failed': total_count - successful_count,
            'success_rate': (successful_count / total_count) if total_count > 0 else 0
        },
        'statistics': {
            'viability': viability_stats,
            'recovery': recovery_stats
        },
        'trends': [viability_trend]
    }
    
    # Add comparisons if requested
    if data.get('compare_with'):
        analysis_result['comparisons'] = {
            'improvement': '15%',
            'significance': 'p < 0.05',
            'detailed_metrics': {
                'viability_change': '+10.5%',
                'recovery_change': '+12.3%'
            }
        }
    
    return jsonify(analysis_result)

@experiment_bp.route('/search', methods=['GET'])
def search_experiments():
    """Search for experiments"""
    query = request.args.get('query', '').lower()
    
    # Filter experiments by query
    filtered_experiments = [
        e for e in db['experiments'] 
        if query in e['name'].lower() or 
        (e.get('description') and query in e['description'].lower()) or
        (e.get('notes') and query in e['notes'].lower()) or
        query in e['researcher'].lower()
    ]
    
    return jsonify({
        'data': filtered_experiments,
        'total': len(filtered_experiments)
    })

# Register the blueprint with a URL prefix
app.register_blueprint(experiment_bp, url_prefix='/api/experiments')

# Home route
@app.route('/')
def home():
    """Home route"""
    return jsonify({
        'message': 'Experiment API Test Server',
        'endpoints': {
            'GET /api/experiments': 'List all experiments',
            'POST /api/experiments': 'Create a new experiment',
            'GET /api/experiments/<id>': 'Get experiment details',
            'PATCH /api/experiments/<id>': 'Update an experiment',
            'DELETE /api/experiments/<id>': 'Delete an experiment',
            'GET /api/experiments/<id>/results': 'Get experiment results',
            'POST /api/experiments/<id>/results': 'Add a result to an experiment',
            'GET /api/experiments/results/<id>': 'Get a specific result',
            'PATCH /api/experiments/results/<id>': 'Update a result',
            'DELETE /api/experiments/results/<id>': 'Delete a result',
            'GET /api/experiments/<id>/timeseries': 'Get experiment time series data',
            'POST /api/experiments/analyze': 'Analyze experiments',
            'GET /api/experiments/search': 'Search for experiments'
        },
        'test_protocol_id': get_test_protocol_id(),
        'test_tissue_type_id': get_test_tissue_type_id()
    })

# Helper endpoint to get test IDs
@app.route('/api/test-ids')
def test_ids():
    """Get test protocol and tissue type IDs"""
    return jsonify({
        'protocol_id': get_test_protocol_id(),
        'tissue_type_id': get_test_tissue_type_id()
    })

if __name__ == '__main__':
    # Create data directory if it doesn't exist
    os.makedirs('data', exist_ok=True)
    
    # Run the app
    app.run(debug=True, port=5000, host='0.0.0.0')