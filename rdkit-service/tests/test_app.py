"""
Tests for the RDKit microservice API.
"""

import json
import pytest
import sys
import os

# Add the app directory to the path so we can import modules
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from app.app import app

@pytest.fixture
def client():
    app.config['TESTING'] = True
    with app.test_client() as client:
        yield client

def test_index(client):
    """Test the index route."""
    response = client.get('/')
    data = json.loads(response.data)
    
    assert response.status_code == 200
    assert data['status'] == 'ok'
    assert 'service' in data
    assert 'version' in data

def test_health_check(client):
    """Test the health check route."""
    response = client.get('/health')
    data = json.loads(response.data)
    
    assert response.status_code == 200
    assert data['status'] == 'ok'
    assert 'rdkit_version' in data
    assert 'timestamp' in data

def test_calculate_properties(client):
    """Test calculating molecular properties."""
    # Test with valid molecule
    response = client.post(
        '/api/calculate-properties',
        data=json.dumps({'molecule_data': 'CCO'}),
        content_type='application/json'
    )
    data = json.loads(response.data)
    
    assert response.status_code == 200
    assert data['status'] == 'success'
    assert 'data' in data
    assert 'logp' in data['data']
    assert 'hydrogen_bonding' in data['data']
    
    # Test with invalid molecule
    response = client.post(
        '/api/calculate-properties',
        data=json.dumps({'molecule_data': 'invalid_smiles'}),
        content_type='application/json'
    )
    data = json.loads(response.data)
    
    assert response.status_code == 400
    assert 'error' in data

def test_visualization(client):
    """Test generating a molecule visualization."""
    response = client.post(
        '/api/visualization',
        data=json.dumps({
            'molecule_data': 'CCO',
            'width': 300,
            'height': 200
        }),
        content_type='application/json'
    )
    data = json.loads(response.data)
    
    assert response.status_code == 200
    assert data['status'] == 'success'
    assert 'data' in data
    assert 'svg' in data['data']
    assert data['data']['width'] == 300
    assert data['data']['height'] == 200

def test_substructure_search(client):
    """Test substructure search functionality."""
    response = client.post(
        '/api/substructure-search',
        data=json.dumps({
            'query_mol_data': '[OH]',
            'target_mol_data': 'CCO',
            'query_format': 'smarts',
            'target_format': 'smiles'
        }),
        content_type='application/json'
    )
    data = json.loads(response.data)
    
    assert response.status_code == 200
    assert data['status'] == 'success'
    assert 'data' in data
    assert data['data']['match'] is True

def test_similarity(client):
    """Test molecular similarity calculation."""
    response = client.post(
        '/api/similarity',
        data=json.dumps({
            'mol1_data': 'CCO',
            'mol2_data': 'CC(=O)O',
            'fingerprint_type': 'morgan'
        }),
        content_type='application/json'
    )
    data = json.loads(response.data)
    
    assert response.status_code == 200
    assert data['status'] == 'success'
    assert 'data' in data
    assert 'tanimoto' in data['data']
    assert 'dice' in data['data']
    assert 'fingerprint_type' in data['data']