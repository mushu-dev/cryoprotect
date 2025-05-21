#!/usr/bin/env python3
"""
Script to test the Experiment API integration with real data.

This script tests the following:
1. Creating a new experiment
2. Adding experiment results
3. Retrieving experiment details
4. Updating experiment data
5. Running experiment analysis
"""

import requests
import json
import uuid
from datetime import datetime, timedelta
import sys

# Base URL for API
BASE_URL = 'http://localhost:5000/api'  # Change this to your API endpoint
EXPERIMENTS_ENDPOINT = f'{BASE_URL}/experiments'

# Authentication headers (adjust based on your authentication method)
# If using JWT or other token-based auth, you'll need to add the token here
HEADERS = {
    'Content-Type': 'application/json',
    'Accept': 'application/json'
}

# Test data for an experiment
def generate_test_experiment():
    """Generate test data for an experiment"""
    # Get test IDs from the API
    test_ids_response = requests.get(f"{BASE_URL}/test-ids")
    test_ids = test_ids_response.json()
    
    protocol_id = test_ids['protocol_id']
    tissue_type_id = test_ids['tissue_type_id']
    
    return {
        'name': f'Test Experiment {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}',
        'description': 'This is a test experiment created by the integration test script',
        'protocol_id': protocol_id,
        'tissue_type_id': tissue_type_id,
        'experiment_type': 'vitrification',
        'start_date': datetime.now().strftime('%Y-%m-%d'),
        'status': 'in_progress',
        'researcher': 'Test Researcher',
        'lab_id': 'Test Lab',
        'equipment': ['Microscope', 'Centrifuge'],
        'environmental_conditions': {
            'temperature': 23.5,
            'humidity': 45.2,
            'pressure': 1013.2
        },
        'notes': 'Test notes',
        'tags': ['test', 'integration', 'api']
    }

# Test data for experiment result
def generate_test_result(experiment_id, tissue_type_id):
    """Generate test data for an experiment result"""
    # In a real implementation, these would be actual IDs from your database
    molecule_id = str(uuid.uuid4())  # Replace with a real molecule ID
    
    return {
        'experiment_id': experiment_id,
        'tissue_type_id': tissue_type_id,
        'molecule_id': molecule_id,
        'concentration': 5.0,
        'concentration_unit': 'mM',
        'viability_percentage': 78.5,
        'recovery_rate': 82.3,
        'functionality_score': 75.0,
        'uncertainty': {
            'viability_percentage': {
                'value': 2.5,
                'type': 'standard',
                'confidence': 0.95,
                'distribution': 'normal'
            }
        },
        'result_details': {
            'temperature': 22.5,
            'duration': 120,
            'duration_unit': 'minutes'
        },
        'notes': 'Test result notes'
    }

def test_create_experiment():
    """Test creating a new experiment"""
    print("Testing experiment creation...")
    
    experiment_data = generate_test_experiment()
    response = requests.post(EXPERIMENTS_ENDPOINT, headers=HEADERS, json=experiment_data)
    
    if response.status_code == 201:
        print("✅ Experiment created successfully")
        experiment = response.json()
        print(f"   ID: {experiment['id']}")
        print(f"   Name: {experiment['name']}")
        return experiment
    else:
        print(f"❌ Failed to create experiment: {response.status_code} - {response.text}")
        return None

def test_get_experiment(experiment_id):
    """Test retrieving an experiment by ID"""
    print(f"\nTesting experiment retrieval for ID: {experiment_id}...")
    
    response = requests.get(f"{EXPERIMENTS_ENDPOINT}/{experiment_id}", headers=HEADERS)
    
    if response.status_code == 200:
        print("✅ Experiment retrieved successfully")
        experiment = response.json()
        print(f"   Name: {experiment['name']}")
        print(f"   Status: {experiment['status']}")
        return experiment
    else:
        print(f"❌ Failed to retrieve experiment: {response.status_code} - {response.text}")
        return None

def test_add_experiment_result(experiment_id, tissue_type_id):
    """Test adding a result to an experiment"""
    print(f"\nTesting result addition for experiment ID: {experiment_id}...")
    
    result_data = generate_test_result(experiment_id, tissue_type_id)
    response = requests.post(f"{EXPERIMENTS_ENDPOINT}/{experiment_id}/results", headers=HEADERS, json=result_data)
    
    if response.status_code == 201:
        print("✅ Result added successfully")
        result = response.json()
        print(f"   Result ID: {result['id']}")
        print(f"   Viability: {result['viability_percentage']}%")
        return result
    else:
        print(f"❌ Failed to add result: {response.status_code} - {response.text}")
        return None

def test_update_experiment(experiment_id):
    """Test updating an experiment"""
    print(f"\nTesting experiment update for ID: {experiment_id}...")
    
    update_data = {
        'status': 'completed',
        'end_date': datetime.now().strftime('%Y-%m-%d'),
        'notes': 'Updated notes from integration test'
    }
    
    response = requests.patch(f"{EXPERIMENTS_ENDPOINT}/{experiment_id}", headers=HEADERS, json=update_data)
    
    if response.status_code == 200:
        print("✅ Experiment updated successfully")
        experiment = response.json()
        print(f"   Status: {experiment['status']}")
        print(f"   End Date: {experiment['end_date']}")
        return experiment
    else:
        print(f"❌ Failed to update experiment: {response.status_code} - {response.text}")
        return None

def test_analyze_experiments(experiment_ids):
    """Test analyzing experiments"""
    print(f"\nTesting experiment analysis for IDs: {experiment_ids}...")
    
    analysis_data = {
        'experiment_ids': experiment_ids,
        'analysis_type': ['viability', 'recovery']
    }
    
    response = requests.post(f"{EXPERIMENTS_ENDPOINT}/analyze", headers=HEADERS, json=analysis_data)
    
    if response.status_code == 200:
        print("✅ Analysis completed successfully")
        analysis = response.json()
        print(f"   Success Rate: {analysis['summary']['success_rate'] * 100:.1f}%")
        print(f"   Mean Viability: {analysis['statistics']['viability']['mean']:.1f}%")
        return analysis
    else:
        print(f"❌ Failed to analyze experiments: {response.status_code} - {response.text}")
        return None

def test_search_experiments(query):
    """Test searching for experiments"""
    print(f"\nTesting experiment search for query: '{query}'...")
    
    response = requests.get(f"{EXPERIMENTS_ENDPOINT}/search?query={query}", headers=HEADERS)
    
    if response.status_code == 200:
        print("✅ Search completed successfully")
        results = response.json()
        print(f"   Found {results['total']} experiments")
        for i, exp in enumerate(results['data'][:3]):  # Show up to 3 results
            print(f"   {i+1}. {exp['name']} - {exp['status']}")
        return results
    else:
        print(f"❌ Failed to search experiments: {response.status_code} - {response.text}")
        return None

def test_list_experiments():
    """Test listing experiments with pagination and filtering"""
    print("\nTesting experiment listing...")
    
    params = {
        'page': 1,
        'per_page': 5,
        'sort_by': 'start_date',
        'sort_order': 'desc'
    }
    
    response = requests.get(EXPERIMENTS_ENDPOINT, headers=HEADERS, params=params)
    
    if response.status_code == 200:
        print("✅ Experiment listing successful")
        results = response.json()
        print(f"   Total: {results['total']} experiments")
        print(f"   Page: {results['page']} of {results['total_pages']}")
        for i, exp in enumerate(results['data']):
            print(f"   {i+1}. {exp['name']} - {exp['status']}")
        return results
    else:
        print(f"❌ Failed to list experiments: {response.status_code} - {response.text}")
        return None

def main():
    """Main test function that runs all the tests"""
    print("==== Experiment API Integration Test ====\n")
    
    # First, create a new experiment
    experiment = test_create_experiment()
    if not experiment:
        print("Cannot continue testing without a valid experiment")
        sys.exit(1)
    
    experiment_id = experiment['id']
    tissue_type_id = experiment['tissue_type_id']
    
    # Test retrieving the experiment
    retrieved_experiment = test_get_experiment(experiment_id)
    
    # Test adding a result to the experiment
    result = test_add_experiment_result(experiment_id, tissue_type_id)
    
    # Test updating the experiment
    updated_experiment = test_update_experiment(experiment_id)
    
    # Test analyzing the experiment
    analysis = test_analyze_experiments([experiment_id])
    
    # Test searching for experiments
    search_results = test_search_experiments("Test")
    
    # Test listing experiments
    list_results = test_list_experiments()
    
    print("\n==== All Tests Completed ====")

if __name__ == "__main__":
    main()