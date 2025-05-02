#!/usr/bin/env python
"""
CryoProtect API Examples - Python

This file contains examples of how to use the CryoProtect API with Python.
It demonstrates common operations such as authentication, working with molecules,
mixtures, predictions, RDKit integration, scoring, and export functionality.
"""

import requests
import time
import random
import os
from datetime import datetime
import json

# Base URL for the API
BASE_URL = "http://localhost:5000"

# Authentication token (would be obtained through login)
TOKEN = None

# Helper function to get authentication token
def get_token():
    """Get the authentication token."""
    return TOKEN

# Authentication Examples
def sign_up(email, password):
    """Sign up a new user."""
    try:
        response = requests.post(
            f"{BASE_URL}/auth/register",
            headers={"Content-Type": "application/json"},
            json={
                "email": email,
                "password": password
            }
        )
        
        response.raise_for_status()
        return response.json()
    except requests.exceptions.RequestException as e:
        print(f"Error signing up: {str(e)}")
        raise

def sign_in(email, password):
    """Sign in an existing user and get an authentication token."""
    try:
        response = requests.post(
            f"{BASE_URL}/auth/login",
            headers={"Content-Type": "application/json"},
            json={
                "email": email,
                "password": password
            }
        )
        
        response.raise_for_status()
        data = response.json()
        
        # Store token for future requests
        global TOKEN
        TOKEN = data.get("access_token")
        
        return data
    except requests.exceptions.RequestException as e:
        print(f"Error signing in: {str(e)}")
        raise

def sign_out():
    """Sign out the current user."""
    try:
        response = requests.post(
            f"{BASE_URL}/auth/logout",
            headers={
                "Content-Type": "application/json",
                "Authorization": f"Bearer {get_token()}"
            }
        )
        
        response.raise_for_status()
        
        # Clear token
        global TOKEN
        TOKEN = None
        
        return {"success": True}
    except requests.exceptions.RequestException as e:
        print(f"Error signing out: {str(e)}")
        raise

# Molecule Examples
def get_molecules():
    """Get a list of all molecules."""
    try:
        response = requests.get(
            f"{BASE_URL}/api/v1/molecules",
            headers={
                "Content-Type": "application/json",
                "Authorization": f"Bearer {get_token()}"
            }
        )
        
        response.raise_for_status()
        return response.json()
    except requests.exceptions.RequestException as e:
        print(f"Error fetching molecules: {str(e)}")
        raise

def create_molecule(molecule_data):
    """Create a new molecule."""
    try:
        response = requests.post(
            f"{BASE_URL}/api/v1/molecules",
            headers={
                "Content-Type": "application/json",
                "Authorization": f"Bearer {get_token()}"
            },
            json=molecule_data
        )
        
        response.raise_for_status()
        return response.json()
    except requests.exceptions.RequestException as e:
        print(f"Error creating molecule: {str(e)}")
        raise

def update_molecule(molecule_id, updates):
    """Update an existing molecule."""
    try:
        response = requests.put(
            f"{BASE_URL}/api/v1/molecules/{molecule_id}",
            headers={
                "Content-Type": "application/json",
                "Authorization": f"Bearer {get_token()}"
            },
            json=updates
        )
        
        response.raise_for_status()
        return response.json()
    except requests.exceptions.RequestException as e:
        print(f"Error updating molecule: {str(e)}")
        raise

def delete_molecule(molecule_id):
    """Delete a molecule."""
    try:
        response = requests.delete(
            f"{BASE_URL}/api/v1/molecules/{molecule_id}",
            headers={
                "Content-Type": "application/json",
                "Authorization": f"Bearer {get_token()}"
            }
        )
        
        response.raise_for_status()
        return {"success": True}
    except requests.exceptions.RequestException as e:
        print(f"Error deleting molecule: {str(e)}")
        raise

# Mixture Examples
def get_mixtures():
    """Get a list of all mixtures."""
    try:
        response = requests.get(
            f"{BASE_URL}/api/v1/mixtures",
            headers={
                "Content-Type": "application/json",
                "Authorization": f"Bearer {get_token()}"
            }
        )
        
        response.raise_for_status()
        return response.json()
    except requests.exceptions.RequestException as e:
        print(f"Error fetching mixtures: {str(e)}")
        raise

def create_mixture(mixture_data, components):
    """Create a new mixture with components."""
    try:
        # Create the mixture
        response = requests.post(
            f"{BASE_URL}/api/v1/mixtures",
            headers={
                "Content-Type": "application/json",
                "Authorization": f"Bearer {get_token()}"
            },
            json=mixture_data
        )
        
        response.raise_for_status()
        mixture = response.json()
        
        # Add components
        mixture_id = mixture.get("id")
        
        for component in components:
            component["mixture_id"] = mixture_id
            component_response = requests.post(
                f"{BASE_URL}/api/v1/mixture_components",
                headers={
                    "Content-Type": "application/json",
                    "Authorization": f"Bearer {get_token()}"
                },
                json=component
            )
            component_response.raise_for_status()
        
        # Get the complete mixture with components
        complete_response = requests.get(
            f"{BASE_URL}/api/v1/mixtures/{mixture_id}",
            headers={
                "Content-Type": "application/json",
                "Authorization": f"Bearer {get_token()}"
            }
        )
        
        complete_response.raise_for_status()
        return complete_response.json()
    except requests.exceptions.RequestException as e:
        print(f"Error creating mixture: {str(e)}")
        raise

# Prediction Examples
def create_prediction(prediction_data):
    """Create a new prediction."""
    try:
        response = requests.post(
            f"{BASE_URL}/api/v1/predictions",
            headers={
                "Content-Type": "application/json",
                "Authorization": f"Bearer {get_token()}"
            },
            json=prediction_data
        )
        
        response.raise_for_status()
        return response.json()
    except requests.exceptions.RequestException as e:
        print(f"Error creating prediction: {str(e)}")
        raise

def get_predictions(mixture_id=None, molecule_id=None):
    """Get predictions, optionally filtered by mixture or molecule ID."""
    try:
        params = {}
        if mixture_id:
            params["mixture_id"] = mixture_id
        elif molecule_id:
            params["molecule_id"] = molecule_id
        
        response = requests.get(
            f"{BASE_URL}/api/v1/predictions",
            headers={
                "Content-Type": "application/json",
                "Authorization": f"Bearer {get_token()}"
            },
            params=params
        )
        
        response.raise_for_status()
        return response.json()
    except requests.exceptions.RequestException as e:
        print(f"Error fetching predictions: {str(e)}")
        raise

# RDKit Integration Examples
def calculate_properties(molecule_data, input_format='smiles'):
    """Calculate molecular properties using RDKit."""
    try:
        response = requests.post(
            f"{BASE_URL}/api/v1/rdkit/properties",
            headers={
                "Content-Type": "application/json",
                "Authorization": f"Bearer {get_token()}"
            },
            json={
                "molecule_data": molecule_data,
                "input_format": input_format
            }
        )
        
        response.raise_for_status()
        return response.json()
    except requests.exceptions.RequestException as e:
        print(f"Error calculating molecular properties: {str(e)}")
        raise

def generate_visualization(molecule_data, input_format='smiles', width=400, height=300, highlight_atoms=None):
    """Generate a visualization of a molecule."""
    try:
        response = requests.post(
            f"{BASE_URL}/api/v1/rdkit/visualization",
            headers={
                "Content-Type": "application/json",
                "Authorization": f"Bearer {get_token()}"
            },
            json={
                "molecule_data": molecule_data,
                "input_format": input_format,
                "width": width,
                "height": height,
                "highlight_atoms": highlight_atoms or []
            }
        )
        
        response.raise_for_status()
        return response.json()
    except requests.exceptions.RequestException as e:
        print(f"Error generating molecular visualization: {str(e)}")
        raise

def perform_substructure_search(query_mol_data, target_mol_data, query_format='smarts', target_format='smiles'):
    """Search for a substructure within a molecule."""
    try:
        response = requests.post(
            f"{BASE_URL}/api/v1/rdkit/substructure",
            headers={
                "Content-Type": "application/json",
                "Authorization": f"Bearer {get_token()}"
            },
            json={
                "query_mol_data": query_mol_data,
                "target_mol_data": target_mol_data,
                "query_format": query_format,
                "target_format": target_format
            }
        )
        
        response.raise_for_status()
        return response.json()
    except requests.exceptions.RequestException as e:
        print(f"Error performing substructure search: {str(e)}")
        raise

def calculate_similarity(mol1_data, mol2_data, mol1_format='smiles', mol2_format='smiles', fingerprint_type='morgan'):
    """Calculate similarity between two molecules."""
    try:
        response = requests.post(
            f"{BASE_URL}/api/v1/rdkit/similarity",
            headers={
                "Content-Type": "application/json",
                "Authorization": f"Bearer {get_token()}"
            },
            json={
                "mol1_data": mol1_data,
                "mol2_data": mol2_data,
                "mol1_format": mol1_format,
                "mol2_format": mol2_format,
                "fingerprint_type": fingerprint_type
            }
        )
        
        response.raise_for_status()
        return response.json()
    except requests.exceptions.RequestException as e:
        print(f"Error calculating molecular similarity: {str(e)}")
        raise

# Scoring Examples
def score_molecule(molecule_data, input_format='smiles', store_result=False):
    """Score a molecule based on its cryoprotection effectiveness."""
    try:
        response = requests.post(
            f"{BASE_URL}/api/v1/scoring/molecules",
            headers={
                "Content-Type": "application/json",
                "Authorization": f"Bearer {get_token()}"
            },
            json={
                "molecule_data": molecule_data,
                "input_format": input_format,
                "store_result": store_result
            }
        )
        
        response.raise_for_status()
        return response.json()
    except requests.exceptions.RequestException as e:
        print(f"Error scoring molecule: {str(e)}")
        raise

def score_molecule_by_id(molecule_id, store_result=True):
    """Score a molecule by ID."""
    try:
        response = requests.post(
            f"{BASE_URL}/api/v1/molecules/{molecule_id}/score",
            headers={
                "Content-Type": "application/json",
                "Authorization": f"Bearer {get_token()}"
            },
            json={
                "store_result": store_result
            }
        )
        
        response.raise_for_status()
        return response.json()
    except requests.exceptions.RequestException as e:
        print(f"Error scoring molecule by ID: {str(e)}")
        raise

def score_mixture(mixture_id, store_result=True):
    """Score a mixture based on its cryoprotection effectiveness."""
    try:
        response = requests.post(
            f"{BASE_URL}/api/v1/mixtures/{mixture_id}/score",
            headers={
                "Content-Type": "application/json",
                "Authorization": f"Bearer {get_token()}"
            },
            json={
                "store_result": store_result
            }
        )
        
        response.raise_for_status()
        return response.json()
    except requests.exceptions.RequestException as e:
        print(f"Error scoring mixture: {str(e)}")
        raise

def batch_score(entity_type, ids, store_results=True):
    """Perform batch scoring of multiple entities."""
    try:
        response = requests.post(
            f"{BASE_URL}/api/v1/scoring/batch",
            headers={
                "Content-Type": "application/json",
                "Authorization": f"Bearer {get_token()}"
            },
            json={
                "entity_type": entity_type,
                "ids": ids,
                "store_results": store_results
            }
        )
        
        response.raise_for_status()
        return response.json()
    except requests.exceptions.RequestException as e:
        print(f"Error performing batch scoring: {str(e)}")
        raise

# Export and Sharing Examples
def export_data(data_type, format, id=None, include_related=False):
    """Export data in various formats."""
    try:
        response = requests.post(
            f"{BASE_URL}/api/v1/export",
            headers={
                "Content-Type": "application/json",
                "Authorization": f"Bearer {get_token()}"
            },
            json={
                "data_type": data_type,
                "format": format,
                "id": id,
                "include_related": include_related
            },
            stream=True  # Stream the response for file downloads
        )
        
        response.raise_for_status()
        
        # Get filename from Content-Disposition header or create a default one
        if 'Content-Disposition' in response.headers:
            filename = response.headers['Content-Disposition'].split('filename=')[1].replace('"', '')
        else:
            filename = f"{data_type}_export.{format}"
        
        # Save the file
        with open(filename, 'wb') as f:
            for chunk in response.iter_content(chunk_size=8192):
                f.write(chunk)
        
        return {'success': True, 'filename': filename, 'path': os.path.abspath(filename)}
    except requests.exceptions.RequestException as e:
        print(f"Error exporting data: {str(e)}")
        raise

def share_results(share_type, data_type, id, **options):
    """Share results via link, email, or embed code."""
    try:
        request_body = {
            'share_type': share_type,
            'data_type': data_type,
            'id': id,
            **options
        }
        
        response = requests.post(
            f"{BASE_URL}/api/v1/share",
            headers={
                "Content-Type": "application/json",
                "Authorization": f"Bearer {get_token()}"
            },
            json=request_body
        )
        
        response.raise_for_status()
        return response.json()
    except requests.exceptions.RequestException as e:
        print(f"Error sharing results: {str(e)}")
        raise

def share_via_link(data_type, id, password_protected=False, password=None, expiration=86400):
    """Share results via a link."""
    return share_results('link', data_type, id, 
                        password_protected=password_protected,
                        password=password,
                        expiration=expiration)

def share_via_email(data_type, id, recipients, message='', password_protected=False, password=None):
    """Share results via email."""
    return share_results('email', data_type, id,
                        recipients=recipients,
                        message=message,
                        password_protected=password_protected,
                        password=password)

def get_embed_code(data_type, id):
    """Get an embed code for sharing results."""
    result = share_results('embed', data_type, id)
    return result.get('embed_code')

# Rate Limiting Examples
def handle_rate_limited_request(url, method='GET', headers=None, json=None, params=None):
    """Handle a request with rate limiting."""
    try:
        if method.upper() == 'GET':
            response = requests.get(url, headers=headers, params=params)
        elif method.upper() == 'POST':
            response = requests.post(url, headers=headers, json=json)
        else:
            raise ValueError(f"Unsupported method: {method}")
        
        # Check for rate limit headers
        rate_limit_limit = response.headers.get('X-RateLimit-Limit')
        rate_limit_remaining = response.headers.get('X-RateLimit-Remaining')
        rate_limit_reset = response.headers.get('X-RateLimit-Reset')
        
        if rate_limit_remaining and int(rate_limit_remaining) < 5:
            reset_time = datetime.fromtimestamp(int(rate_limit_reset)).strftime('%H:%M:%S')
            print(f"Rate limit warning: {rate_limit_remaining}/{rate_limit_limit} requests remaining. Resets at {reset_time}")
        
        if response.status_code == 429:
            # Rate limit exceeded
            retry_after = int(response.headers.get('Retry-After', 60))
            print(f"Rate limit exceeded. Retrying after {retry_after} seconds.")
            
            # Wait for the specified time and retry
            time.sleep(retry_after)
            return handle_rate_limited_request(url, method, headers, json, params)
        
        response.raise_for_status()
        return response.json()
    except requests.exceptions.RequestException as e:
        print(f"Request error: {str(e)}")
        raise

def fetch_with_backoff(url, method='GET', headers=None, json=None, params=None, max_retries=5):
    """Fetch data with exponential backoff for rate limiting."""
    retries = 0
    
    while retries < max_retries:
        try:
            if method.upper() == 'GET':
                response = requests.get(url, headers=headers, params=params)
            elif method.upper() == 'POST':
                response = requests.post(url, headers=headers, json=json)
            else:
                raise ValueError(f"Unsupported method: {method}")
            
            if response.status_code == 429:
                # Rate limit exceeded
                retry_after = int(response.headers.get('Retry-After', 60))
                backoff_time = retry_after * (2 ** retries)
                # Add jitter to avoid thundering herd problem
                jitter = random.uniform(0, 0.1 * backoff_time)
                backoff_time += jitter
                
                print(f"Rate limit exceeded. Retrying after {backoff_time:.1f} seconds (retry {retries + 1}/{max_retries}).")
                time.sleep(backoff_time)
                retries += 1
                continue
            
            response.raise_for_status()
            return response.json()
        except requests.exceptions.RequestException as e:
            if retries >= max_retries - 1:
                print(f"Max retries reached: {str(e)}")
                raise
            
            # For other errors, also implement backoff
            backoff_time = (2 ** retries) + random.uniform(0, 1)
            print(f"Request failed. Retrying after {backoff_time:.1f} seconds (retry {retries + 1}/{max_retries}).")
            time.sleep(backoff_time)
            retries += 1

# Example usage
if __name__ == "__main__":
    # This is a simple example of how to use the API
    # In a real application, you would handle errors and use proper authentication
    
    print("CryoProtect API Examples - Python")
    print("=================================")
    
    # Sign in
    print("\nSigning in...")
    try:
        # Replace with your credentials
        auth_data = sign_in("user@example.com", "password")
        print(f"Signed in successfully. Token: {TOKEN[:10]}...")
    except Exception as e:
        print(f"Failed to sign in: {str(e)}")
        exit(1)
    
    # Calculate properties for ethanol
    print("\nCalculating properties for ethanol...")
    try:
        properties = calculate_properties("CCO")
        print(f"LogP: {properties.get('logp')}")
        print(f"TPSA: {properties.get('tpsa')}")
        print(f"H-Bond Donors: {properties.get('hydrogen_bonding', {}).get('donors')}")
        print(f"H-Bond Acceptors: {properties.get('hydrogen_bonding', {}).get('acceptors')}")
    except Exception as e:
        print(f"Failed to calculate properties: {str(e)}")
    
    # Score ethanol
    print("\nScoring ethanol...")
    try:
        score = score_molecule("CCO")
        print(f"Overall Score: {score.get('overall_score')}")
        print(f"Component Scores: {json.dumps(score.get('component_scores', {}), indent=2)}")
    except Exception as e:
        print(f"Failed to score molecule: {str(e)}")
    
    # Sign out
    print("\nSigning out...")
    try:
        sign_out()
        print("Signed out successfully.")
    except Exception as e:
        print(f"Failed to sign out: {str(e)}")