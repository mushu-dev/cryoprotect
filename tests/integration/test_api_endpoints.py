"""
Integration tests for API endpoints.

This module contains integration tests for the main API endpoints
to ensure they function correctly with the database and application context.
"""

import pytest
import requests

BASE_URL = "http://localhost:5000"

@pytest.mark.integration
def test_get_root():
    """
    Test the root endpoint returns a 200 status and expected message.
    """
    response = requests.get(f"{BASE_URL}/")
    assert response.status_code == 200
    assert "CryoProtect API" in response.text

@pytest.mark.integration
def test_get_molecules():
    """
    Test the /molecules endpoint returns a list of molecules.
    """
    response = requests.get(f"{BASE_URL}/molecules")
    assert response.status_code == 200
    data = response.json()
    assert isinstance(data, list)
    assert len(data) > 0
    assert "molecule_id" in data[0]
    assert "name" in data[0]

@pytest.mark.integration
def test_get_molecule_by_id():
    """
    Test retrieving a single molecule by ID.
    """
    # First, get a valid molecule ID
    response = requests.get(f"{BASE_URL}/molecules")
    assert response.status_code == 200
    data = response.json()
    molecule_id = data[0]["molecule_id"]

    # Now, get the molecule by ID
    response = requests.get(f"{BASE_URL}/molecules/{molecule_id}")
    assert response.status_code == 200
    molecule = response.json()
    assert molecule["molecule_id"] == molecule_id
    assert "name" in molecule

@pytest.mark.integration
def test_post_new_molecule():
    """
    Test creating a new molecule via POST.
    """
    new_molecule = {
        "name": "IntegrationTestMolecule",
        "formula": "C2H6O",
        "description": "Test molecule for integration testing"
    }
    response = requests.post(f"{BASE_URL}/molecules", json=new_molecule)
    assert response.status_code == 201
    created = response.json()
    assert created["name"] == new_molecule["name"]
    assert "molecule_id" in created

    # Cleanup: delete the created molecule if DELETE endpoint exists
    # requests.delete(f"{BASE_URL}/molecules/{created['molecule_id']}")

@pytest.mark.integration
def test_invalid_endpoint_returns_404():
    """
    Test that an invalid endpoint returns a 404 status code.
    """
    response = requests.get(f"{BASE_URL}/nonexistent")
    assert response.status_code == 404