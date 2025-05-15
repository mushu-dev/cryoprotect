"""
RDKit Service Client for Heroku Application

This module provides a client for communicating with the RDKit microservice
deployed on Fly.io from the main Heroku application.
"""

import os
import requests
import json
import logging
import time
from functools import lru_cache
from typing import Dict, Any, Optional, Union

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger('rdkit-client')

class RDKitServiceClient:
    """Client for the RDKit microservice."""
    
    def __init__(self, service_url: Optional[str] = None, api_key: Optional[str] = None, 
                 timeout: int = 30, max_retries: int = 3):
        """
        Initialize the RDKit service client.
        
        Args:
            service_url: URL of the RDKit service. If None, uses RDKIT_SERVICE_URL env var.
            api_key: API key for authentication. If None, uses RDKIT_API_KEY env var.
            timeout: Request timeout in seconds.
            max_retries: Maximum number of retries for failed requests.
        """
        self.service_url = service_url or os.environ.get('RDKIT_SERVICE_URL')
        if not self.service_url:
            logger.warning("RDKit service URL not provided. Using mock implementation.")
            self.use_mock = True
        else:
            self.use_mock = False
            
        self.api_key = api_key or os.environ.get('RDKIT_API_KEY', '')
        self.timeout = timeout
        self.max_retries = max_retries
        
    def _get_headers(self) -> Dict[str, str]:
        """Get the headers for API requests."""
        headers = {'Content-Type': 'application/json'}
        if self.api_key:
            headers['X-API-Key'] = self.api_key
        return headers
        
    def _make_request(self, endpoint: str, data: Dict[str, Any], method: str = 'POST') -> Dict[str, Any]:
        """
        Make a request to the RDKit service.
        
        Args:
            endpoint: API endpoint (without leading slash).
            data: Request payload.
            method: HTTP method (GET, POST, etc.)
            
        Returns:
            Response data as dictionary.
            
        Raises:
            Exception if the request fails after retries.
        """
        if self.use_mock:
            return self._mock_response(endpoint, data)
            
        url = f"{self.service_url.rstrip('/')}/api/{endpoint}"
        headers = self._get_headers()
        
        for attempt in range(self.max_retries):
            try:
                response = requests.request(
                    method=method,
                    url=url,
                    headers=headers,
                    json=data,
                    timeout=self.timeout
                )
                response.raise_for_status()
                return response.json()
            except requests.exceptions.RequestException as e:
                logger.warning(f"Request attempt {attempt + 1} failed: {str(e)}")
                if attempt == self.max_retries - 1:
                    logger.error(f"All {self.max_retries} request attempts failed")
                    raise
                time.sleep(1 * (attempt + 1))  # Exponential backoff
    
    def _mock_response(self, endpoint: str, data: Dict[str, Any]) -> Dict[str, Any]:
        """
        Generate mock responses when the service is unavailable.
        
        Args:
            endpoint: API endpoint.
            data: Request payload.
            
        Returns:
            Mock response data.
        """
        logger.info(f"Using mock implementation for endpoint: {endpoint}")
        
        if endpoint == 'calculate-properties':
            return {
                "status": "success",
                "data": {
                    "hydrogen_bonding": {"donors": 1, "acceptors": 1, "total": 2},
                    "logp": 0.14,
                    "tpsa": 20.23,
                    "molecular_properties": {
                        "molecular_weight": 46.07,
                        "heavy_atom_count": 3,
                        "rotatable_bond_count": 1
                    },
                    "functional_groups": {"alcohol": 1, "hydroxyl": 1},
                    "smiles": data.get('molecule_data', 'CCO'),
                    "inchi": "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3",
                    "inchi_key": "LFQSCWFLJHTTHZ-UHFFFAOYSA-N"
                }
            }
        elif endpoint == 'visualization':
            width = data.get('width', 400)
            height = data.get('height', 300)
            return {
                "status": "success",
                "data": {
                    "svg": f'<svg width="{width}" height="{height}" xmlns="http://www.w3.org/2000/svg"><text x="10" y="20">Mock molecule visualization</text></svg>',
                    "width": width,
                    "height": height
                }
            }
        elif endpoint == 'substructure-search':
            return {
                "status": "success",
                "data": {
                    "match": True,
                    "match_count": 1,
                    "matches": [[2]]
                }
            }
        elif endpoint == 'similarity':
            return {
                "status": "success",
                "data": {
                    "tanimoto": 0.42,
                    "dice": 0.59,
                    "fingerprint_type": data.get('fingerprint_type', 'morgan')
                }
            }
        else:
            return {"status": "error", "message": f"Unknown endpoint: {endpoint}"}
    
    @lru_cache(maxsize=128)
    def calculate_properties(self, molecule_data: str, input_format: str = 'smiles') -> Dict[str, Any]:
        """
        Calculate properties for a molecule.
        
        Args:
            molecule_data: The molecule data (SMILES, etc.)
            input_format: Format of the input data ('smiles', 'mol', 'sdf')
            
        Returns:
            Dictionary of molecular properties.
        """
        return self._make_request('calculate-properties', {
            'molecule_data': molecule_data,
            'input_format': input_format
        })
    
    def generate_visualization(self, molecule_data: str, input_format: str = 'smiles',
                              width: int = 400, height: int = 300,
                              highlight_atoms: Optional[list] = None) -> Dict[str, Any]:
        """
        Generate a visualization of a molecule.
        
        Args:
            molecule_data: The molecule data (SMILES, etc.)
            input_format: Format of the input data ('smiles', 'mol', 'sdf')
            width: Image width in pixels
            height: Image height in pixels
            highlight_atoms: List of atom indices to highlight
            
        Returns:
            Dictionary with SVG data and dimensions.
        """
        data = {
            'molecule_data': molecule_data,
            'input_format': input_format,
            'width': width,
            'height': height
        }
        if highlight_atoms:
            data['highlight_atoms'] = highlight_atoms
            
        return self._make_request('visualization', data)
    
    def substructure_search(self, query_mol_data: str, target_mol_data: str,
                           query_format: str = 'smarts', target_format: str = 'smiles') -> Dict[str, Any]:
        """
        Perform a substructure search.
        
        Args:
            query_mol_data: The query molecule data
            target_mol_data: The target molecule data
            query_format: Format of the query data ('smarts', 'smiles')
            target_format: Format of the target data ('smiles', 'mol', 'sdf')
            
        Returns:
            Dictionary with search results.
        """
        return self._make_request('substructure-search', {
            'query_mol_data': query_mol_data,
            'target_mol_data': target_mol_data,
            'query_format': query_format,
            'target_format': target_format
        })
    
    def calculate_similarity(self, mol1_data: str, mol2_data: str,
                           mol1_format: str = 'smiles', mol2_format: str = 'smiles',
                           fingerprint_type: str = 'morgan') -> Dict[str, Any]:
        """
        Calculate similarity between two molecules.
        
        Args:
            mol1_data: The first molecule data
            mol2_data: The second molecule data
            mol1_format: Format of the first molecule data
            mol2_format: Format of the second molecule data
            fingerprint_type: Type of fingerprint to use ('morgan', 'maccs', 'topological')
            
        Returns:
            Dictionary with similarity metrics.
        """
        return self._make_request('similarity', {
            'mol1_data': mol1_data,
            'mol2_data': mol2_data,
            'mol1_format': mol1_format,
            'mol2_format': mol2_format,
            'fingerprint_type': fingerprint_type
        })