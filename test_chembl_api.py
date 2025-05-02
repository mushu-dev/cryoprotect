#!/usr/bin/env python3
"""
Simple test script for the ChEMBL API.

This script makes a direct request to the ChEMBL API to verify connectivity
and response format.
"""

import requests
import json
import logging

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s"
)
logger = logging.getLogger(__name__)

def test_chembl_api():
    """Test direct connection to the ChEMBL API."""
    # Test URLs
    urls = [
        "https://www.ebi.ac.uk/chembl/api/data/molecule/CHEMBL25",
        "https://www.ebi.ac.uk/chembl/api/data/molecule?limit=1",
        "https://www.ebi.ac.uk/chembl/webservices/chembl_api/v1.5/molecule/CHEMBL25"
    ]
    
    for url in urls:
        logger.info(f"Testing URL: {url}")
        try:
            # Set Accept header to request JSON format
            headers = {
                "Accept": "application/json"
            }
            response = requests.get(url, headers=headers)
            logger.info(f"Status code: {response.status_code}")
            logger.info(f"Headers: {response.headers}")
            
            # Try to parse as JSON
            try:
                data = response.json()
                logger.info(f"Successfully parsed JSON response")
                logger.info(f"Response keys: {list(data.keys()) if isinstance(data, dict) else 'Not a dictionary'}")
            except json.JSONDecodeError as e:
                logger.error(f"Failed to parse JSON: {str(e)}")
                logger.info(f"Response content (first 500 chars): {response.text[:500]}")
        
        except requests.RequestException as e:
            logger.error(f"Request failed: {str(e)}")

if __name__ == "__main__":
    test_chembl_api()