"""
Simple PubChem API client for CryoProtect v2.

This module provides a simplified client for interacting with the PubChem API,
designed to work with the enhanced PubChem importer.
"""

import logging
import requests
from typing import Dict, Any, List, Optional

logger = logging.getLogger(__name__)

class PubChemClient:
    """
    Simple PubChem API client for the enhanced importer.
    
    This is a simplified version of ResilientPubChemClient for use with
    the enhanced PubChem importer, which handles its own caching and rate limiting.
    """
    
    def __init__(self):
        """Initialize the client."""
        self.base_url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
        
    def get_compound(self, cid: str) -> Dict[str, Any]:
        """
        Get compound data from PubChem.
        
        Args:
            cid: PubChem CID
            
        Returns:
            Dictionary containing compound data
            
        Raises:
            Exception: If the request fails
        """
        url = f"{self.base_url}/compound/cid/{cid}/JSON"
        
        try:
            response = requests.get(url)
            response.raise_for_status()
            data = response.json()
            
            # Extract and format compound data
            compound = data.get('PC_Compounds', [])[0] if 'PC_Compounds' in data else {}
            
            # Get basic information
            result = {
                'id': compound.get('id', {}),
                'name': self._get_compound_name(cid),
                'synonyms': self._get_compound_synonyms(cid),
                'molecular_formula': '',
                'molecular_weight': 0.0,
                'inchi': '',
                'inchi_key': '',
                'smiles': '',
                'props': []
            }
            
            # Extract properties
            if 'props' in compound:
                result['props'] = compound['props']
                
                for prop in compound['props']:
                    if prop.get('urn', {}).get('label') == 'IUPAC Name':
                        result['name'] = prop.get('value', {}).get('sval', '')
                    elif prop.get('urn', {}).get('label') == 'InChI':
                        result['inchi'] = prop.get('value', {}).get('sval', '')
                    elif prop.get('urn', {}).get('label') == 'InChIKey':
                        result['inchi_key'] = prop.get('value', {}).get('sval', '')
                    elif prop.get('urn', {}).get('label') == 'SMILES':
                        result['smiles'] = prop.get('value', {}).get('sval', '')
                    elif prop.get('urn', {}).get('label') == 'Molecular Formula':
                        result['molecular_formula'] = prop.get('value', {}).get('sval', '')
                    elif prop.get('urn', {}).get('label') == 'Molecular Weight':
                        result['molecular_weight'] = prop.get('value', {}).get('fval', 0.0)
            
            return result
            
        except Exception as e:
            logger.error(f"Error fetching compound {cid}: {str(e)}")
            raise
    
    def _get_compound_name(self, cid: str) -> str:
        """
        Get compound name from PubChem.
        
        Args:
            cid: PubChem CID
            
        Returns:
            Compound name
        """
        url = f"{self.base_url}/compound/cid/{cid}/description/JSON"
        
        try:
            response = requests.get(url)
            response.raise_for_status()
            data = response.json()
            
            if 'InformationList' in data and 'Information' in data['InformationList']:
                info = data['InformationList']['Information'][0]
                if 'Title' in info:
                    return info['Title']
            
            return ""
            
        except Exception as e:
            logger.debug(f"Error fetching compound name for {cid}: {str(e)}")
            return ""
    
    def _get_compound_synonyms(self, cid: str) -> List[str]:
        """
        Get compound synonyms from PubChem.
        
        Args:
            cid: PubChem CID
            
        Returns:
            List of synonyms
        """
        url = f"{self.base_url}/compound/cid/{cid}/synonyms/JSON"
        
        try:
            response = requests.get(url)
            response.raise_for_status()
            data = response.json()
            
            if 'InformationList' in data and 'Information' in data['InformationList']:
                info = data['InformationList']['Information'][0]
                if 'Synonym' in info:
                    return info['Synonym']
            
            return []
            
        except Exception as e:
            logger.debug(f"Error fetching synonyms for {cid}: {str(e)}")
            return []