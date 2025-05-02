"""
CryoProtect Analyzer - Mock Supabase RPC

This module provides mock implementations of Supabase RPC functions.
"""

import uuid
from datetime import datetime
from .data import _mock_data, _mock_rpcs, DEFAULT_USER_ID
from .query import MockResponse

# Mock RPC
class MockRPC:
    """Mock RPC (Remote Procedure Call) handler."""
    
    def __init__(self, function_name, params=None):
        self.function_name = function_name
        self.params = params or {}
    
    def execute(self):
        """Execute the RPC function."""
        try:
            # Check if the function is registered
            if self.function_name in _mock_rpcs:
                # Call the function with parameters
                result = _mock_rpcs[self.function_name](self.params)
                return MockResponse(result)
            
            # Handle known RPCs with default implementations
            if self.function_name == 'import_molecule_from_pubchem':
                return self._import_molecule_from_pubchem()
            elif self.function_name == 'get_shared_item_data':
                return self._get_shared_item_data()
            elif self.function_name == 'calculate_mixture_score':
                return self._calculate_mixture_score()
            elif self.function_name == 'get_mixture_with_components':
                return self._get_mixture_with_components()
            
            # Unknown function
            return MockResponse(
                None, 
                {"message": f"RPC function '{self.function_name}' not implemented in mock"}
            )
        
        except Exception as e:
            return MockResponse(None, {"message": str(e)})
    
    def _import_molecule_from_pubchem(self):
        """Mock implementation of import_molecule_from_pubchem RPC."""
        cid = self.params.get('p_cid')
        user_id = self.params.get('p_user_id', DEFAULT_USER_ID)
        
        # Create a new molecule with PubChem data
        molecule_id = str(uuid.uuid4())
        
        # Simple mock data based on CID
        molecule = {
            'id': molecule_id,
            'cid': cid,
            'name': f'PubChem Compound {cid}',
            'molecular_formula': 'C10H20O5',  # Example formula
            'smiles': 'CC(C)CCOC(=O)C',  # Example SMILES
            'created_at': datetime.now().isoformat(),
            'updated_at': datetime.now().isoformat(),
            'created_by': user_id
        }
        
        _mock_data['molecules'].append(molecule)
        
        # Add some properties
        property_types = _mock_data['property_types']
        if property_types:
            mw_property_id = property_types[0]['id']
            logp_property_id = property_types[1]['id'] if len(property_types) > 1 else None
            
            if mw_property_id:
                _mock_data['molecular_properties'].append({
                    'id': str(uuid.uuid4()),
                    'molecule_id': molecule_id,
                    'property_type_id': mw_property_id,
                    'numeric_value': 188.26,  # Example value
                    'created_at': datetime.now().isoformat(),
                    'updated_at': datetime.now().isoformat(),
                    'created_by': user_id
                })
            
            if logp_property_id:
                _mock_data['molecular_properties'].append({
                    'id': str(uuid.uuid4()),
                    'molecule_id': molecule_id,
                    'property_type_id': logp_property_id,
                    'numeric_value': 2.45,  # Example value
                    'created_at': datetime.now().isoformat(),
                    'updated_at': datetime.now().isoformat(),
                    'created_by': user_id
                })
        
        return MockResponse(molecule_id)
    
    def _get_shared_item_data(self):
        """Mock implementation of get_shared_item_data RPC."""
        share_id = self.params.get('p_share_id')
        
        # Find the share
        share = next((s for s in _mock_data['shares'] if s['id'] == share_id), None)
        if not share:
            return MockResponse(None, {"message": "Share not found"})
        
        # Get the item based on data_type
        data_type = share.get('data_type')
        item_id = share.get('item_id')
        
        if data_type == 'molecules':
            item = next((m for m in _mock_data['molecules'] if m['id'] == item_id), None)
        elif data_type == 'mixtures':
            item = next((m for m in _mock_data['mixtures'] if m['id'] == item_id), None)
        else:
            item = None
        
        if not item:
            return MockResponse(None, {"message": "Shared item not found"})
        
        return MockResponse(item)
    
    def _calculate_mixture_score(self):
        """Mock implementation of calculate_mixture_score RPC."""
        mixture_id = self.params.get('p_mixture_id')
        
        # Find the mixture
        mixture = next((m for m in _mock_data['mixtures'] if m['id'] == mixture_id), None)
        if not mixture:
            return MockResponse(None, {"message": "Mixture not found"})
        
        # Calculate a mock score (in a real implementation, this would use actual algorithms)
        # For testing, we'll just return a random score between 0 and 100
        import random
        score = round(random.uniform(60, 95), 1)  # Most cryoprotectants should score well
        
        return MockResponse(score)
    
    def _get_mixture_with_components(self):
        """Mock implementation of get_mixture_with_components RPC."""
        mixture_id = self.params.get('p_mixture_id')
        
        # Find the mixture
        mixture = next((m for m in _mock_data['mixtures'] if m['id'] == mixture_id), None)
        if not mixture:
            return MockResponse(None, {"message": "Mixture not found"})
        
        # Find components
        components = [c for c in _mock_data['mixture_components'] if c['mixture_id'] == mixture_id]
        
        # Enrich components with molecule data
        enriched_components = []
        for component in components:
            molecule_id = component['molecule_id']
            molecule = next((m for m in _mock_data['molecules'] if m['id'] == molecule_id), None)
            
            if molecule:
                enriched_component = {
                    **component,
                    'molecule_name': molecule['name'],
                    'molecule_formula': molecule['molecular_formula'],
                    'molecule_smiles': molecule['smiles']
                }
                enriched_components.append(enriched_component)
        
        # Create the result
        result = {
            **mixture,
            'components': enriched_components
        }
        
        return MockResponse(result)


def register_rpc_function(name, func):
    """
    Register a custom RPC function.
    
    Args:
        name: Name of the RPC function
        func: Function to call when the RPC is executed
    """
    _mock_rpcs[name] = func