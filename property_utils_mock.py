#!/usr/bin/env python3
"""
Simplified mock version of property_utils.py for testing purposes.
"""

import logging
import uuid
from typing import Dict, Any, Optional, List, Tuple
from datetime import datetime

from db_connection_utils_mock import get_db_connection
from transaction_utils_mock import safe_transaction, with_transaction_retry, execute_in_transaction

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger('property_utils_mock')

class PropertyManager:
    """
    Mock implementation of PropertyManager for testing.
    """
    
    def __init__(self, connection=None):
        """
        Initialize the PropertyManager.
        """
        self.db = connection or get_db_connection()
        self._property_types_cache = {}
        self._properties_store = {}  # Mock in-memory store for properties
        
        # Pre-populate cache with critical property types
        self._property_types_cache = {
            'logP': {'id': str(uuid.uuid4()), 'data_type': 'numeric'},
            'h_bond_donors': {'id': str(uuid.uuid4()), 'data_type': 'numeric'},
            'h_bond_acceptors': {'id': str(uuid.uuid4()), 'data_type': 'numeric'},
            'alogp': {'id': str(uuid.uuid4()), 'data_type': 'numeric'},
            'molecular_weight': {'id': str(uuid.uuid4()), 'data_type': 'numeric'},
            'rotatable_bonds': {'id': str(uuid.uuid4()), 'data_type': 'numeric'},
            'psa': {'id': str(uuid.uuid4()), 'data_type': 'numeric'},
            'heavy_atoms': {'id': str(uuid.uuid4()), 'data_type': 'numeric'}
        }
        
        logger.info(f"Mock PropertyManager initialized with {len(self._property_types_cache)} property types")
    
    def get_property_type_id(self, property_name: str, data_type: str = 'numeric') -> str:
        """
        Get the property type ID for a given property name.
        
        Args:
            property_name: Name of the property
            data_type: Data type of the property ('numeric', 'text', 'boolean')
            
        Returns:
            UUID of the property type
        """
        if property_name in self._property_types_cache:
            logger.info(f"Found property type {property_name} in cache")
            return self._property_types_cache[property_name]['id']
        
        # Create new property type
        property_id = str(uuid.uuid4())
        self._property_types_cache[property_name] = {
            'id': property_id,
            'data_type': data_type
        }
        
        logger.info(f"Created new property type: {property_name} ({data_type}) with ID {property_id}")
        return property_id
    
    def set_property(self, molecule_id: str, property_name: str,
                    property_value: Any, created_by: Optional[str] = None) -> bool:
        """
        Set a property value for a molecule.
        
        Args:
            molecule_id: UUID of the molecule
            property_name: Name of the property
            property_value: Value of the property
            created_by: UUID of the user creating the property
            
        Returns:
            True if successful, False otherwise
        """
        if property_value is None:
            return True
        
        try:
            # Get property type ID
            property_type_id = self.get_property_type_id(property_name)
            
            # Store property in mock store
            if molecule_id not in self._properties_store:
                self._properties_store[molecule_id] = {}
            
            self._properties_store[molecule_id][property_name] = {
                'value': property_value,
                'property_type_id': property_type_id,
                'created_at': datetime.now().isoformat(),
                'created_by': created_by
            }
            
            logger.info(f"Set property {property_name}={property_value} for molecule {molecule_id}")
            return True
            
        except Exception as e:
            logger.error(f"Error setting property {property_name}={property_value} for molecule {molecule_id}: {str(e)}")
            return False
    
    def set_properties(self, molecule_id: str, properties: Dict[str, Any],
                      created_by: Optional[str] = None) -> Tuple[int, int]:
        """
        Set multiple properties for a molecule.
        
        Args:
            molecule_id: UUID of the molecule
            properties: Dictionary of property name -> value
            created_by: UUID of the user creating the properties
            
        Returns:
            Tuple of (success_count, total_count)
        """
        if not properties:
            return 0, 0
        
        success_count = 0
        total_count = len(properties)
        
        for property_name, property_value in properties.items():
            if self.set_property(molecule_id, property_name, property_value, created_by):
                success_count += 1
        
        logger.info(f"Set {success_count}/{total_count} properties for molecule {molecule_id}")
        return success_count, total_count
    
    def get_properties(self, molecule_id: str, property_names: Optional[List[str]] = None) -> Dict[str, Any]:
        """
        Get properties for a molecule.
        
        Args:
            molecule_id: UUID of the molecule
            property_names: Optional list of property names to retrieve
            
        Returns:
            Dictionary of property name -> value
        """
        result = {}
        
        if molecule_id not in self._properties_store:
            return result
        
        if property_names:
            # Get specific properties
            for prop_name in property_names:
                if prop_name in self._properties_store[molecule_id]:
                    result[prop_name] = self._properties_store[molecule_id][prop_name]['value']
        else:
            # Get all properties
            for prop_name, prop_data in self._properties_store[molecule_id].items():
                result[prop_name] = prop_data['value']
        
        logger.info(f"Retrieved {len(result)} properties for molecule {molecule_id}")
        return result
    
    def batch_get_properties(self, molecule_ids: List[str],
                           property_names: Optional[List[str]] = None) -> Dict[str, Dict[str, Any]]:
        """
        Get properties for multiple molecules.
        
        Args:
            molecule_ids: List of molecule UUIDs
            property_names: Optional list of property names to retrieve
            
        Returns:
            Dictionary of molecule_id -> {property_name: value}
        """
        result = {}
        
        for molecule_id in molecule_ids:
            result[molecule_id] = self.get_properties(molecule_id, property_names)
        
        return result
    
    def delete_property(self, molecule_id: str, property_name: str) -> bool:
        """
        Delete a property for a molecule.
        
        Args:
            molecule_id: UUID of the molecule
            property_name: Name of the property
            
        Returns:
            True if successful, False otherwise
        """
        if molecule_id not in self._properties_store:
            return False
        
        if property_name not in self._properties_store[molecule_id]:
            return False
        
        del self._properties_store[molecule_id][property_name]
        logger.info(f"Deleted property {property_name} for molecule {molecule_id}")
        return True