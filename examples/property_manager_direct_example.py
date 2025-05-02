#!/usr/bin/env python3
"""
Example script demonstrating the usage of the enhanced PropertyManager
with direct PostgreSQL connections.

This script shows how to:
1. Initialize the PropertyManager
2. Set individual properties
3. Set multiple properties at once
4. Batch set properties for multiple molecules
5. Retrieve properties
6. Batch retrieve properties for multiple molecules
7. Delete properties
"""

import os
import uuid
import logging
from dotenv import load_dotenv

# Import the enhanced PropertyManager
from property_utils import PropertyManager

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger('property_manager_example')

def main():
    """Main function demonstrating PropertyManager usage."""
    # Load environment variables for database connection
    load_dotenv()
    
    # Verify required environment variables are set
    required_vars = ['SUPABASE_DB_HOST', 'SUPABASE_DB_PASSWORD']
    missing_vars = [var for var in required_vars if not os.getenv(var)]
    if missing_vars:
        logger.error(f"Missing required environment variables: {', '.join(missing_vars)}")
        logger.error("Please set these variables in your .env file or environment")
        return
    
    logger.info("Initializing PropertyManager with direct PostgreSQL connection")
    property_manager = PropertyManager()
    
    # Create some test molecule IDs (in a real scenario, these would be actual molecule IDs from the database)
    molecule_id1 = uuid.uuid4()
    molecule_id2 = uuid.uuid4()
    molecule_id3 = uuid.uuid4()
    
    logger.info(f"Created test molecule IDs: {molecule_id1}, {molecule_id2}, {molecule_id3}")
    
    # Example 1: Set individual properties
    logger.info("Example 1: Setting individual properties")
    
    # Set a numeric property
    success = property_manager.set_property(
        molecule_id1, 
        'logP', 
        2.45
    )
    logger.info(f"Set logP for molecule {molecule_id1}: {'Success' if success else 'Failed'}")
    
    # Set a text property
    success = property_manager.set_property(
        molecule_id1, 
        'smiles', 
        'C1=CC=CC=C1'
    )
    logger.info(f"Set SMILES for molecule {molecule_id1}: {'Success' if success else 'Failed'}")
    
    # Set a boolean property
    success = property_manager.set_property(
        molecule_id1, 
        'is_cryoprotectant', 
        True
    )
    logger.info(f"Set is_cryoprotectant for molecule {molecule_id1}: {'Success' if success else 'Failed'}")
    
    # Example 2: Set multiple properties at once
    logger.info("Example 2: Setting multiple properties at once")
    
    properties = {
        'logP': 3.21,
        'molecular_weight': 122.16,
        'smiles': 'C1=CC=CC=C1O',
        'h_bond_donors': 1,
        'h_bond_acceptors': 1,
        'is_cryoprotectant': False
    }
    
    success_count, total_count = property_manager.set_properties(
        molecule_id2, 
        properties
    )
    
    logger.info(f"Set {success_count}/{total_count} properties for molecule {molecule_id2}")
    
    # Example 3: Batch set properties for multiple molecules
    logger.info("Example 3: Batch setting properties for multiple molecules")
    
    batch_data = [
        {
            'molecule_id': molecule_id1,
            'properties': {
                'rotatable_bonds': 0,
                'tpsa': 0.0
            }
        },
        {
            'molecule_id': molecule_id2,
            'properties': {
                'rotatable_bonds': 1,
                'tpsa': 20.23
            }
        },
        {
            'molecule_id': molecule_id3,
            'properties': {
                'logP': 1.97,
                'molecular_weight': 92.09,
                'smiles': 'C1=CC=CC=C1N',
                'h_bond_donors': 1,
                'h_bond_acceptors': 1,
                'rotatable_bonds': 0,
                'tpsa': 15.79,
                'is_cryoprotectant': False
            }
        }
    ]
    
    success_count, total_count = property_manager.batch_set_properties(batch_data)
    logger.info(f"Batch set {success_count}/{total_count} properties for {len(batch_data)} molecules")
    
    # Example 4: Retrieve properties for a single molecule
    logger.info("Example 4: Retrieving properties for a single molecule")
    
    # Get all properties
    properties = property_manager.get_properties(molecule_id1)
    logger.info(f"Properties for molecule {molecule_id1}:")
    for name, value in properties.items():
        logger.info(f"  {name}: {value}")
    
    # Get specific properties
    specific_properties = property_manager.get_properties(
        molecule_id2, 
        ['logP', 'molecular_weight', 'smiles']
    )
    
    logger.info(f"Specific properties for molecule {molecule_id2}:")
    for name, value in specific_properties.items():
        logger.info(f"  {name}: {value}")
    
    # Example 5: Batch retrieve properties for multiple molecules
    logger.info("Example 5: Batch retrieving properties for multiple molecules")
    
    batch_properties = property_manager.batch_get_properties(
        [molecule_id1, molecule_id2, molecule_id3],
        ['logP', 'smiles', 'is_cryoprotectant']
    )
    
    for molecule_id, props in batch_properties.items():
        logger.info(f"Properties for molecule {molecule_id}:")
        for name, value in props.items():
            logger.info(f"  {name}: {value}")
    
    # Example 6: Delete a property
    logger.info("Example 6: Deleting a property")
    
    success = property_manager.delete_property(molecule_id1, 'tpsa')
    logger.info(f"Delete TPSA for molecule {molecule_id1}: {'Success' if success else 'Failed'}")
    
    # Example 7: Clear the property types cache
    logger.info("Example 7: Clearing the property types cache")
    property_manager.clear_cache()
    logger.info("Property types cache cleared")
    
    logger.info("PropertyManager example completed successfully")

if __name__ == "__main__":
    main()