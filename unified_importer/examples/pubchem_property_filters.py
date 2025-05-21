#!/usr/bin/env python3
"""
Example script demonstrating the use of PubChem property filters.

This script shows how to configure and use property filters with
the PubChemDataSource to retrieve compounds matching specific criteria.
"""

import os
import json
import logging
import asyncio
from typing import Dict, List, Any, Optional

# Import the unified importer modules
from unified_importer.sources.pubchem_source import PubChemDataSource
from unified_importer.core.config import load_config
from unified_importer.core.database import DatabaseOperations
from unified_importer.core.checkpoint import CheckpointManager
from unified_importer.core.progress import ProgressTracker


# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger('pubchem_property_filters')


class MockDatabase(DatabaseOperations):
    """Mock database operations class for demonstration."""
    
    async def insert_compound(self, compound_data: Dict[str, Any]) -> str:
        """Mock insert compound."""
        logger.info(f"Would insert compound: {compound_data.get('name', 'Unknown')}")
        return "mock-id"
        
    async def insert_properties(self, compound_id: str, properties: List[Dict[str, Any]]) -> List[str]:
        """Mock insert properties."""
        logger.info(f"Would insert {len(properties)} properties for compound {compound_id}")
        return ["mock-prop-id"] * len(properties)
        
    async def insert_batch(self, compounds: List[Dict[str, Any]]) -> List[str]:
        """Mock batch insert."""
        logger.info(f"Would batch insert {len(compounds)} compounds")
        return ["mock-id"] * len(compounds)


async def main():
    """Main function demonstrating property filters."""
    # Create example config with property filters
    config = {
        "database": {
            "url": "mock-url",
            "key": "mock-key",
            "use_connection_pool": True
        },
        "logging": {
            "level": "INFO",
            "format": "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
        },
        "checkpoints": {
            "directory": "checkpoints",
            "enabled": False
        },
        "sources": {
            "pubchem": {
                "api_url": "https://pubchem.ncbi.nlm.nih.gov/rest/pug",
                "api_delay": 1.0,
                "timeout": 30,
                "property_filters": [
                    {
                        "name": "cryoprotectants", 
                        "description": "Known cryoprotectants and related compounds",
                        "terms": ["glycerol", "dmso", "ethylene glycol", "trehalose", "sucrose", "cryoprotect"]
                    },
                    {
                        "name": "small_molecules",
                        "description": "Small molecules suitable for cryoprotection",
                        "molecular_weight_max": 200,
                        "logp_max": 1.0,
                        "rotatable_bonds_max": 5
                    },
                    {
                        "name": "hydrogen_bond_rich",
                        "description": "Molecules with many hydrogen bond donors/acceptors",
                        "hbond_donors_min": 3,
                        "hbond_acceptors_min": 3
                    }
                ]
            }
        }
    }
    
    # Save config to a file for demonstration
    os.makedirs('examples', exist_ok=True)
    with open('examples/pubchem_filters_config.json', 'w') as f:
        json.dump(config, f, indent=4)
    
    # Initialize components
    db = MockDatabase()
    checkpoint_manager = CheckpointManager()
    progress_tracker = ProgressTracker()
    
    # Create the PubChem data source
    pubchem = PubChemDataSource(
        db_operations=db,
        checkpoint_manager=checkpoint_manager,
        progress_tracker=progress_tracker,
        config=config,
        logger=logger
    )
    
    # ------------- Example 1: Search by Filter -------------
    logger.info("=== Example 1: Search by Filter ===")
    
    # Search for compounds matching the 'cryoprotectants' filter
    logger.info("Searching for compounds matching the 'cryoprotectants' filter...")
    cryoprotectants = await pubchem.search_compounds_by_filter(
        'cryoprotectants', 
        max_results=5
    )
    
    logger.info(f"Found {len(cryoprotectants)} matching compounds: {cryoprotectants}")
    
    # Fetch details for the first compound
    if cryoprotectants:
        compound_data = await pubchem.fetch_compound(cryoprotectants[0])
        transformed = await pubchem.transform_compound_data(compound_data)
        properties = await pubchem.get_property_data(compound_data)
        
        logger.info(f"Compound details: {transformed['name']}")
        logger.info(f"Properties: {len(properties)} properties found")
    
    # ------------- Example 2: Apply Filters to Search Query -------------
    logger.info("\n=== Example 2: Apply Filters to Search Query ===")
    
    # Search for compounds matching 'alcohol' with filters applied
    logger.info("Searching for 'alcohol' with filters applied...")
    alcohols_filtered = await pubchem.search_compounds(
        'alcohol', 
        max_results=5, 
        apply_filters=True
    )
    
    logger.info(f"Found {len(alcohols_filtered)} matching compounds with filters: {alcohols_filtered}")
    
    # Search without filters for comparison
    logger.info("Searching for 'alcohol' without filters...")
    alcohols_unfiltered = await pubchem.search_compounds(
        'alcohol', 
        max_results=5, 
        apply_filters=False
    )
    
    logger.info(f"Found {len(alcohols_unfiltered)} matching compounds without filters: {alcohols_unfiltered}")
    
    # ------------- Example 3: Stream Compounds with Filters -------------
    logger.info("\n=== Example 3: Stream Compounds with Filters ===")
    
    # Stream compounds matching 'sugar' with filters
    logger.info("Streaming compounds matching 'sugar' with filters...")
    compound_count = 0
    async for batch in pubchem.stream_compound_identifiers('sugar', batch_size=5, apply_filters=True):
        compound_count += len(batch)
        logger.info(f"Got batch of {len(batch)} compounds: {batch}")
        
        # Stop after a few batches for demonstration
        if compound_count >= 15:
            break
    
    logger.info(f"Streamed a total of {compound_count} compounds")
    
    # Clean up
    await pubchem._close_session()
    logger.info("Done!")


if __name__ == "__main__":
    asyncio.run(main())