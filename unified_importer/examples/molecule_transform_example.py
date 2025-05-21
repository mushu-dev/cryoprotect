#!/usr/bin/env python3
"""
Example of using the molecule transform module.

This script demonstrates how to use the MoleculeTransformer class
for standardizing molecular data, calculating properties, and
handling cross-references between databases.
"""

import asyncio
import json
import logging
from typing import Dict, Any, List

from unified_importer.transforms.molecule_transform import MoleculeTransformer


# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger('molecule_transform_example')


# Sample molecules for demonstration
EXAMPLE_MOLECULES = [
    {
        'name': 'Aspirin',
        'smiles': 'CC(=O)OC1=CC=CC=C1C(=O)O',
        'chembl_id': 'CHEMBL25',
        'data_source': 'ChEMBL'
    },
    {
        'name': 'Glucose',
        'smiles': 'C(C1C(C(C(C(O1)O)O)O)O)O',
        'chembl_id': 'CHEMBL444',
        'data_source': 'ChEMBL'
    },
    {
        'name': 'DMSO (Dimethyl sulfoxide)',
        'smiles': 'CS(=O)C',
        'pubchem_cid': '679',
        'data_source': 'PubChem'
    },
    {
        'name': 'Sodium Chloride Solution',
        'smiles': '[Na+].[Cl-].[H]O[H]',
        'data_source': 'Example',
        'is_mixture': True
    }
]


async def demonstrate_standardization():
    """Demonstrate basic molecule standardization."""
    logger.info("DEMONSTRATING BASIC STANDARDIZATION")
    
    # Create transformer with default settings
    transformer = MoleculeTransformer(logger=logger)
    
    # Process a simple molecule
    molecule = EXAMPLE_MOLECULES[0].copy()
    result = await transformer.standardize_molecule(molecule, resolve_ids=False)
    
    # Display results
    logger.info(f"Standardized molecule: {result['name']}")
    logger.info(f"SMILES: {result['smiles']}")
    logger.info(f"InChI: {result['inchi']}")
    logger.info(f"Formula: {result['formula']}")
    logger.info(f"Molecular Weight: {result['molecular_weight']:.2f}")
    
    # Display calculated properties
    logger.info("Calculated properties:")
    for prop_name, value in result['properties'].items():
        if isinstance(value, (int, float)):
            logger.info(f"  {prop_name}: {value:.2f}")
        else:
            logger.info(f"  {prop_name}: {value}")
    
    return result


async def demonstrate_batch_processing():
    """Demonstrate batch processing of multiple molecules."""
    logger.info("\nDEMONSTRATING BATCH PROCESSING")
    
    # Create transformer with custom batch settings
    config = {
        'batch_size': 10,
        'resolve_cross_references': False
    }
    transformer = MoleculeTransformer(logger=logger, config=config)
    
    # Process all example molecules as a batch
    results = await transformer.standardize_molecules_batch(
        [mol.copy() for mol in EXAMPLE_MOLECULES], 
        resolve_ids=False
    )
    
    # Display summary
    logger.info(f"Processed {len(results)} molecules in batch")
    for i, result in enumerate(results):
        logger.info(f"  {i+1}. {result['name']} - {result.get('formula', 'N/A')}")
    
    return results


async def demonstrate_mixture_handling():
    """Demonstrate handling of molecular mixtures."""
    logger.info("\nDEMONSTRATING MIXTURE HANDLING")
    
    # Create transformer with mixture handling enabled
    config = {'handle_mixtures': True}
    transformer = MoleculeTransformer(logger=logger, config=config)
    
    # Process a mixture
    mixture = EXAMPLE_MOLECULES[3].copy()
    result = await transformer.standardize_molecule(mixture, resolve_ids=False)
    
    # Display mixture information
    logger.info(f"Mixture: {result['name']}")
    logger.info(f"SMILES: {result['smiles']}")
    logger.info(f"Is mixture: {result['is_mixture']}")
    
    # Display component information
    if result['is_mixture'] and 'mixture_components' in result:
        logger.info(f"Number of components: {len(result['mixture_components'])}")
        for i, component in enumerate(result['mixture_components']):
            logger.info(f"Component {i+1}: {component['name']}")
            logger.info(f"  Formula: {component['formula']}")
            logger.info(f"  SMILES: {component['smiles']}")
            logger.info(f"  Molecular Weight: {component['molecular_weight']:.2f}")
    
    return result


async def demonstrate_cross_reference_resolution():
    """Demonstrate resolution of cross-database references."""
    logger.info("\nDEMONSTRATING CROSS-REFERENCE RESOLUTION")
    logger.info("Note: This requires internet connection and may take a moment...")
    
    # Create transformer with cross-reference resolution enabled
    config = {'resolve_cross_references': True}
    transformer = MoleculeTransformer(logger=logger, config=config)
    
    # Process a ChEMBL molecule
    chembl_mol = EXAMPLE_MOLECULES[0].copy()
    result_chembl = await transformer.standardize_molecule(chembl_mol)
    
    # Process a PubChem molecule
    pubchem_mol = EXAMPLE_MOLECULES[2].copy()
    result_pubchem = await transformer.standardize_molecule(pubchem_mol)
    
    # Display results
    logger.info("ChEMBL -> PubChem:")
    logger.info(f"  ChEMBL ID: {result_chembl['chembl_id']}")
    logger.info(f"  Resolved PubChem CID: {result_chembl.get('pubchem_cid', 'Not resolved')}")
    
    logger.info("PubChem -> ChEMBL:")
    logger.info(f"  PubChem CID: {result_pubchem['pubchem_cid']}")
    logger.info(f"  Resolved ChEMBL ID: {result_pubchem.get('chembl_id', 'Not resolved')}")
    
    return result_chembl, result_pubchem


async def demonstrate_merge_operation():
    """Demonstrate merging of molecule data from different sources."""
    logger.info("\nDEMONSTRATING MERGE OPERATION")
    
    # Create two overlapping molecule records
    chembl_data = {
        'name': 'Aspirin',
        'smiles': 'CC(=O)OC1=CC=CC=C1C(=O)O',
        'chembl_id': 'CHEMBL25',
        'molecular_weight': 180.16,
        'data_source': 'ChEMBL',
        'synonyms': ['Aspirin', 'Acetylsalicylic acid', 'ASA']
    }
    
    pubchem_data = {
        'name': 'Acetylsalicylic acid',
        'smiles': 'CC(=O)OC1=CC=CC=C1C(=O)O',
        'pubchem_cid': '2244',
        'formula': 'C9H8O4',
        'data_source': 'PubChem',
        'synonyms': ['Aspirin', '2-acetoxybenzoic acid', 'Acetylsalicylic acid']
    }
    
    # Create transformer
    transformer = MoleculeTransformer(logger=logger)
    
    # Merge the records, preferring ChEMBL data
    merged = await transformer.merge_molecule_data(chembl_data, pubchem_data, prefer_source='ChEMBL')
    
    # Display merged result
    logger.info("Merged molecule record:")
    logger.info(f"  Name: {merged['name']}")
    logger.info(f"  ChEMBL ID: {merged.get('chembl_id', 'None')}")
    logger.info(f"  PubChem CID: {merged.get('pubchem_cid', 'None')}")
    logger.info(f"  Formula: {merged.get('formula', 'None')}")
    logger.info(f"  Synonyms: {', '.join(merged.get('synonyms', []))}")
    
    return merged


async def main():
    """Run all demonstrations."""
    logger.info("MOLECULE TRANSFORMATION EXAMPLES")
    logger.info("================================")
    
    # Run each demonstration
    await demonstrate_standardization()
    await demonstrate_batch_processing()
    await demonstrate_mixture_handling()
    
    # Skip cross-reference resolution by default (requires internet)
    # Uncomment to run this demonstration
    # await demonstrate_cross_reference_resolution()
    
    await demonstrate_merge_operation()
    
    logger.info("\nAll demonstrations completed successfully")


if __name__ == "__main__":
    asyncio.run(main())