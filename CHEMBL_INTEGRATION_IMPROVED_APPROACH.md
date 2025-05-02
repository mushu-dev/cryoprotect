# ChEMBL Integration Improved Approach

## Comparison: Custom Client vs. chembl_webresource_client

### Current Approach: Custom ChEMBL Client
Our current implementation uses a custom-built ChEMBL client:
- **File**: `/mnt/c/Users/1edwa/Documents/CryoProtect v2/chembl/client.py`
- **Class**: `ResilientChEMBLClient`
- **Features**: 
  - Adaptive rate limiting
  - Multi-level caching
  - Circuit breaker pattern
  - Error handling with backoff

### Alternative Approach: chembl_webresource_client
Using the official Python client would provide several advantages:

```python
# Installation
pip install chembl_webresource_client

# Basic usage
from chembl_webresource_client.new_client import new_client
molecule = new_client.molecule
compounds = molecule.filter(pref_name__icontains="aspirin")
```

## Advantages of Using chembl_webresource_client

1. **Official Support**
   - Maintained by the ChEMBL team
   - Automatically updated for API changes
   - Compatible with all ChEMBL endpoints

2. **Simplified Implementation**
   - No need to manually handle HTTP requests
   - Automatic pagination
   - Built-in retry logic

3. **Query Language Support**
   - Django-like query syntax
   - Supports complex filtering operations
   - Chainable query methods

4. **Complete API Coverage**
   - Access to all ChEMBL resources
   - Support for all query parameters
   - Built-in data normalization

## Integration Strategy Using chembl_webresource_client

### 1. Installation
Add to `requirements.txt`:
```
chembl_webresource_client==0.10.10
```

### 2. Enhanced Implementation Example

```python
#!/usr/bin/env python3
"""
Enhanced ChEMBL data acquisition script using the official chembl_webresource_client
"""

import os
import json
import time
import logging
import argparse
from datetime import datetime
from typing import List, Dict, Any, Optional

# Import the official ChEMBL client
from chembl_webresource_client.new_client import new_client

# Import Supabase and other utils
from dotenv import load_dotenv
from service_role_helper import get_supabase_client, get_user_id, ensure_user_profile

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler("chembl_integration.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

# Load environment variables
load_dotenv()

def execute_sql_through_mcp(query: str) -> List[Dict[str, Any]]:
    """Execute SQL through MCP and return results."""
    try:
        # Use the service_role_helper to execute SQL
        # This would be implemented to use the MCP at .roo/mcp.json
        pass
    except Exception as e:
        logger.error(f"Error executing SQL: {str(e)}")
        raise


def fetch_cryoprotectant_compounds(limit: int = 1000, checkpoint_interval: int = 100) -> List[Dict[str, Any]]:
    """
    Fetch cryoprotectant compounds from ChEMBL.
    
    Args:
        limit: Maximum number of compounds to fetch
        checkpoint_interval: Save checkpoint after this many compounds
        
    Returns:
        List of compounds with properties
    """
    logger.info(f"Fetching up to {limit} cryoprotectant compounds from ChEMBL")
    
    # Initialize ChEMBL clients
    molecule = new_client.molecule
    compound_property = new_client.compound_property
    
    # Construct a query for potential cryoprotectants
    # Cryoprotectants often have specific properties:
    # - Multiple hydrogen bond donors/acceptors
    # - Moderate LogP values (water solubility)
    # - Specific molecular weight range
    
    # Start with a broad search for common cryoprotectants
    query_terms = [
        "cryoprotect",
        "glycerol",
        "dmso",
        "dimethyl sulfoxide",
        "ethylene glycol",
        "propylene glycol",
        "trehalose",
        "sucrose",
        "glucose",
        "formamide",
        "acetamide",
        "methanol",
        "polyvinyl alcohol"
    ]
    
    all_compounds = []
    checkpoint_file = "chembl_checkpoint.json"
    
    # Load checkpoint if it exists
    start_index = 0
    if os.path.exists(checkpoint_file):
        try:
            with open(checkpoint_file, 'r') as f:
                checkpoint = json.load(f)
                all_compounds = checkpoint.get("compounds", [])
                start_index = checkpoint.get("next_index", 0)
                logger.info(f"Loaded checkpoint with {len(all_compounds)} compounds, starting from index {start_index}")
        except Exception as e:
            logger.warning(f"Error loading checkpoint: {str(e)}. Starting from scratch.")
    
    # Process each query term
    for i, term in enumerate(query_terms[start_index:], start=start_index):
        try:
            logger.info(f"Searching for '{term}' ({i+1}/{len(query_terms)})")
            
            # Search for compounds matching the term
            results = molecule.filter(
                pref_name__icontains=term
            ).only(
                'molecule_chembl_id',
                'pref_name',
                'molecule_structures'
            )
            
            # Process results
            for compound in results:
                # Skip if we already have enough compounds
                if len(all_compounds) >= limit:
                    break
                    
                # Skip if no structures available
                if not compound.get('molecule_structures'):
                    continue
                    
                # Get compound details
                compound_id = compound.get('molecule_chembl_id')
                try:
                    # Get full compound details
                    full_details = molecule.get(compound_id)
                    
                    # Get compound properties
                    properties = compound_property.filter(
                        molecule_chembl_id=compound_id
                    )
                    
                    # Add properties to compound
                    full_details['properties'] = list(properties)
                    
                    # Add to results
                    all_compounds.append(full_details)
                    
                    # Log progress
                    if len(all_compounds) % 10 == 0:
                        logger.info(f"Fetched {len(all_compounds)}/{limit} compounds")
                    
                    # Save checkpoint if needed
                    if len(all_compounds) % checkpoint_interval == 0:
                        checkpoint = {
                            "compounds": all_compounds,
                            "next_index": i,
                            "timestamp": datetime.now().isoformat()
                        }
                        with open(checkpoint_file, 'w') as f:
                            json.dump(checkpoint, f)
                        logger.info(f"Saved checkpoint with {len(all_compounds)} compounds")
                        
                except Exception as e:
                    logger.warning(f"Error fetching details for {compound_id}: {str(e)}")
                    continue
                    
                # Slight delay to be gentle on the API
                time.sleep(0.2)
            
            # Break if we have enough compounds
            if len(all_compounds) >= limit:
                logger.info(f"Reached limit of {limit} compounds.")
                break
                
        except Exception as e:
            logger.error(f"Error processing term '{term}': {str(e)}")
            
            # Save checkpoint before continuing to next term
            checkpoint = {
                "compounds": all_compounds,
                "next_index": i,
                "timestamp": datetime.now().isoformat()
            }
            with open(checkpoint_file, 'w') as f:
                json.dump(checkpoint, f)
            logger.info(f"Saved checkpoint with {len(all_compounds)} compounds after error")
            
            # Continue with next term
            continue
    
    # Save final checkpoint
    checkpoint = {
        "compounds": all_compounds,
        "next_index": len(query_terms),
        "timestamp": datetime.now().isoformat(),
        "completed": True
    }
    with open(checkpoint_file, 'w') as f:
        json.dump(checkpoint, f)
    logger.info(f"Saved final checkpoint with {len(all_compounds)} compounds")
    
    return all_compounds


def transform_chembl_to_molecule(compound: Dict[str, Any], user_profile_id: str) -> Dict[str, Any]:
    """
    Transform ChEMBL compound data to match our molecule table schema.
    
    Args:
        compound: ChEMBL compound data
        user_profile_id: User profile ID for created_by field
        
    Returns:
        Dictionary matching our molecule table schema
    """
    # Extract structures
    structures = compound.get('molecule_structures', {})
    
    # Extract properties from the compound properties list
    properties_dict = {}
    for prop in compound.get('properties', []):
        prop_name = prop.get('property_name', '').lower()
        if prop_name:
            properties_dict[prop_name] = prop.get('value')
    
    # Transform to our schema
    return {
        "name": compound.get('pref_name') or compound.get('molecule_chembl_id'),
        "smiles": structures.get('canonical_smiles'),
        "inchi": structures.get('standard_inchi'),
        "inchikey": structures.get('standard_inchi_key'),
        "formula": compound.get('molecule_properties', {}).get('full_molformula'),
        "pubchem_cid": None,  # Not directly available from ChEMBL
        "molecular_weight": compound.get('molecule_properties', {}).get('full_mwt'),
        "created_by": user_profile_id,
        "data_source": f"ChEMBL ID: {compound.get('molecule_chembl_id')}",
        "version": 1,
        "modification_history": json.dumps([{
            "timestamp": datetime.now().isoformat(),
            "action": "created",
            "user_id": user_profile_id
        }])
    }


def transform_chembl_to_properties(compound: Dict[str, Any], molecule_id: str, user_profile_id: str, property_type_map: Dict[str, str]) -> List[Dict[str, Any]]:
    """
    Transform ChEMBL compound properties to match our molecular_properties table schema.
    
    Args:
        compound: ChEMBL compound data
        molecule_id: Molecule ID from our database
        user_profile_id: User profile ID for created_by field
        property_type_map: Map of property names to property type IDs
        
    Returns:
        List of dictionaries matching our molecular_properties table schema
    """
    properties = []
    
    # Process standard molecule properties
    mol_props = compound.get('molecule_properties', {})
    property_mappings = {
        'alogp': 'LogP',
        'full_mwt': 'Molecular Weight',
        'hba': 'Hydrogen Bond Acceptor Count',
        'hbd': 'Hydrogen Bond Donor Count',
        'psa': 'Topological Polar Surface Area',
        'rtb': 'Rotatable Bond Count',
        'cx_logp': 'LogP',
        'cx_logd': 'LogD',
        'aromatic_rings': 'Aromatic Ring Count',
        'heavy_atoms': 'Heavy Atom Count',
        'num_ro5_violations': 'Rule of Five Violations'
    }
    
    # Add properties from molecule_properties
    for prop_key, prop_name in property_mappings.items():
        if prop_key in mol_props and mol_props[prop_key] is not None:
            property_type_id = property_type_map.get(prop_name.lower())
            if property_type_id:
                properties.append({
                    "id": str(uuid.uuid4()),
                    "molecule_id": molecule_id,
                    "property_type_id": property_type_id,
                    "numeric_value": float(mol_props[prop_key]),
                    "text_value": None,
                    "boolean_value": None,
                    "created_by": user_profile_id,
                    "data_source": f"ChEMBL: {compound.get('molecule_chembl_id')}",
                    "version": 1,
                    "modification_history": json.dumps([{
                        "timestamp": datetime.now().isoformat(),
                        "action": "created",
                        "user_id": user_profile_id
                    }])
                })
    
    # Add properties from compound properties
    for prop in compound.get('properties', []):
        prop_name = prop.get('property_name')
        if prop_name and prop.get('value') is not None:
            property_type_id = property_type_map.get(prop_name.lower())
            if property_type_id:
                # Determine value type
                value = prop.get('value')
                numeric_value = None
                text_value = None
                boolean_value = None
                
                try:
                    numeric_value = float(value)
                except (ValueError, TypeError):
                    if isinstance(value, bool):
                        boolean_value = value
                    else:
                        text_value = str(value)
                
                properties.append({
                    "id": str(uuid.uuid4()),
                    "molecule_id": molecule_id,
                    "property_type_id": property_type_id,
                    "numeric_value": numeric_value,
                    "text_value": text_value,
                    "boolean_value": boolean_value,
                    "created_by": user_profile_id,
                    "data_source": f"ChEMBL: {compound.get('molecule_chembl_id')}",
                    "version": 1,
                    "modification_history": json.dumps([{
                        "timestamp": datetime.now().isoformat(),
                        "action": "created",
                        "user_id": user_profile_id
                    }])
                })
    
    return properties


def import_compounds_to_database(compounds: List[Dict[str, Any]], batch_size: int = 50) -> Dict[str, int]:
    """
    Import compounds to database with proper error handling and batch processing.
    
    Args:
        compounds: List of ChEMBL compounds
        batch_size: Number of compounds to process in each batch
        
    Returns:
        Dictionary with import statistics
    """
    logger.info(f"Importing {len(compounds)} compounds to database")
    
    # Get Supabase client and authenticated user
    supabase = get_supabase_client()
    auth_user_id = get_user_id()
    user_profile_id = ensure_user_profile(supabase)
    
    if not user_profile_id:
        logger.error("Failed to get user profile ID")
        return {"error": "Failed to get user profile ID"}
    
    # Get property type map
    property_type_map = {}
    try:
        property_types = execute_sql_through_mcp("""
            SELECT id, name FROM property_types;
        """)
        property_type_map = {pt['name'].lower(): pt['id'] for pt in property_types}
    except Exception as e:
        logger.error(f"Error getting property types: {str(e)}")
        return {"error": f"Failed to get property types: {str(e)}"}
    
    # Prepare statistics
    stats = {
        "total_compounds": len(compounds),
        "processed": 0,
        "molecules_inserted": 0,
        "molecules_skipped": 0,
        "properties_inserted": 0,
        "errors": 0
    }
    
    # Disable RLS temporarily
    try:
        execute_sql_through_mcp("""
            ALTER TABLE molecules DISABLE ROW LEVEL SECURITY;
            ALTER TABLE molecular_properties DISABLE ROW LEVEL SECURITY;
        """)
        logger.info("Temporarily disabled RLS for import")
    except Exception as e:
        logger.error(f"Error disabling RLS: {str(e)}")
        return {"error": f"Failed to disable RLS: {str(e)}"}
    
    try:
        # Process compounds in batches
        for i in range(0, len(compounds), batch_size):
            batch = compounds[i:i+batch_size]
            logger.info(f"Processing batch {i//batch_size + 1}/{(len(compounds)-1)//batch_size + 1} ({len(batch)} compounds)")
            
            # Start transaction
            execute_sql_through_mcp("BEGIN;")
            
            try:
                # Process each compound in the batch
                for compound in batch:
                    stats["processed"] += 1
                    
                    # Transform compound to molecule
                    molecule_data = transform_chembl_to_molecule(compound, user_profile_id)
                    
                    # Check if molecule already exists
                    inchikey = molecule_data.get('inchikey')
                    if not inchikey:
                        logger.warning(f"Skipping compound without InChIKey: {molecule_data.get('name')}")
                        stats["molecules_skipped"] += 1
                        continue
                    
                    existing = execute_sql_through_mcp(f"""
                        SELECT id FROM molecules WHERE inchikey = '{inchikey}';
                    """)
                    
                    if existing:
                        # Molecule already exists
                        molecule_id = existing[0]['id']
                        logger.debug(f"Molecule already exists: {inchikey} (ID: {molecule_id})")
                        stats["molecules_skipped"] += 1
                    else:
                        # Insert new molecule
                        result = execute_sql_through_mcp(f"""
                            INSERT INTO molecules (
                                name, smiles, inchi, inchikey, formula, 
                                molecular_weight, created_by, data_source, version, 
                                modification_history
                            ) VALUES (
                                '{molecule_data['name']}',
                                '{molecule_data['smiles'] or ''}',
                                '{molecule_data['inchi'] or ''}',
                                '{molecule_data['inchikey']}',
                                '{molecule_data['formula'] or ''}',
                                {molecule_data['molecular_weight'] or 'NULL'},
                                '{molecule_data['created_by']}',
                                '{molecule_data['data_source']}',
                                {molecule_data['version']},
                                '{molecule_data['modification_history']}'
                            )
                            RETURNING id;
                        """)
                        
                        if result:
                            molecule_id = result[0]['id']
                            logger.debug(f"Inserted molecule: {inchikey} (ID: {molecule_id})")
                            stats["molecules_inserted"] += 1
                        else:
                            logger.warning(f"Failed to insert molecule: {inchikey}")
                            stats["errors"] += 1
                            continue
                    
                    # Transform and insert properties
                    properties = transform_chembl_to_properties(
                        compound, molecule_id, user_profile_id, property_type_map
                    )
                    
                    # Insert properties
                    for prop in properties:
                        try:
                            execute_sql_through_mcp(f"""
                                INSERT INTO molecular_properties (
                                    id, molecule_id, property_type_id, 
                                    numeric_value, text_value, boolean_value,
                                    created_by, data_source, version, 
                                    modification_history
                                ) VALUES (
                                    '{prop['id']}',
                                    '{prop['molecule_id']}',
                                    '{prop['property_type_id']}',
                                    {prop['numeric_value'] if prop['numeric_value'] is not None else 'NULL'},
                                    {f"'{prop['text_value']}'" if prop['text_value'] is not None else 'NULL'},
                                    {str(prop['boolean_value']).lower() if prop['boolean_value'] is not None else 'NULL'},
                                    '{prop['created_by']}',
                                    '{prop['data_source']}',
                                    {prop['version']},
                                    '{prop['modification_history']}'
                                )
                                ON CONFLICT DO NOTHING;
                            """)
                            stats["properties_inserted"] += 1
                        except Exception as e:
                            logger.warning(f"Error inserting property for molecule {molecule_id}: {str(e)}")
                            stats["errors"] += 1
                
                # Commit transaction
                execute_sql_through_mcp("COMMIT;")
                logger.info(f"Successfully processed batch {i//batch_size + 1}")
                
            except Exception as e:
                # Rollback transaction on error
                execute_sql_through_mcp("ROLLBACK;")
                logger.error(f"Error processing batch {i//batch_size + 1}: {str(e)}")
                stats["errors"] += 1
                
            # Sleep to avoid overloading the database
            time.sleep(1)
    
    finally:
        # Re-enable RLS
        try:
            execute_sql_through_mcp("""
                ALTER TABLE molecules ENABLE ROW LEVEL SECURITY;
                ALTER TABLE molecular_properties ENABLE ROW LEVEL SECURITY;
            """)
            logger.info("Re-enabled RLS after import")
        except Exception as e:
            logger.error(f"Error re-enabling RLS: {str(e)}")
    
    return stats


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(description="Import data from ChEMBL to CryoProtect")
    parser.add_argument("--limit", type=int, default=1000, help="Maximum number of compounds to import")
    parser.add_argument("--batch-size", type=int, default=50, help="Batch size for database operations")
    parser.add_argument("--checkpoint-interval", type=int, default=100, help="Interval for saving checkpoints")
    parser.add_argument("--dry-run", action="store_true", help="Don't actually insert data, just simulate")
    args = parser.parse_args()
    
    try:
        # Fetch compounds from ChEMBL
        compounds = fetch_cryoprotectant_compounds(limit=args.limit, checkpoint_interval=args.checkpoint_interval)
        
        # Import to database if not a dry run
        if not args.dry_run:
            stats = import_compounds_to_database(compounds, batch_size=args.batch_size)
            logger.info(f"Import statistics: {json.dumps(stats, indent=2)}")
        else:
            logger.info(f"Dry run completed. {len(compounds)} compounds would be imported.")
        
        logger.info("Done.")
        
    except Exception as e:
        logger.error(f"Error in main: {str(e)}")
        return 1
    
    return 0


if __name__ == "__main__":
    sys.exit(main())
```

### 3. Benefits of This Approach

1. **Simplified Code**
   - The official client handles many low-level details
   - No need to maintain our own HTTP client
   - Reduces potential for bugs in API interaction

2. **More Robust Data Fetching**
   - Automatic handling of API pagination
   - Better error handling through official client
   - Access to additional API features

3. **Improved Maintenance**
   - Official client is maintained by ChEMBL team
   - Automatic compatibility with API updates
   - Reduced maintenance burden for our team

## Implementation Plan

### Step 1: Add Dependency
```bash
pip install chembl_webresource_client
echo "chembl_webresource_client==0.10.10" >> requirements.txt
```

### Step 2: Create New Integration Script
Create `/mnt/c/Users/1edwa/Documents/CryoProtect v2/ChEMBL_Integrated_Import.py` with the implementation above.

### Step 3: Update Roadmap
Update the `CHEMBL_INTEGRATION_ROADMAP.md` to reference this new approach.

### Step 4: Test With Small Dataset
```bash
python ChEMBL_Integrated_Import.py --limit 10 --dry-run
```

### Step 5: Execute Full Import
```bash
python ChEMBL_Integrated_Import.py --limit 1000 --batch-size 50 --checkpoint-interval 100
```

## Conclusion

The `chembl_webresource_client` approach offers significant advantages over our custom implementation. It provides a more robust, maintainable, and feature-rich interface to the ChEMBL API while reducing the complexity of our codebase. This approach will likely resolve many of the issues encountered with our current implementation and provide a more reliable path forward for ChEMBL integration.