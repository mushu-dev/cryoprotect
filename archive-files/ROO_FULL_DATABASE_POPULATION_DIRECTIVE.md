# ROO FULL DATABASE POPULATION DIRECTIVE

## TASK OVERVIEW
Complete the population of the CryoProtect database with comprehensive scientific data by integrating data from multiple sources, enhancing existing data quality, and ensuring proper cross-referencing between different identifier systems.

## SUCCESS CRITERIA
- ChEMBL data for at least 500 cryoprotectant compounds is imported with complete properties
- All imported compounds have comprehensive property data including logP, H-bond donors/acceptors, etc.
- All reference compounds are present with complete property information
- Cross-references between PubChem CIDs and ChEMBL IDs are established for all possible compounds
- Molecular property visualization works for all imported molecules
- Database queries perform within acceptable parameters (<100ms average response time)

## FILE REFERENCES

### Key Implementation Files:
- `/mnt/c/Users/1edwa/Documents/CryoProtect v2/import_full_chembl.py`
- `/mnt/c/Users/1edwa/Documents/CryoProtect v2/chembl/reference_compounds.py`
- `/mnt/c/Users/1edwa/Documents/CryoProtect v2/PubChem_CryoProtectants_Supabase_Enhanced.py`
- `/mnt/c/Users/1edwa/Documents/CryoProtect v2/cryoprotectant_identifiers.py`
- `/mnt/c/Users/1edwa/Documents/CryoProtect v2/reconcile_chembl_properties.py`

### Key Reference Files:
- `/mnt/c/Users/1edwa/Documents/CryoProtect v2/reports/database_quality_report.md/database_quality_report.md`
- `/mnt/c/Users/1edwa/Documents/CryoProtect v2/reports/reference_compounds_verification.md`
- `/mnt/c/Users/1edwa/Documents/CryoProtect v2/reports/pubchem_data_quality.md`
- `/mnt/c/Users/1edwa/Documents/CryoProtect v2/reports/chembl_data_quality_report.md`

## IMPLEMENTATION TASKS

### Task 1: Complete Reference Compound Import
Ensure all reference cryoprotectants are properly imported with complete property data.

**Implementation Instructions:**

1. Update the reference compounds list in `chembl/reference_compounds.py`:

```python
def get_reference_compound_ids() -> List[str]:
    """
    Get list of reference cryoprotectant compound IDs.
    
    Returns:
        List of ChEMBL IDs for reference compounds
    """
    return [
        "CHEMBL388978",     # Glycerol
        "CHEMBL1098659",    # DMSO
        "CHEMBL66195",      # beta-Alanine
        "CHEMBL500033",     # tert-Butanol
        "CHEMBL1487",       # Urea
        "CHEMBL6196",       # Ethylene glycol
        "CHEMBL967",        # Propylene glycol
        "CHEMBL262548",     # Trehalose
        "CHEMBL6752"        # Glycine
    ]
```

2. Create a dedicated reference compounds import script that ensures complete property data:

```python
#!/usr/bin/env python3
"""
Reference compounds import script with complete property data import.

This script ensures that all reference cryoprotectant compounds are imported with
complete property data from both ChEMBL and PubChem sources.
"""

import os
import sys
import logging
import argparse
import json
from typing import Dict, List, Any, Optional
from datetime import datetime

# Import custom modules
from chembl.client import ChEMBLClient
from pubchem.client import PubChemClient
from cryoprotectant_identifiers import CryoprotectantIdentifierManager
from connection_pool_wrapper import ConnectionManager
from chembl.reference_compounds import get_reference_compound_ids

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

def import_reference_compounds(output_report: Optional[str] = None) -> Dict:
    """
    Import reference cryoprotectant compounds with complete property data.
    
    Args:
        output_report: Optional path to save the import report
        
    Returns:
        Import results dictionary
    """
    # Initialize clients
    chembl_client = ChEMBLClient()
    pubchem_client = PubChemClient()
    id_manager = CryoprotectantIdentifierManager.get_instance()
    
    # Get reference compound IDs
    chembl_ids = get_reference_compound_ids()
    logger.info(f"Importing {len(chembl_ids)} reference compounds")
    
    results = {
        'timestamp': datetime.now().isoformat(),
        'total_compounds': len(chembl_ids),
        'imported': 0,
        'updated': 0,
        'failed': 0,
        'details': {}
    }
    
    # Process each reference compound
    for chembl_id in chembl_ids:
        try:
            logger.info(f"Processing reference compound: {chembl_id}")
            
            # Step 1: Get compound data from ChEMBL
            chembl_data = chembl_client.get_compound(chembl_id)
            if not chembl_data:
                logger.error(f"Failed to fetch data for ChEMBL ID: {chembl_id}")
                results['failed'] += 1
                results['details'][chembl_id] = {'status': 'failed', 'error': 'ChEMBL data not found'}
                continue
                
            # Step 2: Get internal ID or create new one
            internal_id, is_new = id_manager.resolve_identifier(chembl_id=chembl_id)
            if not internal_id:
                internal_id = f"CRYO{len(id_manager.identifiers) + 1:04d}"
                is_new = True
            
            # Step 3: Get PubChem data if possible
            pubchem_cid = None
            pubchem_data = None
            if chembl_data.get('cross_references', {}).get('pubchem_cid'):
                pubchem_cid = chembl_data['cross_references']['pubchem_cid']
                pubchem_data = pubchem_client.get_compound(pubchem_cid)
            
            # Step 4: Merge properties from both sources
            properties = {
                'logs': {
                    'imported_at': datetime.now().isoformat(),
                    'source': 'reference_import'
                },
                'chembl': chembl_data,
                'basic': {
                    'name': chembl_data.get('pref_name') or chembl_data.get('molecule_properties', {}).get('full_molformula'),
                    'molecular_formula': chembl_data.get('molecule_properties', {}).get('full_molformula'),
                    'molecular_weight': chembl_data.get('molecule_properties', {}).get('full_mwt')
                },
                'identifiers': {
                    'chembl_id': chembl_id,
                    'pubchem_cid': pubchem_cid,
                    'inchi': chembl_data.get('molecule_structures', {}).get('standard_inchi'),
                    'inchi_key': chembl_data.get('molecule_structures', {}).get('standard_inchi_key'),
                    'smiles': chembl_data.get('molecule_structures', {}).get('canonical_smiles')
                },
                'properties': {
                    'alogp': chembl_data.get('molecule_properties', {}).get('alogp'),
                    'h_bond_donors': chembl_data.get('molecule_properties', {}).get('hbd'),
                    'h_bond_acceptors': chembl_data.get('molecule_properties', {}).get('hba'),
                    'rotatable_bonds': chembl_data.get('molecule_properties', {}).get('rtb'),
                    'psa': chembl_data.get('molecule_properties', {}).get('psa'),
                    'heavy_atoms': chembl_data.get('molecule_properties', {}).get('heavy_atoms')
                }
            }
            
            # Add PubChem properties if available
            if pubchem_data:
                properties['pubchem'] = pubchem_data
                
                # Extract and add PubChem-specific properties
                pubchem_props = {}
                if 'props' in pubchem_data:
                    for prop in pubchem_data.get('props', []):
                        if prop.get('urn', {}).get('label') == 'LogP':
                            pubchem_props['logP'] = prop.get('value', {}).get('sval')
                        elif prop.get('urn', {}).get('label') == 'Water Solubility':
                            pubchem_props['water_solubility'] = prop.get('value', {}).get('sval')
                        elif prop.get('urn', {}).get('label') == 'Melting Point':
                            pubchem_props['melting_point'] = prop.get('value', {}).get('sval')
                        elif prop.get('urn', {}).get('label') == 'Boiling Point':
                            pubchem_props['boiling_point'] = prop.get('value', {}).get('sval')
                
                # Merge with existing properties
                properties['properties'].update(pubchem_props)
            
            # Step 5: Update the identifier manager
            id_manager.add_molecule(internal_id, {
                "internal_id": internal_id,
                "chembl_id": chembl_id,
                "pubchem_cid": pubchem_cid,
                "names": [properties['basic']['name']],
                "inchi_key": properties['identifiers']['inchi_key'],
                "smiles": properties['identifiers']['smiles'],
                "formula": properties['basic']['molecular_formula'],
                "molecular_weight": properties['basic']['molecular_weight'],
                "category": "reference"
            })
            id_manager.save_identifiers()
            
            # Step 6: Insert or update in database
            with ConnectionManager() as conn:
                with conn.cursor() as cursor:
                    if is_new:
                        # Insert new molecule
                        insert_query = """
                        INSERT INTO molecules 
                        (id, name, properties, created_at, updated_at) 
                        VALUES (%s, %s, %s, NOW(), NOW())
                        """
                        cursor.execute(insert_query, 
                                      (internal_id, properties['basic']['name'], json.dumps(properties)))
                        
                        logger.info(f"Inserted new reference molecule {internal_id} (ChEMBL ID: {chembl_id})")
                        results['imported'] += 1
                    else:
                        # Update existing molecule
                        update_query = """
                        UPDATE molecules 
                        SET name = %s, properties = properties || %s, updated_at = NOW() 
                        WHERE id = %s
                        """
                        cursor.execute(update_query, 
                                      (properties['basic']['name'], json.dumps(properties), internal_id))
                        
                        logger.info(f"Updated reference molecule {internal_id} (ChEMBL ID: {chembl_id})")
                        results['updated'] += 1
                    
                    # Commit transaction
                    conn.commit()
            
            # Record success
            results['details'][chembl_id] = {
                'status': 'success',
                'internal_id': internal_id,
                'is_new': is_new,
                'has_pubchem': pubchem_cid is not None
            }
            
        except Exception as e:
            logger.error(f"Error processing reference compound {chembl_id}: {str(e)}")
            results['failed'] += 1
            results['details'][chembl_id] = {'status': 'failed', 'error': str(e)}
    
    # Generate report
    if output_report:
        os.makedirs(os.path.dirname(output_report), exist_ok=True)
        with open(output_report, 'w') as f:
            json.dump(results, f, indent=2)
        logger.info(f"Import report saved to {output_report}")
    
    logger.info(f"Reference compound import completed: {results['imported']} imported, " +
               f"{results['updated']} updated, {results['failed']} failed")
    
    return results

def main():
    """CLI entry point for reference compound import."""
    parser = argparse.ArgumentParser(description='Import reference cryoprotectant compounds')
    parser.add_argument('--report', 
                      default=f"reports/reference_import_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json",
                      help='Output file path for report')
    
    args = parser.parse_args()
    
    try:
        import_reference_compounds(args.report)
        return 0
    except Exception as e:
        logger.error(f"Import failed: {str(e)}")
        return 1

if __name__ == "__main__":
    sys.exit(main())
```

3. Run the reference compounds import script:

```bash
python import_reference_compounds.py
```

### Task 2: Enhance ChEMBL Data Import

Extend the full ChEMBL import to include more cryoprotectant compounds and ensure complete property data.

**Implementation Instructions:**

1. Create a search terms configuration file that includes more specific cryoprotectant categories:

```python
# chembl/search_terms.py
"""
Search terms for ChEMBL data import.

This module defines search terms for finding cryoprotectant-related compounds
in the ChEMBL database.
"""

CRYOPROTECTANT_CATEGORIES = [
    {
        "category": "polyols",
        "terms": ["glycerol", "sorbitol", "mannitol", "xylitol", "erythritol"]
    },
    {
        "category": "sugars",
        "terms": ["trehalose", "sucrose", "glucose", "lactose", "maltose", "dextran"]
    },
    {
        "category": "alcohols",
        "terms": ["ethanol", "methanol", "propanol", "butanol", "ethylene glycol", "propylene glycol"]
    },
    {
        "category": "amides",
        "terms": ["formamide", "acetamide", "urea", "hydroxyethyl starch"]
    },
    {
        "category": "sulfoxides",
        "terms": ["dimethyl sulfoxide", "DMSO"]
    },
    {
        "category": "amino_acids",
        "terms": ["glycine", "alanine", "proline", "glutamine", "lysine", "arginine"]
    },
]

def get_all_search_terms():
    """Get flattened list of all search terms."""
    terms = []
    for category in CRYOPROTECTANT_CATEGORIES:
        terms.extend(category["terms"])
    return terms

def get_terms_by_category(category):
    """Get search terms for a specific category."""
    for cat in CRYOPROTECTANT_CATEGORIES:
        if cat["category"] == category:
            return cat["terms"]
    return []
```

2. Modify the import_full_chembl.py script to use these search terms and ensure property import:

```python
# In the import_full_chembl.py script, add:

from chembl.search_terms import get_all_search_terms

def expand_search_results(chembl_client, search_results, limit=500):
    """
    Expand search results by adding similar compounds.
    
    Args:
        chembl_client: ChEMBL client instance
        search_results: Initial search results (ChEMBL IDs)
        limit: Maximum number of compounds to return
        
    Returns:
        Expanded list of ChEMBL IDs
    """
    expanded_results = set(search_results)
    
    # For each compound in the initial results
    for chembl_id in search_results:
        # Find similar compounds
        similar_compounds = chembl_client.get_similar_compounds(chembl_id, similarity=80)
        
        # Add to expanded results
        expanded_results.update(similar_compounds)
        
        # Stop if limit reached
        if len(expanded_results) >= limit:
            break
    
    return list(expanded_results)[:limit]
```

3. Add code to enhance the property extraction in ChEMBL import:

```python
def _extract_properties(self, compound_data: Dict[str, Any]) -> Dict[str, Any]:
    """
    Extract comprehensive property data from ChEMBL compound data.
    
    Args:
        compound_data: ChEMBL compound data
        
    Returns:
        Dictionary of extracted properties
    """
    properties = {}
    
    # Basic molecular properties
    molecule_props = compound_data.get('molecule_properties', {})
    if molecule_props:
        properties.update({
            'alogp': molecule_props.get('alogp'),
            'hba': molecule_props.get('hba'),
            'hbd': molecule_props.get('hbd'),
            'psa': molecule_props.get('psa'),
            'rtb': molecule_props.get('rtb'),
            'aromatic_rings': molecule_props.get('aromatic_rings'),
            'heavy_atoms': molecule_props.get('heavy_atoms'),
            'qed_weighted': molecule_props.get('qed_weighted'),
            'full_mwt': molecule_props.get('full_mwt'),
            'mw_freebase': molecule_props.get('mw_freebase'),
            'full_molformula': molecule_props.get('full_molformula'),
            'ro3_pass': molecule_props.get('ro3_pass') == 'Y'
        })
    
    # Physical properties (from bioactivities if available)
    bioactivities = compound_data.get('bioactivities', [])
    if bioactivities:
        # Extract properties from bioactivities
        for activity in bioactivities:
            if activity.get('type') == 'LogP':
                properties['experimental_logp'] = activity.get('value')
            elif activity.get('type') == 'Solubility':
                properties['solubility'] = activity.get('value')
    
    return properties
```

4. Run the enhanced ChEMBL import:

```bash
python import_full_chembl.py --search-terms-all --expand-similar --limit 500
```

### Task 3: Reconcile Cross-References

Create a script to establish cross-references between PubChem and ChEMBL identifiers for the same molecules.

**Implementation Instructions:**

1. Create `reconcile_chembl_properties.py`:

```python
#!/usr/bin/env python3
"""
Reconcile ChEMBL and PubChem properties and cross-references.

This script identifies molecules that exist in both ChEMBL and PubChem databases
and reconciles their properties and cross-references.
"""

import os
import sys
import logging
import argparse
import json
from typing import Dict, List, Any, Optional
from datetime import datetime

# Import custom modules
from connection_pool_wrapper import ConnectionManager
from cryoprotectant_identifiers import CryoprotectantIdentifierManager

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

def reconcile_cross_references(output_report: Optional[str] = None) -> Dict:
    """
    Reconcile cross-references between ChEMBL and PubChem.
    
    Args:
        output_report: Optional path to save the reconciliation report
        
    Returns:
        Reconciliation results dictionary
    """
    id_manager = CryoprotectantIdentifierManager.get_instance()
    
    results = {
        'timestamp': datetime.now().isoformat(),
        'molecules_updated': 0,
        'cross_references_added': 0,
        'conflicts_resolved': 0,
        'details': {}
    }
    
    # Step 1: Find molecules with PubChem or ChEMBL IDs
    with ConnectionManager() as conn:
        with conn.cursor() as cursor:
            # Get molecules with PubChem IDs
            cursor.execute("""
                SELECT id, name, properties->'identifiers'->>'pubchem_cid' as pubchem_cid
                FROM molecules
                WHERE properties->'identifiers'->>'pubchem_cid' IS NOT NULL
            """)
            pubchem_molecules = {row['id']: {'name': row['name'], 'pubchem_cid': row['pubchem_cid']} 
                               for row in cursor.fetchall()}
            
            # Get molecules with ChEMBL IDs
            cursor.execute("""
                SELECT id, name, properties->'identifiers'->>'chembl_id' as chembl_id
                FROM molecules
                WHERE properties->'identifiers'->>'chembl_id' IS NOT NULL
            """)
            chembl_molecules = {row['id']: {'name': row['name'], 'chembl_id': row['chembl_id']} 
                              for row in cursor.fetchall()}
    
    logger.info(f"Found {len(pubchem_molecules)} molecules with PubChem IDs")
    logger.info(f"Found {len(chembl_molecules)} molecules with ChEMBL IDs")
    
    # Step 2: Find molecules that should be linked based on shared identifiers
    inchi_key_map = {}
    with ConnectionManager() as conn:
        with conn.cursor() as cursor:
            # Get molecules with InChI Keys
            cursor.execute("""
                SELECT id, properties->'identifiers'->>'inchi_key' as inchi_key
                FROM molecules
                WHERE properties->'identifiers'->>'inchi_key' IS NOT NULL
            """)
            for row in cursor.fetchall():
                inchi_key = row['inchi_key']
                if inchi_key not in inchi_key_map:
                    inchi_key_map[inchi_key] = []
                inchi_key_map[inchi_key].append(row['id'])
    
    # Find molecules that share the same InChI Key
    reconciliation_candidates = []
    for inchi_key, molecule_ids in inchi_key_map.items():
        if len(molecule_ids) > 1:
            reconciliation_candidates.append({
                'inchi_key': inchi_key,
                'molecule_ids': molecule_ids
            })
    
    logger.info(f"Found {len(reconciliation_candidates)} InChI Keys with multiple molecules")
    
    # Step 3: Process reconciliation candidates
    for candidate in reconciliation_candidates:
        molecule_ids = candidate['molecule_ids']
        inchi_key = candidate['inchi_key']
        
        # Find molecules with PubChem and ChEMBL IDs in this group
        pubchem_id = None
        chembl_id = None
        pubchem_molecule = None
        chembl_molecule = None
        
        for mol_id in molecule_ids:
            if mol_id in pubchem_molecules:
                pubchem_id = mol_id
                pubchem_molecule = pubchem_molecules[mol_id]
            if mol_id in chembl_molecules:
                chembl_id = mol_id
                chembl_molecule = chembl_molecules[mol_id]
        
        # If we have one of each, we can reconcile
        if pubchem_id and chembl_id:
            logger.info(f"Reconciling molecules with InChI Key {inchi_key}: " +
                       f"{pubchem_id} (PubChem) and {chembl_id} (ChEMBL)")
            
            # Step 4: Update molecule properties and cross-references
            with ConnectionManager() as conn:
                with conn.cursor() as cursor:
                    # Update PubChem molecule with ChEMBL ID
                    cursor.execute("""
                        UPDATE molecules
                        SET properties = jsonb_set(
                            properties,
                            '{identifiers,chembl_id}',
                            %s::jsonb,
                            true
                        )
                        WHERE id = %s
                    """, (json.dumps(chembl_molecule['chembl_id']), pubchem_id))
                    
                    # Update ChEMBL molecule with PubChem CID
                    cursor.execute("""
                        UPDATE molecules
                        SET properties = jsonb_set(
                            properties,
                            '{identifiers,pubchem_cid}',
                            %s::jsonb,
                            true
                        )
                        WHERE id = %s
                    """, (json.dumps(pubchem_molecule['pubchem_cid']), chembl_id))
                    
                    # Commit transaction
                    conn.commit()
                    
                    results['cross_references_added'] += 2
                    results['molecules_updated'] += 2
                    results['details'][inchi_key] = {
                        'pubchem_molecule': pubchem_id,
                        'chembl_molecule': chembl_id,
                        'action': 'cross_references_added'
                    }
            
            # Step 5: Update identifier manager
            pubchem_cid = pubchem_molecule['pubchem_cid']
            chembl_id_str = chembl_molecule['chembl_id']
            
            # Update PubChem molecule's record
            pubchem_record = id_manager.get_molecule_by_internal_id(pubchem_id)
            if pubchem_record:
                pubchem_record['chembl_id'] = chembl_id_str
                id_manager.add_molecule(pubchem_id, pubchem_record)
            
            # Update ChEMBL molecule's record
            chembl_record = id_manager.get_molecule_by_internal_id(chembl_id)
            if chembl_record:
                chembl_record['pubchem_cid'] = pubchem_cid
                id_manager.add_molecule(chembl_id, chembl_record)
            
            id_manager.save_identifiers()
    
    # Generate report
    if output_report:
        os.makedirs(os.path.dirname(output_report), exist_ok=True)
        with open(output_report, 'w') as f:
            json.dump(results, f, indent=2)
        logger.info(f"Reconciliation report saved to {output_report}")
    
    logger.info(f"Reconciliation completed: {results['molecules_updated']} molecules updated, " +
               f"{results['cross_references_added']} cross-references added, " +
               f"{results['conflicts_resolved']} conflicts resolved")
    
    return results

def main():
    """CLI entry point for reconciliation."""
    parser = argparse.ArgumentParser(description='Reconcile ChEMBL and PubChem properties and cross-references')
    parser.add_argument('--report', 
                      default=f"reports/reconciliation_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json",
                      help='Output file path for report')
    
    args = parser.parse_args()
    
    try:
        reconcile_cross_references(args.report)
        return 0
    except Exception as e:
        logger.error(f"Reconciliation failed: {str(e)}")
        return 1

if __name__ == "__main__":
    sys.exit(main())
```

2. Run the reconciliation script:

```bash
python reconcile_chembl_properties.py
```

### Task 4: Enhance PubChem Import with Property Data

Extend the PubChem importer to ensure all molecules have complete property data.

**Implementation Instructions:**

1. Create a script to enhance PubChem molecules with missing properties:

```python
#!/usr/bin/env python3
"""
Enhance PubChem molecules with missing properties.

This script identifies PubChem molecules with missing properties and retrieves
the missing data from the PubChem API.
"""

import os
import sys
import logging
import argparse
import json
import time
from typing import Dict, List, Any, Optional
from datetime import datetime
import concurrent.futures

# Import custom modules
from connection_pool_wrapper import ConnectionManager
from pubchem.client import PubChemClient
from pubchem.rate_limiter import RateLimiter

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Constants
BATCH_SIZE = 10  # Number of compounds to process in parallel
MAX_RETRIES = 3  # Maximum number of retries for failed operations

def enhance_pubchem_properties(output_report: Optional[str] = None, batch_size: int = BATCH_SIZE) -> Dict:
    """
    Enhance PubChem molecules with missing properties.
    
    Args:
        output_report: Optional path to save the enhancement report
        batch_size: Number of compounds to process in parallel
        
    Returns:
        Enhancement results dictionary
    """
    pubchem_client = PubChemClient()
    rate_limiter = RateLimiter(requests_per_minute=5)
    
    results = {
        'timestamp': datetime.now().isoformat(),
        'total_molecules': 0,
        'molecules_enhanced': 0,
        'molecules_skipped': 0,
        'molecules_failed': 0,
        'details': {}
    }
    
    # Step 1: Find PubChem molecules with missing properties
    with ConnectionManager() as conn:
        with conn.cursor() as cursor:
            cursor.execute("""
                SELECT id, name, properties->'identifiers'->>'pubchem_cid' as pubchem_cid
                FROM molecules
                WHERE properties->'identifiers'->>'pubchem_cid' IS NOT NULL
                AND (
                    properties->'properties'->>'logP' IS NULL OR
                    properties->'properties'->>'h_bond_donors' IS NULL OR
                    properties->'properties'->>'h_bond_acceptors' IS NULL
                )
            """)
            molecules_to_enhance = cursor.fetchall()
    
    results['total_molecules'] = len(molecules_to_enhance)
    logger.info(f"Found {results['total_molecules']} PubChem molecules with missing properties")
    
    # Step 2: Process molecules in batches
    for i in range(0, len(molecules_to_enhance), batch_size):
        batch = molecules_to_enhance[i:i+batch_size]
        logger.info(f"Processing batch {i//batch_size + 1}/{(len(molecules_to_enhance) + batch_size - 1)//batch_size}")
        
        # Process batch in parallel
        with concurrent.futures.ThreadPoolExecutor(max_workers=batch_size) as executor:
            future_to_molecule = {
                executor.submit(_enhance_molecule, molecule, pubchem_client, rate_limiter): molecule
                for molecule in batch
            }
            
            # Process results as they complete
            for future in concurrent.futures.as_completed(future_to_molecule):
                molecule = future_to_molecule[future]
                try:
                    success, is_enhanced = future.result()
                    
                    if success:
                        if is_enhanced:
                            results['molecules_enhanced'] += 1
                            results['details'][molecule['id']] = {
                                'status': 'enhanced',
                                'name': molecule['name'],
                                'pubchem_cid': molecule['pubchem_cid']
                            }
                        else:
                            results['molecules_skipped'] += 1
                            results['details'][molecule['id']] = {
                                'status': 'skipped',
                                'name': molecule['name'],
                                'pubchem_cid': molecule['pubchem_cid']
                            }
                    else:
                        results['molecules_failed'] += 1
                        results['details'][molecule['id']] = {
                            'status': 'failed',
                            'name': molecule['name'],
                            'pubchem_cid': molecule['pubchem_cid']
                        }
                        
                except Exception as e:
                    logger.error(f"Error processing molecule {molecule['id']}: {str(e)}")
                    results['molecules_failed'] += 1
                    results['details'][molecule['id']] = {
                        'status': 'failed',
                        'name': molecule['name'],
                        'pubchem_cid': molecule['pubchem_cid'],
                        'error': str(e)
                    }
        
        # Avoid overtaxing the PubChem API
        time.sleep(1)
    
    # Generate report
    if output_report:
        os.makedirs(os.path.dirname(output_report), exist_ok=True)
        with open(output_report, 'w') as f:
            json.dump(results, f, indent=2)
        logger.info(f"Enhancement report saved to {output_report}")
    
    logger.info(f"Enhancement completed: {results['molecules_enhanced']} molecules enhanced, " +
               f"{results['molecules_skipped']} molecules skipped, " +
               f"{results['molecules_failed']} molecules failed")
    
    return results

def _enhance_molecule(molecule: Dict, pubchem_client: PubChemClient, rate_limiter: RateLimiter) -> (bool, bool):
    """
    Enhance a single molecule with missing properties.
    
    Args:
        molecule: Molecule data (id, name, pubchem_cid)
        pubchem_client: PubChem client instance
        rate_limiter: Rate limiter instance
        
    Returns:
        Tuple of (success, is_enhanced)
    """
    molecule_id = molecule['id']
    pubchem_cid = molecule['pubchem_cid']
    
    logger.debug(f"Enhancing molecule {molecule_id} (CID {pubchem_cid})")
    
    try:
        # Wait for rate limiter
        rate_limiter.wait()
        
        # Get compound data from PubChem
        compound_data = pubchem_client.get_compound(pubchem_cid)
        if not compound_data:
            logger.warning(f"Failed to fetch data for PubChem CID: {pubchem_cid}")
            return False, False
        
        # Extract properties
        properties = {}
        if 'props' in compound_data:
            for prop in compound_data.get('props', []):
                if prop.get('urn', {}).get('label') == 'LogP':
                    properties['logP'] = prop.get('value', {}).get('sval')
                elif prop.get('urn', {}).get('label') == 'Water Solubility':
                    properties['water_solubility'] = prop.get('value', {}).get('sval')
                elif prop.get('urn', {}).get('label') == 'Melting Point':
                    properties['melting_point'] = prop.get('value', {}).get('sval')
                elif prop.get('urn', {}).get('label') == 'Boiling Point':
                    properties['boiling_point'] = prop.get('value', {}).get('sval')
                elif prop.get('urn', {}).get('label') == 'Complexity':
                    properties['complexity'] = prop.get('value', {}).get('ival')
                elif prop.get('urn', {}).get('label') == 'H-Bond Donor':
                    properties['h_bond_donors'] = prop.get('value', {}).get('ival')
                elif prop.get('urn', {}).get('label') == 'H-Bond Acceptor':
                    properties['h_bond_acceptors'] = prop.get('value', {}).get('ival')
                elif prop.get('urn', {}).get('label') == 'Rotatable Bond':
                    properties['rotatable_bonds'] = prop.get('value', {}).get('ival')
                elif prop.get('urn', {}).get('label') == 'Heavy Atom':
                    properties['heavy_atoms'] = prop.get('value', {}).get('ival')
        
        # If no properties extracted, skip
        if not properties:
            logger.warning(f"No properties found for PubChem CID: {pubchem_cid}")
            return True, False
        
        # Update molecule in database
        with ConnectionManager() as conn:
            with conn.cursor() as cursor:
                cursor.execute("""
                    UPDATE molecules
                    SET properties = jsonb_set(
                        properties,
                        '{properties}',
                        properties->'properties' || %s::jsonb,
                        true
                    )
                    WHERE id = %s
                """, (json.dumps(properties), molecule_id))
                
                # Commit transaction
                conn.commit()
        
        logger.debug(f"Enhanced molecule {molecule_id} with properties: {properties.keys()}")
        return True, True
    
    except Exception as e:
        logger.error(f"Error enhancing molecule {molecule_id} (CID {pubchem_cid}): {str(e)}")
        return False, False

def main():
    """CLI entry point for property enhancement."""
    parser = argparse.ArgumentParser(description='Enhance PubChem molecules with missing properties')
    parser.add_argument('--batch-size', type=int, default=BATCH_SIZE,
                      help='Number of molecules to process in parallel')
    parser.add_argument('--report', 
                      default=f"reports/pubchem_enhancement_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json",
                      help='Output file path for report')
    
    args = parser.parse_args()
    
    try:
        enhance_pubchem_properties(args.report, args.batch_size)
        return 0
    except Exception as e:
        logger.error(f"Enhancement failed: {str(e)}")
        return 1

if __name__ == "__main__":
    sys.exit(main())
```

2. Run the enhancement script:

```bash
python enhance_pubchem_properties.py
```

### Task 5: Performance Optimization

Optimize database queries to ensure acceptable performance.

**Implementation Instructions:**

1. Create a script to add performance indexes:

```python
#!/usr/bin/env python3
"""
Add performance indexes to database tables.

This script adds indexes to improve query performance for common operations.
"""

import os
import sys
import logging
import argparse
import time
from typing import Dict, List, Any, Optional
from connection_pool_wrapper import ConnectionManager

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

def add_performance_indexes() -> Dict:
    """
    Add performance indexes to database tables.
    
    Returns:
        Dictionary with results
    """
    results = {
        'indexes_added': 0,
        'tables_modified': set(),
        'details': []
    }
    
    # Define indexes to add
    indexes = [
        {
            'name': 'molecules_name_idx',
            'table': 'molecules',
            'columns': ['name'],
            'description': 'Index on molecule name for faster name-based searches'
        },
        {
            'name': 'molecules_properties_pubchem_cid_idx',
            'table': 'molecules',
            'expression': "((properties->'identifiers'->>'pubchem_cid'))",
            'description': 'Expression index on PubChem CID for faster identifier-based searches'
        },
        {
            'name': 'molecules_properties_chembl_id_idx',
            'table': 'molecules',
            'expression': "((properties->'identifiers'->>'chembl_id'))",
            'description': 'Expression index on ChEMBL ID for faster identifier-based searches'
        },
        {
            'name': 'molecules_properties_inchi_key_idx',
            'table': 'molecules',
            'expression': "((properties->'identifiers'->>'inchi_key'))",
            'description': 'Expression index on InChI Key for faster structure-based searches'
        },
        {
            'name': 'mixture_components_mixture_id_idx',
            'table': 'mixture_components',
            'columns': ['mixture_id'],
            'description': 'Index on mixture_id for faster mixture component lookups'
        },
        {
            'name': 'mixture_components_molecule_id_idx',
            'table': 'mixture_components',
            'columns': ['molecule_id'],
            'description': 'Index on molecule_id for faster compound usage lookups'
        }
    ]
    
    # Add each index
    with ConnectionManager() as conn:
        with conn.cursor() as cursor:
            for index in indexes:
                try:
                    # Check if index already exists
                    cursor.execute("""
                        SELECT 1 FROM pg_indexes
                        WHERE indexname = %s
                    """, (index['name'],))
                    
                    if cursor.fetchone():
                        logger.info(f"Index {index['name']} already exists, skipping")
                        continue
                    
                    # Create the index
                    if 'expression' in index:
                        # Expression index
                        sql = f"""
                            CREATE INDEX {index['name']}
                            ON {index['table']} ({index['expression']})
                        """
                    else:
                        # Regular index
                        sql = f"""
                            CREATE INDEX {index['name']}
                            ON {index['table']} ({', '.join(index['columns'])})
                        """
                    
                    start_time = time.time()
                    cursor.execute(sql)
                    
                    # Add comment with description
                    comment_sql = f"""
                        COMMENT ON INDEX {index['name']}
                        IS %s
                    """
                    cursor.execute(comment_sql, (index['description'],))
                    
                    # Commit transaction
                    conn.commit()
                    
                    elapsed_time = time.time() - start_time
                    logger.info(f"Added index {index['name']} in {elapsed_time:.2f} seconds")
                    
                    results['indexes_added'] += 1
                    results['tables_modified'].add(index['table'])
                    results['details'].append({
                        'name': index['name'],
                        'table': index['table'],
                        'elapsed_time': elapsed_time,
                        'description': index['description']
                    })
                    
                except Exception as e:
                    logger.error(f"Error adding index {index['name']}: {str(e)}")
                    conn.rollback()
    
    results['tables_modified'] = list(results['tables_modified'])
    logger.info(f"Added {results['indexes_added']} indexes to {len(results['tables_modified'])} tables")
    
    return results

def main():
    """CLI entry point for adding performance indexes."""
    parser = argparse.ArgumentParser(description='Add performance indexes to database tables')
    
    try:
        results = add_performance_indexes()
        logger.info(f"Performance optimization completed. Added {results['indexes_added']} indexes.")
        return 0
    except Exception as e:
        logger.error(f"Performance optimization failed: {str(e)}")
        return 1

if __name__ == "__main__":
    sys.exit(main())
```

2. Run the performance optimization script:

```bash
python add_performance_indexes.py
```

## COMPREHENSIVE DATABASE POPULATION EXECUTION PLAN

To successfully populate the database with complete scientific data, follow these steps in sequence:

1. **Connection Pool Setup**
   - Verify connection pool wrapper implementation
   - Test connection pool functionality

2. **Reference Compounds Import**
   - Create and run import_reference_compounds.py script
   - Verify all reference compounds are imported with complete properties

3. **ChEMBL Data Import**
   - Create search_terms.py with expanded cryoprotectant categories
   - Enhance import_full_chembl.py to use search terms and ensure property import
   - Run ChEMBL import with 500 compound limit

4. **PubChem Property Enhancement**
   - Run enhance_pubchem_properties.py to add missing properties to existing molecules

5. **Cross-Reference Reconciliation**
   - Run reconcile_chembl_properties.py to establish cross-references

6. **Performance Optimization**
   - Run add_performance_indexes.py to optimize query performance

7. **Verification**
   - Run verify_imported_data.py to verify database population
   - Generate final data quality report

## VERIFICATION STEPS

After completing the database population, verify its success with these steps:

1. **Count imported molecules:**
   ```sql
   SELECT 
     COUNT(*) AS total_molecules,
     COUNT(properties->'identifiers'->>'pubchem_cid') AS pubchem_molecules,
     COUNT(properties->'identifiers'->>'chembl_id') AS chembl_molecules,
     COUNT(CASE WHEN properties->'identifiers'->>'pubchem_cid' IS NOT NULL AND 
                     properties->'identifiers'->>'chembl_id' IS NOT NULL 
           THEN 1 END) AS cross_referenced_molecules
   FROM molecules;
   ```

2. **Verify reference compounds:**
   ```sql
   SELECT id, name, properties->'identifiers'->>'chembl_id' AS chembl_id
   FROM molecules
   WHERE properties->'identifiers'->>'chembl_id' IN (
     'CHEMBL388978', 'CHEMBL1098659', 'CHEMBL66195', 'CHEMBL500033', 
     'CHEMBL1487', 'CHEMBL6196', 'CHEMBL967', 'CHEMBL262548', 'CHEMBL6752'
   );
   ```

3. **Verify property completeness:**
   ```sql
   SELECT 
     COUNT(*) AS total,
     COUNT(properties->'properties'->>'logP') AS logp_count,
     COUNT(properties->'properties'->>'h_bond_donors') AS hbond_donors_count,
     COUNT(properties->'properties'->>'h_bond_acceptors') AS hbond_acceptors_count
   FROM molecules
   WHERE properties->'identifiers'->>'pubchem_cid' IS NOT NULL;
   ```

4. **Verify query performance:**
   ```sql
   EXPLAIN ANALYZE
   SELECT * FROM molecules
   WHERE properties->'identifiers'->>'pubchem_cid' = '962';
   ```

5. **Run full verification script:**
   ```bash
   python verify_imported_data.py --full-verification
   ```

## TROUBLESHOOTING COMMON ISSUES

1. **Missing property data:**
   - Verify PubChem API responses contain the expected property data
   - Check property extraction functions to ensure they handle all property formats
   - Run the enhance_pubchem_properties.py script with verbose logging

2. **Missing cross-references:**
   - Verify InChI Key generation for consistent molecule identification
   - Check both PubChem and ChEMBL sources provide cross-references
   - Run reconcile_chembl_properties.py with manual matching logic

3. **Performance issues:**
   - Verify indexes have been properly created and are being used (check query plans)
   - Check for high-volume tables that may need partitioning
   - Consider connection pool tuning (min/max connections)

4. **Import failures:**
   - Check checkpoint files for resumable imports
   - Verify rate limiting to avoid API throttling
   - Increase timeout parameters for large data imports