# Database Population Remediation Plan

## Overview

This document outlines the plan to redesign the data import logic for the CryoProtect database population directive. The previous implementation failed validation due to a schema mismatch: the import scripts were written assuming a 'properties' JSONB column in the 'molecules' table, but the actual database uses a normalized schema with separate 'property_types' and 'molecular_properties' tables.

## Current Schema Structure

### Molecules Table
```sql
CREATE TABLE public.molecules (
    id UUID PRIMARY KEY DEFAULT uuid_generate_v4(),
    cid INTEGER UNIQUE NOT NULL,  -- PubChem Compound ID
    name TEXT,
    molecular_formula TEXT,
    smiles TEXT,
    inchi TEXT,
    inchikey CHARACTER VARYING,
    formula CHARACTER VARYING,
    molecular_weight NUMERIC,
    created_at TIMESTAMPTZ NOT NULL DEFAULT NOW(),
    updated_at TIMESTAMPTZ NOT NULL DEFAULT NOW(),
    created_by UUID REFERENCES auth.users(id),
    pubchem_cid INTEGER,
    pubchem_link TEXT GENERATED ALWAYS AS ('https://pubchem.ncbi.nlm.nih.gov/compound/' || pubchem_cid) STORED,
    chembl_id VARCHAR(20)
);
```

### Property Types Table
```sql
CREATE TABLE public.property_types (
    id UUID PRIMARY KEY DEFAULT uuid_generate_v4(),
    name TEXT UNIQUE NOT NULL,
    description TEXT,
    units TEXT,
    data_type TEXT NOT NULL,  -- 'numeric', 'text', 'boolean', etc.
    created_at TIMESTAMPTZ NOT NULL DEFAULT NOW(),
    updated_at TIMESTAMPTZ NOT NULL DEFAULT NOW()
);
```

### Molecular Properties Table
```sql
CREATE TABLE public.molecular_properties (
    id UUID PRIMARY KEY DEFAULT uuid_generate_v4(),
    molecule_id UUID NOT NULL REFERENCES public.molecules(id) ON DELETE CASCADE,
    property_type_id UUID NOT NULL REFERENCES public.property_types(id),
    numeric_value NUMERIC,
    text_value TEXT,
    boolean_value BOOLEAN,
    created_at TIMESTAMPTZ NOT NULL DEFAULT NOW(),
    updated_at TIMESTAMPTZ NOT NULL DEFAULT NOW(),
    created_by UUID REFERENCES auth.users(id),
    UNIQUE(molecule_id, property_type_id)
);
```

## Remediation Strategy

The remediation strategy involves:

1. Creating a common utility module (`property_utils.py`) to handle property operations with the normalized schema
2. Modifying each import script to use this utility module instead of directly manipulating a JSONB column
3. Creating a verification script to validate the database population

## Implementation Plan

### 1. Create Common Utility Module

First, we'll create a `property_utils.py` module that provides a consistent interface for working with molecular properties in the normalized schema. This will encapsulate the complexity of property type management and value insertion.

### 2. Modify Import Scripts

We'll update each of the four import scripts to use the PropertyManager class from our utility module:

- `import_reference_compounds.py`
- `import_full_chembl.py`
- `reconcile_chembl_properties.py`
- `enhance_pubchem_properties.py`

### 3. Create Verification Script
## Detailed Implementation

### Common Utility Module: `property_utils.py`

```python
#!/usr/bin/env python3
"""
Utility functions for working with molecular properties in the normalized database schema.
"""

import logging
from typing import Dict, Any, Optional, Union, Tuple
from uuid import UUID
import json

logger = logging.getLogger(__name__)

class PropertyManager:
    """
    Manages property operations for the normalized database schema.
    """
    
    def __init__(self, connection_manager):
        """
        Initialize the PropertyManager.
        
        Args:
            connection_manager: Database connection manager
        """
        self.connection_manager = connection_manager
        self._property_types_cache = {}
        self._load_property_types()
    
    def _load_property_types(self):
        """Load property types from the database into cache."""
        try:
            with self.connection_manager() as conn:
                with conn.cursor() as cursor:
                    cursor.execute("""
                        SELECT id, name, data_type FROM property_types
                    """)
                    for row in cursor.fetchall():
                        self._property_types_cache[row['name']] = {
                            'id': row['id'],
                            'data_type': row['data_type']
                        }
            logger.debug(f"Loaded {len(self._property_types_cache)} property types from database")
        except Exception as e:
            logger.error(f"Error loading property types: {str(e)}")
            # Continue with empty cache, will create property types as needed
    
    def get_property_type_id(self, property_name: str, data_type: str = 'numeric') -> UUID:
        """
        Get the property type ID for a given property name, creating it if it doesn't exist.
        
        Args:
            property_name: Name of the property
            data_type: Data type of the property ('numeric', 'text', 'boolean')
            
        Returns:
            UUID of the property type
        """
        # Check cache first
        if property_name in self._property_types_cache:
            return self._property_types_cache[property_name]['id']
        
        # Not in cache, try to find in database or create
        try:
            with self.connection_manager() as conn:
                with conn.cursor() as cursor:
                    # Try to find existing property type
                    cursor.execute("""
                        SELECT id, data_type FROM property_types WHERE name = %s
                    """, (property_name,))
                    row = cursor.fetchone()
                    
                    if row:
                        # Found existing property type
                        property_type_id = row['id']
                        self._property_types_cache[property_name] = {
                            'id': property_type_id,
                            'data_type': row['data_type']
                        }
                        return property_type_id
                    
                    # Property type doesn't exist, create it
                    cursor.execute("""
                        INSERT INTO property_types (name, data_type, description)
                        VALUES (%s, %s, %s)
                        RETURNING id
                    """, (property_name, data_type, f"Auto-created property: {property_name}"))
                    
                    property_type_id = cursor.fetchone()['id']
                    conn.commit()
                    
                    # Update cache
                    self._property_types_cache[property_name] = {
                        'id': property_type_id,
                        'data_type': data_type
                    }
                    
                    logger.info(f"Created new property type: {property_name} ({data_type})")
                    return property_type_id
                    
        except Exception as e:
            logger.error(f"Error getting/creating property type for {property_name}: {str(e)}")
            raise
    
    def set_property(self, molecule_id: UUID, property_name: str, 
                    property_value: Any, created_by: Optional[UUID] = None) -> bool:
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
        # Determine data type based on value
        if property_value is None:
            return True  # Skip None values
            
        if isinstance(property_value, (int, float)):
            data_type = 'numeric'
        elif isinstance(property_value, bool):
            data_type = 'boolean'
        else:
            data_type = 'text'
            # Convert to string if not already
            if not isinstance(property_value, str):
                property_value = str(property_value)
        
        try:
            # Get property type ID
            property_type_id = self.get_property_type_id(property_name, data_type)
            
            with self.connection_manager() as conn:
                with conn.cursor() as cursor:
                    # Check if property already exists
                    cursor.execute("""
                        SELECT id FROM molecular_properties
                        WHERE molecule_id = %s AND property_type_id = %s
                    """, (molecule_id, property_type_id))
                    
                    existing_property = cursor.fetchone()
                    
                    if existing_property:
                        # Update existing property
                        if data_type == 'numeric':
                            cursor.execute("""
                                UPDATE molecular_properties
                                SET numeric_value = %s, text_value = NULL, boolean_value = NULL, updated_at = NOW()
                                WHERE molecule_id = %s AND property_type_id = %s
                            """, (property_value, molecule_id, property_type_id))
                        elif data_type == 'boolean':
                            cursor.execute("""
                                UPDATE molecular_properties
                                SET numeric_value = NULL, text_value = NULL, boolean_value = %s, updated_at = NOW()
                                WHERE molecule_id = %s AND property_type_id = %s
                            """, (property_value, molecule_id, property_type_id))
                        else:  # text
                            cursor.execute("""
                                UPDATE molecular_properties
                                SET numeric_value = NULL, text_value = %s, boolean_value = NULL, updated_at = NOW()
                                WHERE molecule_id = %s AND property_type_id = %s
                            """, (property_value, molecule_id, property_type_id))
                    else:
                        # Insert new property
                        if data_type == 'numeric':
                            cursor.execute("""
                                INSERT INTO molecular_properties
                                (molecule_id, property_type_id, numeric_value, created_by)
                                VALUES (%s, %s, %s, %s)
                            """, (molecule_id, property_type_id, property_value, created_by))
                        elif data_type == 'boolean':
                            cursor.execute("""
                                INSERT INTO molecular_properties
                                (molecule_id, property_type_id, boolean_value, created_by)
                                VALUES (%s, %s, %s, %s)
                            """, (molecule_id, property_type_id, property_value, created_by))
                        else:  # text
                            cursor.execute("""
                                INSERT INTO molecular_properties
                                (molecule_id, property_type_id, text_value, created_by)
                                VALUES (%s, %s, %s, %s)
                            """, (molecule_id, property_type_id, property_value, created_by))
                    
                    conn.commit()
                    return True
                    
        except Exception as e:
            logger.error(f"Error setting property {property_name}={property_value} for molecule {molecule_id}: {str(e)}")
            return False
    
    def set_properties(self, molecule_id: UUID, properties: Dict[str, Any], 
                      created_by: Optional[UUID] = None) -> Tuple[int, int]:
        """
        Set multiple properties for a molecule.
        
        Args:
            molecule_id: UUID of the molecule
            properties: Dictionary of property name -> value
            created_by: UUID of the user creating the properties
            
        Returns:
            Tuple of (success_count, total_count)
        """
        success_count = 0
        total_count = len(properties)
        
        for property_name, property_value in properties.items():
            if self.set_property(molecule_id, property_name, property_value, created_by):
                success_count += 1
        
        return success_count, total_count
    
    def get_properties(self, molecule_id: UUID, property_names: Optional[list] = None) -> Dict[str, Any]:
        """
        Get multiple properties for a molecule.
        
        Args:
            molecule_id: UUID of the molecule
            property_names: Optional list of property names to retrieve (all if None)
            
        Returns:
            Dictionary of property name -> value
        """
        result = {}
        
        try:
            with self.connection_manager() as conn:
                with conn.cursor() as cursor:
                    if property_names:
                        # Get specific properties
                        placeholders = ', '.join(['%s'] * len(property_names))
                        cursor.execute(f"""
                            SELECT pt.name, mp.numeric_value, mp.text_value, mp.boolean_value, pt.data_type
                            FROM molecular_properties mp
                            JOIN property_types pt ON mp.property_type_id = pt.id
                            WHERE mp.molecule_id = %s AND pt.name IN ({placeholders})
                        """, [molecule_id] + property_names)
                    else:
                        # Get all properties
                        cursor.execute("""
                            SELECT pt.name, mp.numeric_value, mp.text_value, mp.boolean_value, pt.data_type
                            FROM molecular_properties mp
                            JOIN property_types pt ON mp.property_type_id = pt.id
                            WHERE mp.molecule_id = %s
                        """, (molecule_id,))
                    
                    for row in cursor.fetchall():
                        property_name = row['name']
                        data_type = row['data_type']
                        
                        # Get the appropriate value based on data type
                        if data_type == 'numeric':
                            result[property_name] = row['numeric_value']
                        elif data_type == 'boolean':
                            result[property_name] = row['boolean_value']
                        else:  # text
                            result[property_name] = row['text_value']
                    
        except Exception as e:
            logger.error(f"Error getting properties for molecule {molecule_id}: {str(e)}")
        
        return result
```

We'll create a `verify_imported_data.py` script to validate that the database has been properly populated.

### 4. Testing Plan

We'll test each modified script individually and then run an end-to-end test to ensure all components work together correctly.
### Script Modifications

#### 1. Modify `import_reference_compounds.py`

The key changes required for this script are:

1. Update imports to include the PropertyManager
2. Modify the database insertion logic to store basic molecule data in the molecules table
3. Use PropertyManager to insert properties into the molecular_properties table

```python
# Add to imports
from property_utils import PropertyManager

# Modify the import_reference_compounds function
def import_reference_compounds(output_report: Optional[str] = None) -> Dict:
    # ... existing code ...
    
    # Initialize PropertyManager
    property_manager = PropertyManager(ConnectionManager)
    
    # ... existing code ...
    
    # When inserting/updating molecules, modify the database operation:
    with ConnectionManager() as conn:
        with conn.cursor() as cursor:
            if is_new:
                # Insert new molecule with basic fields only (no properties JSONB)
                insert_query = """
                INSERT INTO molecules 
                (id, name, smiles, inchi, inchikey, formula, molecular_weight, chembl_id, pubchem_cid, created_at, updated_at) 
                VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, NOW(), NOW())
                """
                cursor.execute(insert_query, 
                              (internal_id, 
                               properties['basic']['name'],
                               properties['identifiers']['smiles'],
                               properties['identifiers']['inchi'],
                               properties['identifiers']['inchi_key'],
                               properties['basic']['molecular_formula'],
                               properties['basic']['molecular_weight'],
                               chembl_id,
                               pubchem_cid))
                
                logger.info(f"Inserted new reference molecule {internal_id} (ChEMBL ID: {chembl_id})")
                results['imported'] += 1
            else:
                # Update existing molecule basic fields
                update_query = """
                UPDATE molecules 
                SET name = %s, smiles = %s, inchi = %s, inchikey = %s, 
                    formula = %s, molecular_weight = %s, 
                    chembl_id = %s, pubchem_cid = %s, updated_at = NOW() 
                WHERE id = %s
                """
                cursor.execute(update_query, 
                              (properties['basic']['name'],
                               properties['identifiers']['smiles'],
                               properties['identifiers']['inchi'],
                               properties['identifiers']['inchi_key'],
                               properties['basic']['molecular_formula'],
                               properties['basic']['molecular_weight'],
                               chembl_id,
                               pubchem_cid,
                               internal_id))
                
                logger.info(f"Updated reference molecule {internal_id} (ChEMBL ID: {chembl_id})")
                results['updated'] += 1
            
            # Commit transaction for basic molecule data
            conn.commit()
            
            # Now insert properties using PropertyManager
            property_data = {}
            
            # Add properties from ChEMBL
            if 'properties' in properties:
                for prop_name, prop_value in properties['properties'].items():
                    if prop_value is not None:
                        property_data[prop_name] = prop_value
            
            # Insert properties
            success_count, total_count = property_manager.set_properties(
                internal_id, property_data, created_by=None)
            
            logger.info(f"Added {success_count}/{total_count} properties for molecule {internal_id}")
```

#### 2. Modify `import_full_chembl.py`

The key changes required for this script are:

1. Update the ChEMBLWorker class to use the PropertyManager
2. Modify the database insertion logic to handle the normalized schema

```python
# Add to imports in the worker module
from property_utils import PropertyManager

# Modify the process_compound method in ChEMBLWorker
def process_compound(self, compound_id, dry_run=False):
    # ... existing code ...
    
    # Initialize PropertyManager
    property_manager = PropertyManager(ConnectionManager)
    
    # ... existing code ...
    
    # When inserting molecule data, modify the database operation:
    with ConnectionManager() as conn:
        with conn.cursor() as cursor:
            # Insert or update basic molecule data
            cursor.execute("""
                INSERT INTO molecules 
                (id, name, smiles, inchi, inchikey, formula, molecular_weight, chembl_id, created_at, updated_at) 
                VALUES (%s, %s, %s, %s, %s, %s, %s, %s, NOW(), NOW())
                ON CONFLICT (id) DO UPDATE 
                SET name = EXCLUDED.name, 
                    smiles = EXCLUDED.smiles,
                    inchi = EXCLUDED.inchi,
                    inchikey = EXCLUDED.inchikey,
                    formula = EXCLUDED.formula,
                    molecular_weight = EXCLUDED.molecular_weight,
                    chembl_id = EXCLUDED.chembl_id,
                    updated_at = NOW()
            """, (molecule_id, 
                  molecule_data['name'],
                  molecule_data.get('smiles'),
                  molecule_data.get('inchi'),
                  molecule_data.get('inchi_key'),
                  molecule_data.get('molecular_formula'),
                  molecule_data.get('molecular_weight'),
                  compound_id))
            
            # Commit transaction for basic molecule data
            conn.commit()
    
    # Insert properties using PropertyManager
    if 'properties' in molecule_data:
        property_data = molecule_data['properties']
        success_count, total_count = property_manager.set_properties(
            molecule_id, property_data, created_by=None)
        
        logger.info(f"Added {success_count}/{total_count} properties for molecule {molecule_id}")
```

#### 3. Modify `reconcile_chembl_properties.py`

The key changes required for this script are:

1. Update the queries to work with the normalized schema
2. No need to use PropertyManager as this script only updates cross-references

```python
# Modify the reconcile_cross_references function
def reconcile_cross_references(output_report: Optional[str] = None, dry_run: bool = False) -> Dict:
    # ... existing code ...
    
    # Step 1: Find molecules with PubChem or ChEMBL IDs
    try:
        with ConnectionManager() as conn:
            with conn.cursor() as cursor:
                # Get molecules with PubChem IDs
                cursor.execute("""
                    SELECT id, name, pubchem_cid
                    FROM molecules
                    WHERE pubchem_cid IS NOT NULL
                """)
                pubchem_molecules = {row['id']: {'name': row['name'], 'pubchem_cid': row['pubchem_cid']}
                                   for row in cursor.fetchall()}
                
                # Get molecules with ChEMBL IDs
                cursor.execute("""
                    SELECT id, name, chembl_id
                    FROM molecules
                    WHERE chembl_id IS NOT NULL
                """)
                chembl_molecules = {row['id']: {'name': row['name'], 'chembl_id': row['chembl_id']}
                                  for row in cursor.fetchall()}
    # ... existing code ...
    
    # Step 2: Find molecules that should be linked based on shared identifiers
    if not dry_run:
        inchi_key_map = {}
        try:
            with ConnectionManager() as conn:
                with conn.cursor() as cursor:
                    # Get molecules with InChI Keys
                    cursor.execute("""
                        SELECT id, inchikey
                        FROM molecules
                        WHERE inchikey IS NOT NULL
                    """)
                    for row in cursor.fetchall():
                        inchi_key = row['inchikey']
                        if inchi_key not in inchi_key_map:
                            inchi_key_map[inchi_key] = []
                        inchi_key_map[inchi_key].append(row['id'])
        # ... existing code ...
    
    # Step 4: Update molecule properties and cross-references
    if not dry_run:
        try:
            with ConnectionManager() as conn:
                with conn.cursor() as cursor:
                    # Update PubChem molecule with ChEMBL ID
                    cursor.execute("""
                        UPDATE molecules
                        SET chembl_id = %s, updated_at = NOW()
                        WHERE id = %s
                    """, (chembl_molecule['chembl_id'], pubchem_id))
                    
                    # Update ChEMBL molecule with PubChem CID
                    cursor.execute("""
                        UPDATE molecules
                        SET pubchem_cid = %s, updated_at = NOW()
                        WHERE id = %s
                    """, (pubchem_molecule['pubchem_cid'], chembl_id))
                    
                    # Commit transaction
                    conn.commit()
        # ... existing code ...
```

#### 4. Modify `enhance_pubchem_properties.py`

The key changes required for this script are:

1. Update the script to use the PropertyManager for property insertion
2. Modify the database queries to work with the normalized schema

```python
# Add to imports
from property_utils import PropertyManager

# Modify the enhance_pubchem_properties function
def enhance_pubchem_properties(output_report: Optional[str] = None, batch_size: int = BATCH_SIZE, dry_run: bool = DRY_RUN) -> Dict:
    # ... existing code ...
    
    # Step 1: Find PubChem molecules with missing properties
    if not dry_run:
        # Query database for molecules with missing properties
        with ConnectionManager() as conn:
            with conn.cursor() as cursor:
                # Initialize PropertyManager
                property_manager = PropertyManager(ConnectionManager)
                
                # Get molecules with PubChem CIDs
                cursor.execute("""
                    SELECT m.id, m.name, m.pubchem_cid
                    FROM molecules m
                    WHERE m.pubchem_cid IS NOT NULL
                """)
                all_molecules = cursor.fetchall()
                
                # Filter molecules with missing properties
                molecules_to_enhance = []
                for molecule in all_molecules:
                    # Get existing properties
                    properties = property_manager.get_properties(molecule['id'])
                    
                    # Check if key properties are missing
                    if ('logP' not in properties or 
                        'h_bond_donors' not in properties or 
                        'h_bond_acceptors' not in properties):
                        molecules_to_enhance.append(molecule)
    # ... existing code ...
    
    # Modify the _enhance_molecule function
    def _enhance_molecule(molecule: Tuple, pubchem_client: PubChemClient, rate_limiter: RateLimiter, dry_run: bool = DRY_RUN) -> Tuple[bool, bool]:
        # ... existing code ...
        
        # Update molecule in database
        if not dry_run:
            # Initialize PropertyManager
            property_manager = PropertyManager(ConnectionManager)
            
            # Insert properties using PropertyManager
            success_count, total_count = property_manager.set_properties(
                molecule_id, properties, created_by=None)
            
            if success_count > 0:
                logger.debug(f"Enhanced molecule {molecule_id} with {success_count}/{total_count} properties")
                return True, True
            else:
                logger.warning(f"Failed to add any properties for molecule {molecule_id}")
                return False, False
        else:
            # In dry run mode, just log the properties that would be updated
            logger.info(f"[DRY RUN] Would update molecule {molecule_id} with properties: {list(properties.keys())}")
            return True, True
```
### Verification Script

To verify the database population after implementing the changes, we'll create a verification script:

```python
#!/usr/bin/env python3
"""
Verification script for database population.

This script verifies that the database has been properly populated with
reference compounds, ChEMBL data, and property data.
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
from property_utils import PropertyManager
from chembl.reference_compounds import get_reference_compound_ids

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

def verify_database_population(output_report: Optional[str] = None) -> Dict:
    """
    Verify database population.
    
    Args:
        output_report: Optional path to save the verification report
        
    Returns:
        Verification results dictionary
    """
    results = {
        'timestamp': datetime.now().isoformat(),
        'tests': {},
        'overall_status': 'PASS'
    }
    
    # Initialize PropertyManager
    property_manager = PropertyManager(ConnectionManager)
    
    # Test 1: Count imported molecules
    try:
        with ConnectionManager() as conn:
            with conn.cursor() as cursor:
                cursor.execute("""
                    SELECT 
                      COUNT(*) AS total_molecules,
                      COUNT(pubchem_cid) AS pubchem_molecules,
                      COUNT(chembl_id) AS chembl_molecules,
                      COUNT(CASE WHEN pubchem_cid IS NOT NULL AND chembl_id IS NOT NULL 
                            THEN 1 END) AS cross_referenced_molecules
                    FROM molecules
                """)
                
                counts = cursor.fetchone()
                
                results['tests']['molecule_counts'] = {
                    'total_molecules': counts['total_molecules'],
                    'pubchem_molecules': counts['pubchem_molecules'],
                    'chembl_molecules': counts['chembl_molecules'],
                    'cross_referenced_molecules': counts['cross_referenced_molecules'],
                    'status': 'PASS' if counts['total_molecules'] > 0 else 'FAIL'
                }
                
                if counts['chembl_molecules'] < 500:
                    results['tests']['molecule_counts']['status'] = 'FAIL'
                    results['tests']['molecule_counts']['message'] = f"Expected >=500 ChEMBL molecules, found {counts['chembl_molecules']}"
                    results['overall_status'] = 'FAIL'
    except Exception as e:
        results['tests']['molecule_counts'] = {
            'status': 'ERROR',
            'message': str(e)
        }
        results['overall_status'] = 'FAIL'
    
    # Test 2: Verify reference compounds
    try:
        reference_ids = get_reference_compound_ids()
        
        with ConnectionManager() as conn:
            with conn.cursor() as cursor:
                placeholders = ', '.join(['%s'] * len(reference_ids))
                cursor.execute(f"""
                    SELECT id, name, chembl_id
                    FROM molecules
                    WHERE chembl_id IN ({placeholders})
                """, reference_ids)
                
                reference_compounds = cursor.fetchall()
                
                found_ids = [row['chembl_id'] for row in reference_compounds]
                missing_ids = [ref_id for ref_id in reference_ids if ref_id not in found_ids]
                
                results['tests']['reference_compounds'] = {
                    'found': len(reference_compounds),
                    'expected': len(reference_ids),
                    'missing': missing_ids,
                    'status': 'PASS' if len(reference_compounds) == len(reference_ids) else 'FAIL'
                }
                
                if len(reference_compounds) < len(reference_ids):
                    results['tests']['reference_compounds']['message'] = f"Missing reference compounds: {missing_ids}"
                    results['overall_status'] = 'FAIL'
    except Exception as e:
        results['tests']['reference_compounds'] = {
            'status': 'ERROR',
            'message': str(e)
        }
        results['overall_status'] = 'FAIL'
    
    # Test 3: Verify property completeness
    try:
        with ConnectionManager() as conn:
            with conn.cursor() as cursor:
                # Get molecules with PubChem CIDs
                cursor.execute("""
                    SELECT id, pubchem_cid
                    FROM molecules
                    WHERE pubchem_cid IS NOT NULL
                    LIMIT 100  -- Check a sample for performance
                """)
                
                pubchem_molecules = cursor.fetchall()
                
                # Check properties for each molecule
                property_stats = {
                    'total': len(pubchem_molecules),
                    'with_logp': 0,
                    'with_h_bond_donors': 0,
                    'with_h_bond_acceptors': 0,
                    'complete': 0
                }
                
                for molecule in pubchem_molecules:
                    properties = property_manager.get_properties(molecule['id'])
                    
                    has_logp = 'logP' in properties
                    has_hbd = 'h_bond_donors' in properties
                    has_hba = 'h_bond_acceptors' in properties
                    
                    if has_logp:
                        property_stats['with_logp'] += 1
                    if has_hbd:
                        property_stats['with_h_bond_donors'] += 1
                    if has_hba:
                        property_stats['with_h_bond_acceptors'] += 1
                    
                    if has_logp and has_hbd and has_hba:
                        property_stats['complete'] += 1
                
                # Calculate percentages
                if property_stats['total'] > 0:
                    property_stats['complete_percent'] = (property_stats['complete'] / property_stats['total']) * 100
                    property_stats['logp_percent'] = (property_stats['with_logp'] / property_stats['total']) * 100
                    property_stats['hbd_percent'] = (property_stats['with_h_bond_donors'] / property_stats['total']) * 100
                    property_stats['hba_percent'] = (property_stats['with_h_bond_acceptors'] / property_stats['total']) * 100
                
                results['tests']['property_completeness'] = property_stats
                
                # Pass if at least 80% of molecules have complete properties
                if property_stats['total'] > 0 and property_stats['complete_percent'] >= 80:
                    results['tests']['property_completeness']['status'] = 'PASS'
                else:
                    results['tests']['property_completeness']['status'] = 'FAIL'
                    results['tests']['property_completeness']['message'] = f"Only {property_stats['complete_percent']:.1f}% of molecules have complete properties"
                    results['overall_status'] = 'FAIL'
    except Exception as e:
        results['tests']['property_completeness'] = {
            'status': 'ERROR',
            'message': str(e)
        }
        results['overall_status'] = 'FAIL'
    
    # Test 4: Verify cross-references
    try:
        with ConnectionManager() as conn:
            with conn.cursor() as cursor:
                cursor.execute("""
                    SELECT COUNT(*) AS cross_ref_count
                    FROM molecules
                    WHERE pubchem_cid IS NOT NULL AND chembl_id IS NOT NULL
                """)
                
                cross_ref_count = cursor.fetchone()['cross_ref_count']
                
                results['tests']['cross_references'] = {
                    'count': cross_ref_count,
                    'status': 'PASS' if cross_ref_count > 0 else 'FAIL'
                }
                
                if cross_ref_count == 0:
                    results['tests']['cross_references']['message'] = "No molecules have both PubChem CID and ChEMBL ID"
                    results['overall_status'] = 'FAIL'
    except Exception as e:
        results['tests']['cross_references'] = {
            'status': 'ERROR',
            'message': str(e)
        }
        results['overall_status'] = 'FAIL'
    
    # Test 5: Verify query performance
    try:
        with ConnectionManager() as conn:
            with conn.cursor() as cursor:
                # Time a typical query
                start_time = datetime.now()
                
                cursor.execute("""
                    SELECT m.id, m.name, m.pubchem_cid, m.chembl_id,
                           pt.name as property_name, 
                           mp.numeric_value, mp.text_value, mp.boolean_value
                    FROM molecules m
                    JOIN molecular_properties mp ON m.id = mp.molecule_id
                    JOIN property_types pt ON mp.property_type_id = pt.id
                    WHERE m.pubchem_cid = 962  -- Glycerol
                """)
                
                query_results = cursor.fetchall()
                
                end_time = datetime.now()
                query_time_ms = (end_time - start_time).total_seconds() * 1000
                
                results['tests']['query_performance'] = {
                    'query_time_ms': query_time_ms,
                    'result_count': len(query_results),
                    'status': 'PASS' if query_time_ms < 100 else 'FAIL'
                }
                
                if query_time_ms >= 100:
                    results['tests']['query_performance']['message'] = f"Query took {query_time_ms:.1f}ms, expected <100ms"
                    results['overall_status'] = 'FAIL'
    except Exception as e:
        results['tests']['query_performance'] = {
            'status': 'ERROR',
            'message': str(e)
        }
        results['overall_status'] = 'FAIL'
    
    # Generate report
    if output_report:
        os.makedirs(os.path.dirname(output_report), exist_ok=True)
        with open(output_report, 'w') as f:
            json.dump(results, f, indent=2)
        logger.info(f"Verification report saved to {output_report}")
    
    # Log summary
    logger.info(f"Verification completed with status: {results['overall_status']}")
    for test_name, test_result in results['tests'].items():
        logger.info(f"  {test_name}: {test_result['status']}")
    
    return results

def main():
    """CLI entry point for verification."""
    parser = argparse.ArgumentParser(description='Verify database population')
    parser.add_argument('--report', 
                      default=f"reports/final_data_quality_report.md",
                      help='Output file path for report')
    parser.add_argument('--full-verification', action='store_true',
                      help='Run full verification (may be slower)')
    
    args = parser.parse_args()
    
    try:
        results = verify_database_population(args.report)
        
        # Exit with appropriate status code
        if results['overall_status'] == 'PASS':
            return 0
        else:
            return 1
    except Exception as e:
        logger.error(f"Verification failed: {str(e)}")
        return 1

if __name__ == "__main__":
    sys.exit(main())
```

## Implementation Plan

The implementation will proceed in the following phases:

### Phase 1: Create Common Utility Module

1. Create the `property_utils.py` module with the PropertyManager class
2. Test the PropertyManager with a small sample of data to ensure it works correctly

### Phase 2: Modify Import Scripts

1. Modify `import_reference_compounds.py` and test with a small set of reference compounds
2. Modify `import_full_chembl.py` and test with a limited number of compounds
3. Modify `reconcile_chembl_properties.py` and test with a small set of cross-references
4. Modify `enhance_pubchem_properties.py` and test with a small set of molecules

### Phase 3: Create Verification Script

1. Create the `verify_imported_data.py` script
2. Test the verification script against the current database state

### Phase 4: Full Database Population

1. Run the modified scripts in sequence:
   - `import_reference_compounds.py`
   - `import_full_chembl.py`
   - `enhance_pubchem_properties.py`
   - `reconcile_chembl_properties.py`
2. Run the verification script to ensure all data has been properly imported

## Conclusion

This remediation plan addresses the schema mismatch issue by modifying the import scripts to work with the normalized database schema. The key components of the solution are:

1. A common PropertyManager class that handles the complexity of working with the normalized schema
2. Modified import scripts that use the PropertyManager to insert properties
3. A verification script to ensure the database is properly populated

By implementing this plan, we will ensure that all reference compounds, ChEMBL data, and property data are correctly imported into the database, meeting the success criteria specified in the directive.