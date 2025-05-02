# Minimal Population Workflow Guide

## Simplified Workflow Principles

1. **Focus on Success Metrics**: Target one verification criterion at a time
2. **Readable Over Optimized**: Prioritize code simplicity over performance 
3. **Observable Progress**: Add clear checkpoints and logging
4. **Atomic Operations**: Complete each molecule's full processing in one transaction
5. **Hard-Coded Reference Data**: Use guaranteed working data for reference compounds

## Execution Steps

1. **First verify database connection and schema**
   - Check that essential tables exist (molecules, molecular_properties, property_types)
   - Verify connection with a simple query
   - Create any missing tables or columns with simple DDL

2. **Populate and verify reference compounds (9/9)**
   - Use hard-coded data for all 9 reference compounds
   - Insert base molecule data and properties in one transaction per molecule
   - Verify that all 9 reference compounds have the 3 critical properties
   - Only proceed once 9/9 complete

3. **Import and verify small batch of PubChem data (100-500 molecules)**
   - Start with a small batch size (10-20 compounds)
   - Focus on complete properties over volume
   - Implement file-based checkpointing after each small batch
   - Verify critical properties after each batch
   - Gradually increase batch size once success is proven

4. **Import and verify small batch of ChEMBL data (100-500 molecules)**
   - Use reference compounds as seeds for finding similar compounds
   - Implement same methodical, cautious approach
   - Prioritize property completeness over volume
   - Create detailed progress logs

5. **Once basic functionality is verified:**
   - Expand to full dataset
   - Add performance optimizations
   - Run full verification suite

## Implementation Guidelines

### Database Connection

Use a simple, consistent connection function:

```python
def get_connection():
    """Get a simple database connection without complexity."""
    import os
    import psycopg2
    from psycopg2.extras import RealDictCursor
    
    db_host = os.getenv('SUPABASE_DB_HOST')
    db_port = os.getenv('SUPABASE_DB_PORT', '6543')  # Session pooler port
    db_name = os.getenv('SUPABASE_DB_NAME')
    db_user = os.getenv('SUPABASE_DB_USER')
    db_password = os.getenv('SUPABASE_DB_PASSWORD')
    
    conn = psycopg2.connect(
        host=db_host,
        port=db_port,
        dbname=db_name,
        user=db_user,
        password=db_password,
        cursor_factory=RealDictCursor
    )
    
    return conn
```

### Simple Transaction Handling

Use straightforward try/except blocks:

```python
def with_transaction(conn, func, *args, **kwargs):
    """Execute a function within a transaction with simple error handling."""
    try:
        # Start transaction
        conn.autocommit = False
        
        # Execute function
        result = func(*args, **kwargs)
        
        # Commit transaction
        conn.commit()
        return result
        
    except Exception as e:
        # Rollback on error
        print(f"Transaction error: {str(e)}")
        conn.rollback()
        raise
        
    finally:
        # Ensure connection returns to autocommit mode
        conn.autocommit = True
```

### Property Handling

Focus on the 3 critical properties that affect verification:

```python
def set_critical_properties(conn, molecule_id, properties):
    """Set the 3 critical properties for verification."""
    # Get property type IDs
    cursor = conn.cursor()
    
    # Create property types if they don't exist
    property_types = {}
    for prop_name in ['logP', 'h_bond_donors', 'h_bond_acceptors']:
        cursor.execute("""
            INSERT INTO property_types (name, data_type)
            VALUES (%s, 'numeric')
            ON CONFLICT (name) DO UPDATE SET updated_at = NOW()
            RETURNING id
        """, (prop_name,))
        property_types[prop_name] = cursor.fetchone()['id']
    
    # Set property values
    for prop_name, prop_id in property_types.items():
        if prop_name in properties:
            cursor.execute("""
                INSERT INTO molecular_properties (molecule_id, property_type_id, numeric_value)
                VALUES (%s, %s, %s)
                ON CONFLICT (molecule_id, property_type_id) 
                DO UPDATE SET numeric_value = EXCLUDED.numeric_value
            """, (molecule_id, prop_id, properties[prop_name]))
    
    return True
```

### Verification Functions

Keep verification simple and focused:

```python
def verify_molecule_properties(molecule_id):
    """Verify a single molecule has all required properties."""
    conn = get_connection()
    cursor = conn.cursor()
    
    # Get molecule details
    cursor.execute("""
        SELECT name, pubchem_cid, chembl_id FROM molecules WHERE id = %s
    """, (molecule_id,))
    molecule = cursor.fetchone()
    
    if not molecule:
        print(f"Molecule with ID {molecule_id} not found")
        return False
    
    # Check critical properties
    required_props = ['logP', 'h_bond_donors', 'h_bond_acceptors']
    for prop in required_props:
        cursor.execute("""
            SELECT mp.numeric_value
            FROM molecular_properties mp
            JOIN property_types pt ON mp.property_type_id = pt.id
            WHERE mp.molecule_id = %s AND pt.name = %s
        """, (molecule_id, prop))
        
        result = cursor.fetchone()
        if not result or result['numeric_value'] is None:
            print(f"Molecule {molecule['name']} missing property: {prop}")
            return False
    
    print(f"Molecule {molecule['name']} has all required properties")
    return True
```

## Debugging Process

1. **Start with single, hard-coded molecule insert**
   - Insert a single reference compound with complete properties
   - Verify the insert worked correctly
   - Review database state

2. **Add each stage with careful verification**
   - Only proceed to next stage when current one succeeds
   - Verify after each small batch
   - Create checkpoints frequently

3. **Isolate database operations from API calls**
   - Separate data fetching from database operations
   - Test database operations with mock data
   - Handle API-specific errors separately

4. **Log every database operation result**
   - Log success/failure of each operation
   - Include error details when operations fail
   - Use consistent logging format

5. **Verify each property is stored correctly**
   - Check property values after insertion
   - Compare inserted values with source data
   - Handle type conversions carefully

## Recommended Folder Structure

- `simplified/`: Contains all simplified implementation scripts
  - `reference_compounds.py`: Reference compounds population
  - `pubchem_minimal.py`: Minimal PubChem implementation
  - `chembl_minimal.py`: Minimal ChEMBL implementation
  - `verification.py`: Simple verification utilities
  - `db_helpers.py`: Simplified database utilities
  - `execute_simplified.py`: Main execution script

## Execution Script Example

```python
#!/usr/bin/env python3
"""
Simplified Database Population Script

Focuses on reliable population of CryoProtect database
using a simplified, step-by-step approach.
"""

import argparse
import os
import sys
import time

# Import simplified modules
from simplified.db_helpers import get_connection, check_connection, verify_schema
from simplified.reference_compounds import populate_reference_compounds, verify_reference_compounds
from simplified.pubchem_minimal import import_pubchem_minimal
from simplified.chembl_minimal import import_chembl_minimal
from simplified.verification import verify_all_molecules, verify_performance

def main():
    """Main execution function."""
    parser = argparse.ArgumentParser(description="Populate CryoProtect database with simplified approach")
    parser.add_argument("--step", type=str, choices=["check", "references", "pubchem", "chembl", "all", "verify"],
                      help="Step to execute", required=True)
    parser.add_argument("--limit", type=int, default=500,
                      help="Maximum number of compounds to import (for pubchem/chembl steps)")
    args = parser.parse_args()
    
    # Check connection first
    if args.step == "check" or args.step == "all":
        print("Checking database connection and schema...")
        if not check_connection():
            print("Database connection failed. Please check your credentials.")
            return 1
        if not verify_schema():
            print("Database schema verification failed. Please check your database setup.")
            return 1
        print("Database connection and schema verification successful!")
    
    # Populate reference compounds
    if args.step == "references" or args.step == "all":
        print("Populating reference compounds...")
        success_count = populate_reference_compounds()
        
        if success_count == 9:
            print("Successfully populated all 9 reference compounds!")
        else:
            print(f"Warning: Only populated {success_count}/9 reference compounds.")
        
        # Verify reference compounds
        if verify_reference_compounds():
            print("Reference compound verification passed!")
        else:
            print("Reference compound verification failed!")
            if args.step == "all":
                print("Stopping 'all' process due to reference verification failure.")
                return 1
    
    # Import PubChem data
    if args.step == "pubchem" or args.step == "all":
        print(f"Importing up to {args.limit} PubChem compounds...")
        pubchem_count = import_pubchem_minimal(args.limit)
        print(f"Imported {pubchem_count} PubChem compounds.")
    
    # Import ChEMBL data
    if args.step == "chembl" or args.step == "all":
        print(f"Importing up to {args.limit} ChEMBL compounds...")
        chembl_count = import_chembl_minimal(args.limit)
        print(f"Imported {chembl_count} ChEMBL compounds.")
    
    # Verify all
    if args.step == "verify" or args.step == "all":
        print("Verifying database population...")
        if verify_all_molecules():
            print("Database verification PASSED! All criteria met.")
        else:
            print("Database verification FAILED. See above for details.")
        
        # Check performance
        print("Checking query performance...")
        avg_time = verify_performance()
        print(f"Average query time: {avg_time:.2f}ms")
        if avg_time < 50:
            print("Performance requirement met!")
        else:
            print("Performance requirement not met. Consider adding indexes.")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())
```

## Common Pitfalls to Avoid

1. **Over-optimization too early**
   - Don't worry about performance until basic functionality works
   - Simple queries are easier to debug than complex ones

2. **Batch size too large**
   - Start with very small batches (10-20 compounds)
   - Only increase batch size after proven success

3. **Insufficient error handling**
   - Log all errors in detail
   - Handle different error types differently
   - Don't proceed when critical errors occur

4. **Incomplete property data**
   - Verify property values after insertion
   - Skip molecules with missing critical properties
   - Log property failures separately

5. **Connection instability**
   - Use simple, manual reconnection logic
   - Don't over-engineer connection pooling early
   - Log all connection issues

## Milestone Checklist

- [ ] Connection and schema verification successful
- [ ] 9/9 reference compounds populated with all properties
- [ ] 100+ PubChem compounds with complete properties
- [ ] 100+ ChEMBL compounds with complete properties
- [ ] 500+ total molecules with 90%+ property completeness
- [ ] 5,000+ total molecules with 90%+ property completeness
- [ ] Query performance under 50ms

## References

- `database/connection.py`: Original connection utilities
- `verify_imported_data.py`: Verification script
- `project_state.json`: Project state tracking
- `import_pubchem_data_direct.py`: Original PubChem import script
- `import_full_chembl.py`: Original ChEMBL import script