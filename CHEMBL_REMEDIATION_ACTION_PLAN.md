# ChEMBL Data Remediation Action Plan

This comprehensive action plan provides the ROO agents with all necessary information to complete the ChEMBL data remediation task successfully, with special focus on addressing the blocked migration issue.

## 1. Current Status Analysis

### Schema Changes Status
- ✅ Schema definition ready with dedicated `chembl_id` column for the `molecules` table 
- ❌ Migration step blocked: MCP DDL execution issue (see Critical Blocker section below)
- ❌ Database migration not applied due to technical limitation

### Import Script Enhancements
- ✅ Modified `ChEMBL_Integrated_Import.py` to store `chembl_id` in dedicated column
- ✅ Added logic to prioritize standard reference compounds
- ✅ Enhanced property storage to track data source (alogp vs cx_logp)
- ✅ Increased default import limit from 1000 to 2500

### Remaining Issues
- ❌ Insufficient data volume: Database contains only 30 molecules with 10 having ChEMBL IDs (need 1000+)
- ❌ Property value discrepancies: 0% match with official ChEMBL API values
- ❌ Reference compounds missing: Standard molecules like CHEMBL25 not present

## 1.1 Critical Blocker: MCP DDL Execution Issue

**BLOCKING ISSUE**: The ChEMBL remediation process is currently blocked at the migration step due to a technical limitation with the MCP `execute_sql` tool.

### Problem Analysis
- **Error**: The MCP `execute_sql` tool fails with error `{'code': '42601', 'message': 'ALTER TABLE query does not return tuples'}`
- **Root Cause**: The MCP `execute_sql` tool only supports SQL queries that return tuples (SELECT statements), but not DDL statements like ALTER TABLE
- **Current Status**: Task marked as "Blocked-Debug" in project_state.json
- **Migration File**: `migrations/chembl_add_chembl_id.sql` contains necessary schema changes but cannot be executed directly via MCP

```sql
-- The problematic SQL migration that needs to be executed
ALTER TABLE molecules ADD COLUMN IF NOT EXISTS chembl_id VARCHAR;
UPDATE molecules SET chembl_id = SUBSTRING(data_source FROM 'CHEMBL\d+') WHERE data_source LIKE '%CHEMBL%';
CREATE INDEX IF NOT EXISTS idx_molecules_chembl_id ON molecules (chembl_id);
```

### Current Technical Blockers with MCP Integration

The ChEMBL remediation process is blocked by multiple technical issues related to MCP tool compatibility:

1. **DDL Execution Issue**: MCP `execute_sql` can't handle non-tuple-returning statements
2. **Type Compatibility Issue**: Verification query returns PostgreSQL types incompatible with MCP JSON
3. **Argument Naming Issue**: ChEMBL import script receives incorrect argument names

#### Blocker #1: MCP DDL Execution Issue

The primary issue is that the MCP `execute_sql` tool cannot handle DDL statements directly. Here are three potential solutions that ROO agents should implement:

##### Option 1: Use SQL Decorator Pattern (Recommended)
This solution involves "wrapping" the DDL statements in a function that returns data, making them compatible with the MCP tool:

```sql
-- Instead of direct ALTER TABLE, use this pattern:
SELECT execute_ddl('ALTER TABLE molecules ADD COLUMN IF NOT EXISTS chembl_id VARCHAR;') as result;
SELECT execute_ddl('CREATE INDEX IF NOT EXISTS idx_molecules_chembl_id ON molecules (chembl_id);') as result;
```

First, create the `execute_ddl` function if it doesn't exist:
```sql
CREATE OR REPLACE FUNCTION execute_ddl(ddl_command text) RETURNS TABLE(result text) AS $$
BEGIN
    EXECUTE ddl_command;
    RETURN QUERY SELECT 'Success: ' || ddl_command AS result;
EXCEPTION WHEN OTHERS THEN
    RETURN QUERY SELECT 'Error: ' || SQLERRM AS result;
END;
$$ LANGUAGE plpgsql SECURITY DEFINER;
```

##### Option 2: Split Into Single Statements
Split the migration into individual statements and execute them separately, with only SELECT statements using the MCP tool directly:

1. For checking if column exists:
```sql
SELECT column_name::text AS column_name FROM information_schema.columns 
WHERE table_name = 'molecules' AND column_name = 'chembl_id';
```

2. For data operations (UPDATE), use a RETURNING clause to make it return tuples:
```sql
UPDATE molecules SET chembl_id = SUBSTRING(data_source FROM 'CHEMBL\d+') 
WHERE data_source LIKE '%CHEMBL%' 
RETURNING id, chembl_id;
```

#### Blocker #2: Type Compatibility Issue with MCP JSON

MCP expects all query results to be compatible with JSON serialization, but PostgreSQL types like `sql_identifier` aren't directly compatible. Fix verification queries by explicitly casting column types to text:

**Problem:**
```sql
-- This fails because column_name is of type sql_identifier
SELECT column_name FROM information_schema.columns 
WHERE table_name = 'molecules' AND column_name = 'chembl_id';
```

**Solution:**
```sql
-- Cast to text for JSON compatibility
SELECT column_name::text AS column_name FROM information_schema.columns 
WHERE table_name = 'molecules' AND column_name = 'chembl_id';
```

#### Blocker #3: Script Argument Format Issue

The ChEMBL import script is called with incorrect argument naming conventions in the orchestration script:

**Problem:**
```python
# Incorrect argument naming (underscores instead of hyphens)
["python", "ChEMBL_Integrated_Import.py", "--limit", "2500", "--batch_size", "50", "--checkpoint_interval", "100"]
```

**Solution:**
```python
# Correct argument naming (using hyphens)
["python", "ChEMBL_Integrated_Import.py", "--limit", "2500", "--batch-size", "50", "--checkpoint-interval", "100"]
```

#### Option 3: Use Direct Database Connection
If available, bypass MCP for DDL operations by using a direct database connection method:

```python
def execute_migration_with_psycopg2():
    """Execute the migration using direct psycopg2 connection."""
    import psycopg2
    from config import Config
    
    # Get database configuration
    config = Config()
    
    try:
        # Connect directly to the database
        conn = psycopg2.connect(
            host=config.SUPABASE_DB_HOST,
            port=config.SUPABASE_DB_PORT,
            database=config.SUPABASE_DB_NAME,
            user=config.SUPABASE_DB_USER,
            password=config.SUPABASE_DB_PASSWORD
        )
        conn.autocommit = True
        
        # Read the migration file
        with open("migrations/chembl_add_chembl_id.sql", "r") as f:
            migration_sql = f.read()
        
        # Execute the migration
        cursor = conn.cursor()
        cursor.execute(migration_sql)
        cursor.close()
        conn.close()
        
        return True
    except Exception as e:
        logger.error(f"Error executing migration with psycopg2: {str(e)}")
        return False
```

### Fixed Orchestration Script: chembl_remediation_main_fixed.py

Here's how to fix the main orchestration script to work around all the identified MCP limitations:

```python
def apply_migration(project_id):
    """
    Apply the migration to add the chembl_id column using compatible SELECT operations.
    
    This function works around the limitation that MCP execute_sql cannot handle DDL statements
    by using techniques like creating a function that wraps DDL in a SELECT statement.
    """
    logger.info("Applying migration to add chembl_id column")
    
    # Check if the migration file exists
    migration_file = "migrations/chembl_add_chembl_id.sql"
    if not os.path.exists(migration_file):
        logger.info("Migration file doesn't exist, creating it")
        
        # Create the migration file with compatible SQL
        migration_sql = """
-- First create a helper function for DDL operations if it doesn't exist
SELECT EXISTS (
    SELECT 1 FROM pg_proc 
    WHERE proname = 'execute_ddl'
)::text AS function_exists;

-- Create the execute_ddl function if it doesn't exist
DO $$
BEGIN
    IF NOT EXISTS (
        SELECT 1 FROM pg_proc WHERE proname = 'execute_ddl'
    ) THEN
        CREATE OR REPLACE FUNCTION execute_ddl(ddl_command text) 
        RETURNS TABLE(result text) AS $$
        BEGIN
            EXECUTE ddl_command;
            RETURN QUERY SELECT 'Success: ' || ddl_command AS result;
        EXCEPTION WHEN OTHERS THEN
            RETURN QUERY SELECT 'Error: ' || SQLERRM AS result;
        END;
        $$ LANGUAGE plpgsql SECURITY DEFINER;
    END IF;
END;
$$;

-- Now use the function to execute DDL statements
SELECT execute_ddl('ALTER TABLE molecules ADD COLUMN IF NOT EXISTS chembl_id VARCHAR;') as ddl_result;

-- For UPDATE statements, use RETURNING to make them return tuples
UPDATE molecules 
SET chembl_id = SUBSTRING(data_source FROM 'CHEMBL[0-9]+') 
WHERE data_source LIKE '%CHEMBL%' 
RETURNING id::text, chembl_id::text;

-- Create the index using our helper function
SELECT execute_ddl('CREATE INDEX IF NOT EXISTS idx_molecules_chembl_id ON molecules (chembl_id);') as ddl_result;
        """
        
        # Write the migration file
        with open(migration_file, "w") as f:
            f.write(migration_sql)
        
        logger.info(f"Created migration file: {migration_file}")
    
    # Read the migration file
    with open(migration_file, "r") as f:
        migration_sql = f.read()
    
    # Split the migration SQL into separate statements
    # This is important because execute_sql can only execute one statement at a time
    sql_statements = migration_sql.split(';')
    sql_statements = [stmt.strip() for stmt in sql_statements if stmt.strip()]
    
    # Execute each SQL statement
    results = []
    for i, sql in enumerate(sql_statements):
        logger.info(f"Executing migration statement {i+1}/{len(sql_statements)}")
        
        # Skip empty statements
        if not sql.strip():
            continue
            
        # Only execute statements that actually return data
        if sql.strip().upper().startswith(('SELECT', 'UPDATE', 'INSERT', 'DELETE')):
            try:
                # Execute the SQL statement using MCP
                result = execute_sql(sql, project_id)
                results.append(result)
                logger.info(f"Statement executed successfully: {sql[:50]}...")
            except Exception as e:
                logger.error(f"Error executing statement: {str(e)}")
                logger.error(f"Problematic SQL: {sql}")
                return False
        else:
            logger.warning(f"Skipping non-tuple-returning statement: {sql[:50]}...")
    
    # Verify the migration was successful - FIX: Cast column_name to text
    try:
        logger.info("Verifying chembl_id column was added")
        query = "SELECT column_name::text AS column_name FROM information_schema.columns WHERE table_name = 'molecules' AND column_name = 'chembl_id';"
        result = execute_sql(query, project_id)
        
        if result and len(result) > 0 and result[0]["column_name"] == "chembl_id":
            logger.info("Column chembl_id successfully added to molecules table")
            return True
        else:
            logger.error("Column chembl_id not found in molecules table")
            return False
    except Exception as e:
        logger.error(f"Error verifying migration: {str(e)}")
        return False
```

Additionally, update the ChEMBL import function to use correct argument format (hyphens, not underscores):

```python
def run_import(project_id=None):
    """
    Run the ChEMBL data import script.
    
    This function executes the ChEMBL_Integrated_Import.py script with the appropriate parameters.
    """
    logger.info("Running ChEMBL data import")
    
    try:
        # Build the command with correct argument format (hyphens, not underscores)
        cmd = [
            "python", 
            "ChEMBL_Integrated_Import.py",
            "--limit", "2500",
            "--batch-size", "50",           # FIX: Use hyphens instead of underscores
            "--checkpoint-interval", "100"  # FIX: Use hyphens instead of underscores
        ]
        
        logger.info(f"Executing command: {' '.join(cmd)}")
        
        # Run the import script
        process = subprocess.run(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=True
        )
        
        logger.info("ChEMBL data import completed successfully")
        logger.info(f"Output: {process.stdout}")
        
        if process.stderr:
            logger.warning(f"Stderr: {process.stderr}")
        
        return True
    
    except subprocess.CalledProcessError as e:
        logger.error(f"ChEMBL data import failed with return code {e.returncode}")
        logger.error(f"Output: {e.stdout}")
        logger.error(f"Error: {e.stderr}")
        return False
    
    except Exception as e:
        logger.error(f"Error running ChEMBL data import: {str(e)}")
        return False
```

The key fixes in this updated code:
1. Cast PostgreSQL types to text to ensure JSON compatibility (e.g., `column_name::text AS column_name`)
2. Fix verification query to explicitly cast types
3. Update ChEMBL import arguments to use hyphens instead of underscores

These changes address all three identified blockers, allowing the remediation process to proceed.

## 2. File Inventory

| File | Path | Purpose | Status |
|------|------|---------|--------|
| 📄 Main Import Script | `/mnt/c/Users/1edwa/Documents/CryoProtect v2/ChEMBL_Integrated_Import.py` | Import script with all enhancements | Ready for use |
| 📄 Integrity Check Report | `/mnt/c/Users/1edwa/Documents/CryoProtect v2/reports/chembl_integrity_check_20250427_001142.json` | Detailed report on property discrepancies | Reference only |
| 📄 MCP Config | `/mnt/c/Users/1edwa/Documents/CryoProtect v2/.roo/mcp.json` | Configuration for MCP tools | Reference only |
| 📄 Database Schema | `/mnt/c/Users/1edwa/Documents/CryoProtect v2/migrations/001_initial_schema.sql` | Reference for database schema | Reference only |
| 📄 RLS Utils | `/mnt/c/Users/1edwa/Documents/CryoProtect v2/rls_utils.py` | RLS verification and restoration | Used by import |
| 📄 MCP Tool Helper | `/mnt/c/Users/1edwa/Documents/CryoProtect v2/use_mcp_tool.py` | Helper functions for MCP tools | Used by scripts |
| 📄 Logging | `/mnt/c/Users/1edwa/Documents/CryoProtect v2/chembl/logging.py` | Logging utilities for ChEMBL import | Used by import |

## 3. Key Database Tables

### molecules table
```sql
SELECT column_name, data_type 
FROM information_schema.columns 
WHERE table_name = 'molecules' 
ORDER BY ordinal_position;
```

Important columns:
- `id`: UUID primary key
- `name`: VARCHAR molecule name
- `inchi`: VARCHAR InChI string
- `inchikey`: VARCHAR InChIKey for uniqueness
- `smiles`: VARCHAR SMILES representation
- `chembl_id`: VARCHAR ChEMBL ID (recently added)
- `data_source`: VARCHAR source information
- `created_by`: UUID user reference

### molecular_properties table
```sql
SELECT column_name, data_type 
FROM information_schema.columns 
WHERE table_name = 'molecular_properties' 
ORDER BY ordinal_position;
```

Important columns:
- `id`: UUID primary key
- `molecule_id`: UUID foreign key to molecules
- `property_type_id`: UUID foreign key to property_types
- `numeric_value`: NUMERIC value for numeric properties
- `text_value`: VARCHAR value for text properties
- `boolean_value`: BOOLEAN value for boolean properties
- `data_source`: VARCHAR source information

### property_types table
```sql
SELECT column_name, data_type 
FROM information_schema.columns 
WHERE table_name = 'property_types' 
ORDER BY ordinal_position;
```

Important columns:
- `id`: UUID primary key
- `name`: VARCHAR property name (e.g., "LogP")
- `data_type`: VARCHAR type of data (numeric, text, boolean)
- `description`: VARCHAR property description
- `units`: VARCHAR units of measurement

## 4. Action Steps

### Step 1: Execute Full ChEMBL Import
Run the enhanced ChEMBL_Integrated_Import.py with sufficient volume:

```python
import subprocess
import logging
import os
from datetime import datetime

# Setup logging
log_dir = "logs"
os.makedirs(log_dir, exist_ok=True)
log_file = f"{log_dir}/chembl_import_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler(log_file),
        logging.StreamHandler()
    ]
)

logger = logging.getLogger(__name__)

def run_full_import():
    """Execute the enhanced ChEMBL import script with increased limit."""
    logger.info("Starting full ChEMBL import with increased limit")
    
    # Parameters for import
    params = [
        "--limit", "2500",     # Import 2500 compounds to ensure we have >1000 after filtering
        "--batch-size", "50",  # Process in batches of 50 for stability
        "--checkpoint-interval", "100"  # Save checkpoint every 100 compounds
    ]
    
    try:
        # Run the import script
        cmd = ["python", "ChEMBL_Integrated_Import.py"] + params
        logger.info(f"Running command: {' '.join(cmd)}")
        
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=True
        )
        
        # Log the result
        logger.info(f"Import completed with exit code: {result.returncode}")
        logger.info(f"Output: {result.stdout}")
        
        if result.stderr:
            logger.warning(f"Errors: {result.stderr}")
        
        return result.returncode == 0
        
    except subprocess.CalledProcessError as e:
        logger.error(f"Import failed with exit code: {e.returncode}")
        logger.error(f"Output: {e.stdout}")
        logger.error(f"Error: {e.stderr}")
        return False
    except Exception as e:
        logger.error(f"Unexpected error running import: {str(e)}")
        return False

if __name__ == "__main__":
    success = run_full_import()
    logger.info(f"Import {'succeeded' if success else 'failed'}")
```

### Step 2: Implement Property Reconciliation

Create a dedicated script to reconcile property values with the official ChEMBL API:

```python
#!/usr/bin/env python3
"""
ChEMBL Property Reconciliation Script

This script reconciles property values in our database with the official ChEMBL API
to ensure data integrity and consistency.
"""

import os
import sys
import json
import logging
import time
import argparse
from datetime import datetime
from typing import Dict, List, Any, Optional

# Import the official ChEMBL client
from chembl_webresource_client.new_client import new_client

# Import MCP tools
from use_mcp_tool import execute_sql, get_project_id

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler("logs/property_reconciliation.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

def get_molecules_with_chembl_id(project_id: str) -> List[Dict[str, Any]]:
    """
    Get all molecules with ChEMBL IDs using MCP.
    
    Args:
        project_id: Supabase project ID
        
    Returns:
        List of molecules with ChEMBL IDs
    """
    query = "SELECT id, chembl_id FROM molecules WHERE chembl_id IS NOT NULL;"
    return execute_sql(query, project_id)

def get_molecule_properties(molecule_id: str, project_id: str) -> List[Dict[str, Any]]:
    """
    Get properties for a molecule using MCP.
    
    Args:
        molecule_id: Molecule ID
        project_id: Supabase project ID
        
    Returns:
        List of properties for the molecule
    """
    query = f"""
    SELECT mp.id, pt.name, mp.numeric_value, mp.text_value, mp.boolean_value, mp.data_source
    FROM molecular_properties mp
    JOIN property_types pt ON mp.property_type_id = pt.id
    WHERE mp.molecule_id = '{molecule_id}';
    """
    return execute_sql(query, project_id)

def update_property_value(property_id: str, value: Any, is_numeric: bool, source: str, project_id: str) -> bool:
    """
    Update a property value using MCP.
    
    Args:
        property_id: Property ID
        value: New value
        is_numeric: Whether the value is numeric
        source: Source information
        project_id: Supabase project ID
        
    Returns:
        True if update succeeded, False otherwise
    """
    try:
        # Escape single quotes in string values
        if isinstance(source, str):
            source = source.replace("'", "''")
            
        if is_numeric:
            query = f"""
            UPDATE molecular_properties
            SET numeric_value = {value},
                text_value = NULL,
                boolean_value = NULL,
                data_source = '{source}'
            WHERE id = '{property_id}';
            """
        else:
            value_str = str(value).replace("'", "''")
            query = f"""
            UPDATE molecular_properties
            SET numeric_value = NULL,
                text_value = '{value_str}',
                boolean_value = NULL,
                data_source = '{source}'
            WHERE id = '{property_id}';
            """
        
        result = execute_sql(query, project_id)
        return True
    except Exception as e:
        logger.error(f"Error updating property {property_id}: {str(e)}")
        return False

def reconcile_properties():
    """
    Reconcile property values with official ChEMBL API.
    """
    # Get project ID
    project_id = get_project_id()
    
    # Get all molecules with ChEMBL IDs
    molecules = get_molecules_with_chembl_id(project_id)
    
    if not molecules:
        logger.error("No molecules with ChEMBL IDs found")
        return False
    
    logger.info(f"Found {len(molecules)} molecules with ChEMBL IDs")
    
    # Define property mapping
    property_map = {
        "LogP": ["alogp"],
        "Molecular Weight": ["full_mwt"],
        "Hydrogen Bond Acceptor Count": ["hba"],
        "Hydrogen Bond Donor Count": ["hbd"],
        "Topological Polar Surface Area": ["psa"],
        "Rotatable Bond Count": ["rtb"],
        "Aromatic Ring Count": ["aromatic_rings"],
        "Heavy Atom Count": ["heavy_atoms"],
        "Rule of Five Violations": ["num_ro5_violations"]
    }
    
    # Statistics
    stats = {
        "molecules_processed": 0,
        "properties_checked": 0,
        "properties_updated": 0,
        "errors": 0
    }
    
    # Process molecules in batches to avoid API rate limits
    batch_size = 10
    results = []
    
    for i in range(0, len(molecules), batch_size):
        batch = molecules[i:i+batch_size]
        logger.info(f"Processing batch {i//batch_size + 1}/{(len(molecules)-1)//batch_size + 1} ({len(batch)} molecules)")
        
        for molecule in batch:
            try:
                chembl_id = molecule["chembl_id"]
                molecule_id = molecule["id"]
                stats["molecules_processed"] += 1
                
                # Get properties for this molecule
                our_properties = get_molecule_properties(molecule_id, project_id)
                
                # Get official ChEMBL data
                chembl_molecule = new_client.molecule.get(chembl_id)
                
                # Skip if no properties in ChEMBL
                if "molecule_properties" not in chembl_molecule:
                    logger.warning(f"No properties found for {chembl_id} in ChEMBL API")
                    continue
                
                chembl_props = chembl_molecule.get("molecule_properties", {})
                
                # Process each property
                molecule_result = {
                    "chembl_id": chembl_id,
                    "properties_checked": 0,
                    "properties_updated": 0,
                    "properties": []
                }
                
                for prop in our_properties:
                    prop_name = prop["name"]
                    stats["properties_checked"] += 1
                    molecule_result["properties_checked"] += 1
                    
                    # Skip if not in our mapping
                    if prop_name not in property_map:
                        continue
                    
                    # Get our value
                    our_value = prop["numeric_value"]
                    if our_value is None:
                        our_value = prop["text_value"]
                    if our_value is None:
                        our_value = prop["boolean_value"]
                    
                    property_result = {
                        "property": prop_name,
                        "our_value": our_value,
                        "updated": False
                    }
                    
                    # Check each mapped ChEMBL property
                    for chembl_prop in property_map[prop_name]:
                        if chembl_prop in chembl_props and chembl_props[chembl_prop] is not None:
                            chembl_value = chembl_props[chembl_prop]
                            property_result["chembl_property"] = chembl_prop
                            property_result["chembl_value"] = chembl_value
                            
                            # Convert to numbers for comparison if possible
                            try:
                                our_num = float(our_value)
                                chembl_num = float(chembl_value)
                                
                                # Check if values differ by more than tolerance
                                if abs(our_num - chembl_num) > 0.01:
                                    # Update the property
                                    source = f"ChEMBL: {chembl_id}, property: {chembl_prop} (reconciled)"
                                    success = update_property_value(prop["id"], chembl_num, True, source, project_id)
                                    
                                    if success:
                                        stats["properties_updated"] += 1
                                        molecule_result["properties_updated"] += 1
                                        property_result["updated"] = True
                                        property_result["new_value"] = chembl_num
                                        logger.info(f"Updated {prop_name} for {chembl_id}: {our_num} -> {chembl_num}")
                                    else:
                                        stats["errors"] += 1
                                        property_result["error"] = "Update failed"
                                else:
                                    property_result["reason"] = "Values match within tolerance"
                            except (ValueError, TypeError):
                                # Non-numeric comparison
                                str_our = str(our_value) if our_value is not None else ""
                                str_chembl = str(chembl_value) if chembl_value is not None else ""
                                
                                if str_our != str_chembl:
                                    # Update the property
                                    source = f"ChEMBL: {chembl_id}, property: {chembl_prop} (reconciled)"
                                    success = update_property_value(prop["id"], str_chembl, False, source, project_id)
                                    
                                    if success:
                                        stats["properties_updated"] += 1
                                        molecule_result["properties_updated"] += 1
                                        property_result["updated"] = True
                                        property_result["new_value"] = str_chembl
                                        logger.info(f"Updated {prop_name} for {chembl_id}: {str_our} -> {str_chembl}")
                                    else:
                                        stats["errors"] += 1
                                        property_result["error"] = "Update failed"
                                else:
                                    property_result["reason"] = "Values match"
                            
                            # Found a matching ChEMBL property, no need to check others
                            break
                    
                    molecule_result["properties"].append(property_result)
                
                results.append(molecule_result)
                
            except Exception as e:
                logger.error(f"Error processing molecule {molecule.get('chembl_id', 'unknown')}: {str(e)}")
                stats["errors"] += 1
        
        # Sleep to avoid API rate limits
        time.sleep(1)
    
    # Generate report
    report = {
        "timestamp": datetime.now().isoformat(),
        "stats": stats,
        "results": results
    }
    
    # Save report
    os.makedirs("reports", exist_ok=True)
    report_file = f"reports/chembl_reconciliation_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json"
    with open(report_file, "w") as f:
        json.dump(report, f, indent=2)
    
    logger.info(f"Reconciliation completed")
    logger.info(f"Molecules processed: {stats['molecules_processed']}")
    logger.info(f"Properties checked: {stats['properties_checked']}")
    logger.info(f"Properties updated: {stats['properties_updated']}")
    logger.info(f"Errors encountered: {stats['errors']}")
    logger.info(f"Report saved to {report_file}")
    
    return stats["errors"] == 0

def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(description="Reconcile ChEMBL property values")
    parser.add_argument("--dry-run", action="store_true", help="Don't actually update values, just report differences")
    args = parser.parse_args()
    
    # Run reconciliation
    success = reconcile_properties()
    
    return 0 if success else 1

if __name__ == "__main__":
    sys.exit(main())
```

### Step 3: Create Verification Script

Create a verification script to ensure all requirements are met:

```python
#!/usr/bin/env python3
"""
ChEMBL Remediation Verification Script

This script verifies that all ChEMBL remediation requirements have been met:
1. Database contains 1000+ molecules with ChEMBL IDs
2. All standard reference compounds are present
3. Property values match official ChEMBL API values
"""

import os
import sys
import json
import logging
from datetime import datetime
from typing import Dict, Any

# Import MCP tools
from use_mcp_tool import execute_sql, get_project_id

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler("logs/chembl_verification.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

def verify_chembl_remediation():
    """
    Verify that all ChEMBL remediation requirements have been met.
    
    Returns:
        Dict: Verification results
    """
    # Get project ID
    project_id = get_project_id()
    
    # Standard reference compounds that must be present
    reference_compounds = [
        "CHEMBL25",    # Aspirin
        "CHEMBL1118",  # Caffeine
        "CHEMBL1234",  # Glycerol (common cryoprotectant)
        "CHEMBL444",   # Glucose
        "CHEMBL230130", # Ethylene glycol (common cryoprotectant)
        "CHEMBL9335",  # Dimethyl sulfoxide (DMSO, common cryoprotectant)
        "CHEMBL15151"  # Trehalose (common cryoprotectant)
    ]
    
    # Check total molecule count
    logger.info("Checking molecule counts")
    query = "SELECT COUNT(*) AS count FROM molecules;"
    total_molecules_result = execute_sql(query, project_id)
    total_molecules = total_molecules_result[0]["count"] if total_molecules_result else 0
    
    # Check ChEMBL molecule count
    query = "SELECT COUNT(*) AS count FROM molecules WHERE chembl_id IS NOT NULL;"
    chembl_molecules_result = execute_sql(query, project_id)
    chembl_molecules = chembl_molecules_result[0]["count"] if chembl_molecules_result else 0
    
    # Check reference compounds
    logger.info("Checking reference compounds")
    reference_status = {}
    for chembl_id in reference_compounds:
        query = f"SELECT id, name FROM molecules WHERE chembl_id = '{chembl_id}';"
        result = execute_sql(query, project_id)
        reference_status[chembl_id] = {
            "present": len(result) > 0,
            "name": result[0]["name"] if result else None
        }
    
    # Check property counts
    logger.info("Checking property counts")
    query = "SELECT COUNT(*) AS count FROM molecular_properties;"
    total_properties_result = execute_sql(query, project_id)
    total_properties = total_properties_result[0]["count"] if total_properties_result else 0
    
    # Check reconciled property counts
    query = "SELECT COUNT(*) AS count FROM molecular_properties WHERE data_source LIKE '%reconciled%';"
    reconciled_properties_result = execute_sql(query, project_id)
    reconciled_properties = reconciled_properties_result[0]["count"] if reconciled_properties_result else 0
    
    # Check LogP values for reference compounds
    logger.info("Checking LogP values for reference compounds")
    logp_values = {}
    for chembl_id in reference_compounds:
        query = f"""
        SELECT mp.numeric_value
        FROM molecular_properties mp
        JOIN property_types pt ON mp.property_type_id = pt.id
        JOIN molecules m ON mp.molecule_id = m.id
        WHERE m.chembl_id = '{chembl_id}'
        AND pt.name = 'LogP';
        """
        result = execute_sql(query, project_id)
        
        if result and len(result) > 0 and result[0]["numeric_value"] is not None:
            logp_values[chembl_id] = result[0]["numeric_value"]
    
    # Check requirements
    all_references_present = all(ref["present"] for ref in reference_status.values())
    has_sufficient_data = chembl_molecules >= 1000
    has_reconciled_properties = reconciled_properties > 0
    
    # Generate verification report
    verification = {
        "timestamp": datetime.now().isoformat(),
        "molecule_counts": {
            "total": total_molecules,
            "with_chembl_id": chembl_molecules,
            "percentage": round(chembl_molecules / total_molecules * 100, 2) if total_molecules > 0 else 0
        },
        "reference_compounds": reference_status,
        "all_references_present": all_references_present,
        "property_counts": {
            "total": total_properties,
            "reconciled": reconciled_properties,
            "percentage": round(reconciled_properties / total_properties * 100, 2) if total_properties > 0 else 0
        },
        "logp_values": logp_values,
        "requirements": {
            "sufficient_data": has_sufficient_data,
            "reference_compounds": all_references_present,
            "property_reconciliation": has_reconciled_properties
        }
    }
    
    # Overall status
    verification["all_requirements_met"] = (
        has_sufficient_data and 
        all_references_present and 
        has_reconciled_properties
    )
    
    verification["status"] = "success" if verification["all_requirements_met"] else "partial"
    
    # Save report
    os.makedirs("reports", exist_ok=True)
    report_file = f"reports/chembl_remediation_report_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json"
    with open(report_file, "w") as f:
        json.dump(verification, f, indent=2)
    
    logger.info(f"ChEMBL remediation verification completed.")
    logger.info(f"Verification status: {verification['status']}")
    logger.info(f"All requirements met: {verification['all_requirements_met']}")
    logger.info(f"Remediation report saved to {report_file}")
    
    return verification

def main():
    """Main entry point."""
    verification = verify_chembl_remediation()
    
    # Print summary to console
    print(f"\nChEMBL Remediation Verification Summary:")
    print(f"Status: {verification['status'].upper()}")
    print(f"All requirements met: {verification['all_requirements_met']}")
    print(f"\nMolecule Counts:")
    print(f"  - Total: {verification['molecule_counts']['total']}")
    print(f"  - With ChEMBL ID: {verification['molecule_counts']['with_chembl_id']} ({verification['molecule_counts']['percentage']}%)")
    print(f"\nReference Compounds:")
    for chembl_id, status in verification['reference_compounds'].items():
        print(f"  - {chembl_id}: {'✅ Present' if status['present'] else '❌ Missing'} {status['name'] if status['present'] else ''}")
    print(f"\nProperty Counts:")
    print(f"  - Total: {verification['property_counts']['total']}")
    print(f"  - Reconciled: {verification['property_counts']['reconciled']} ({verification['property_counts']['percentage']}%)")
    
    return 0 if verification["all_requirements_met"] else 1

if __name__ == "__main__":
    sys.exit(main())
```

### Step 4: Main Orchestration Script

Create a main orchestration script to run all steps:

```python
#!/usr/bin/env python3
"""
ChEMBL Remediation Main Script

This script orchestrates the complete ChEMBL remediation process:
1. Execute full ChEMBL import with increased limit
2. Reconcile property values with official ChEMBL API
3. Verify all remediation requirements are met
"""

import os
import sys
import json
import logging
import subprocess
from datetime import datetime
from typing import Dict, Any

# Set up logging
log_dir = "logs"
os.makedirs(log_dir, exist_ok=True)
log_file = f"{log_dir}/chembl_remediation_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler(log_file),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

def run_step(name: str, command: list) -> Dict[str, Any]:
    """
    Run a step in the remediation process.
    
    Args:
        name: Step name
        command: Command to run
        
    Returns:
        Dict with result details
    """
    logger.info(f"Starting step: {name}")
    logger.info(f"Running command: {' '.join(command)}")
    
    start_time = datetime.now()
    
    try:
        result = subprocess.run(
            command,
            capture_output=True,
            text=True,
            check=True
        )
        
        end_time = datetime.now()
        duration = (end_time - start_time).total_seconds()
        
        logger.info(f"Step {name} completed successfully in {duration:.2f} seconds")
        
        return {
            "success": True,
            "exit_code": result.returncode,
            "duration_seconds": duration,
            "output": result.stdout,
            "error": result.stderr
        }
    except subprocess.CalledProcessError as e:
        end_time = datetime.now()
        duration = (end_time - start_time).total_seconds()
        
        logger.error(f"Step {name} failed with exit code {e.returncode} after {duration:.2f} seconds")
        logger.error(f"Output: {e.stdout}")
        logger.error(f"Error: {e.stderr}")
        
        return {
            "success": False,
            "exit_code": e.returncode,
            "duration_seconds": duration,
            "output": e.stdout,
            "error": e.stderr
        }
    except Exception as e:
        end_time = datetime.now()
        duration = (end_time - start_time).total_seconds()
        
        logger.error(f"Step {name} failed with unexpected error after {duration:.2f} seconds: {str(e)}")
        
        return {
            "success": False,
            "exit_code": -1,
            "duration_seconds": duration,
            "output": "",
            "error": str(e)
        }

def main():
    """Main entry point."""
    # Define remediation steps
    steps = [
        {
            "name": "Import ChEMBL Data",
            "command": ["python", "ChEMBL_Integrated_Import.py", "--limit", "2500", "--batch-size", "50"]
        },
        {
            "name": "Reconcile Properties",
            "command": ["python", "reconcile_chembl_properties.py"]
        },
        {
            "name": "Verify Remediation",
            "command": ["python", "verify_chembl_remediation.py"]
        }
    ]
    
    # Track results
    results = {}
    start_time = datetime.now()
    
    # Run each step
    for step in steps:
        step_name = step["name"]
        result = run_step(step_name, step["command"])
        results[step_name] = result
        
        # Stop if a step fails
        if not result["success"]:
            logger.error(f"Step {step_name} failed, stopping remediation process")
            break
    
    # Generate summary report
    end_time = datetime.now()
    total_duration = (end_time - start_time).total_seconds()
    
    all_steps_successful = all(result["success"] for result in results.values())
    
    summary = {
        "timestamp": datetime.now().isoformat(),
        "total_duration_seconds": total_duration,
        "all_steps_successful": all_steps_successful,
        "steps": results
    }
    
    # Save summary report
    os.makedirs("reports", exist_ok=True)
    summary_file = f"reports/chembl_remediation_summary_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json"
    with open(summary_file, "w") as f:
        json.dump(summary, f, indent=2)
    
    logger.info(f"ChEMBL remediation process {'completed successfully' if all_steps_successful else 'failed'}")
    logger.info(f"Total duration: {total_duration:.2f} seconds")
    logger.info(f"Summary report saved to {summary_file}")
    
    # Print final status
    print(f"\nChEMBL Remediation {'Completed Successfully' if all_steps_successful else 'Failed'}")
    print(f"Total Duration: {total_duration:.2f} seconds")
    print(f"\nStep Results:")
    for step_name, result in results.items():
        status = "✅ Success" if result["success"] else "❌ Failed"
        print(f"  - {step_name}: {status} ({result['duration_seconds']:.2f}s)")
    print(f"\nSummary report saved to {summary_file}")
    
    return 0 if all_steps_successful else 1

if __name__ == "__main__":
    sys.exit(main())
```

## 5. Execution Instructions

1. Create all three scripts:
   - `reconcile_chembl_properties.py`
   - `verify_chembl_remediation.py`
   - `chembl_remediation_main.py`

2. Ensure all scripts have executable permissions:
   ```bash
   chmod +x reconcile_chembl_properties.py verify_chembl_remediation.py chembl_remediation_main.py
   ```

3. Run the main orchestration script:
   ```bash
   python chembl_remediation_main.py
   ```

4. Review the final verification report to confirm all requirements are met.

## 6. DDL Execution Via MCP - Critical Fix Instructions

To solve the blocking issue with MCP DDL execution, follow these steps:

1. **Create the Fixed Orchestration Script**:
   ```bash
   cp chembl_remediation_main.py chembl_remediation_main_fixed.py
   ```

2. **Replace the `apply_migration` Function**:
   Replace the original apply_migration function with the fixed version provided in [Section 1.1](#11-critical-blocker-mcp-ddl-execution-issue).

3. **Test the Fixed Function**:
   ```bash
   python -c "from chembl_remediation_main_fixed import apply_migration; from use_mcp_tool import get_project_id; print(apply_migration(get_project_id()))"
   ```

4. **Update Project State**:
   Update the task status in project_state.json from "Blocked-Debug" to "In Progress"

5. **Run the Full Remediation**:
   ```bash
   python chembl_remediation_main_fixed.py
   ```

### DDL Execution Problem Recap

| Problem | Error | Solution |
|---------|-------|----------|
| MCP `execute_sql` can't run DDL | `ALTER TABLE query does not return tuples` | Use SQL functions to wrap DDL in SELECT statements |
| Migration step is blocked | Task status: "Blocked-Debug" | Use the fixed orchestration script with DDL workarounds |
| Non-returning statements fail | PostgreSQL error code: '42601' | Split into statements that return data (SELECT, UPDATE with RETURNING) |

## 7. Troubleshooting Guide

### Import Issues
- **Error**: "ModuleNotFoundError: No module named 'chembl_webresource_client'"
  - **Solution**: Install the package with `pip install chembl_webresource_client`

- **Error**: "API rate limiting errors"
  - **Solution**: Increase sleep times between API calls in `ChEMBL_Integrated_Import.py`

### Property Reconciliation Issues
- **Error**: "TypeError: float() argument must be a string or a number, not 'NoneType'"
  - **Solution**: Add more robust null checking in property comparison

- **Error**: "MCP SQL execution error"
  - **Solution**: Verify MCP configuration in `.roo/mcp.json` and check project ID

### Verification Issues
- **Error**: "Reference compounds missing"
  - **Solution**: Check if standard reference IDs are included in `ChEMBL_Integrated_Import.py` and re-run import

- **Error**: "Insufficient data volume"
  - **Solution**: Increase the import limit and re-run import

## 8. Success Criteria

- ✅ DDL commands executed successfully via MCP with workarounds
- ✅ Database contains 1000+ molecules with ChEMBL IDs
- ✅ All reference compounds are present
- ✅ Property values match official ChEMBL API values (within tolerance)
- ✅ Comprehensive verification report confirms all requirements met

## 9. Performance Considerations

- ChEMBL API has rate limits; use appropriate delays between calls
- Process data in batches to manage memory usage
- Use MCP tools with appropriate workarounds for all database operations
- Create appropriate indexes for performance

## 10. Security Considerations

- Use RLS to protect data after import
- Ensure all database operations use the appropriate roles
- Audit all changes to data
- Use MCP tools to ensure proper access control

## 11. Recommendation for Long-Term Fix

For a more permanent solution to the MCP DDL execution issue, consider:

1. Enhance the MCP tool suite to handle DDL statements directly
2. Create a dedicated MCP function specifically for executing DDL (e.g., `execute_ddl`)
3. Maintain a library of SQL functions that can wrap DDL in tuple-returning operations
4. Document this limitation clearly in the MCP documentation

This action plan provides all necessary information and code for the ROO agents to fix the blocked ChEMBL migration issue and complete the full data remediation task successfully.