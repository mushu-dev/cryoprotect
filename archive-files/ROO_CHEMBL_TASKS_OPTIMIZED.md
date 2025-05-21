# ChEMBL Integration: Optimized Tasks

## Task 1: Direct Supabase Connection

**File**: `supabase_direct.py`  
**Priority**: 1  
**Reference**: `connection_pool_wrapper.py:10-85`

**Implementation**:
```python
class SupabaseDirectConnection:
    _instance = None
    
    @classmethod
    def get_instance(cls):
        # Singleton implementation
        pass
        
    def __init__(self):
        # Initialize connection pool
        # Read from environment variables:
        # SUPABASE_DB_HOST, SUPABASE_DB_PORT, SUPABASE_DB_NAME, 
        # SUPABASE_DB_USER, SUPABASE_DB_PASSWORD
        pass
        
    def execute_query(self, query, params=None):
        # Execute single query with result conversion
        pass
        
    def execute_batch(self, queries, transaction=True):
        # Execute multiple queries in transaction
        pass
```

**Verification**:
```python
db = SupabaseDirectConnection.get_instance()
result = db.execute_query("SELECT 1 as test")
assert result[0]['test'] == 1
```

## Task 2: ChEMBL Client Update

**File**: `chembl/client.py`  
**Priority**: 2  
**Reference**: `PubChem_CryoProtectants_Supabase_Enhanced_MCP.py:850-900`

**Implementation**:
1. Add weekday-specific rate limiting (lines 72-85)
2. Update caching for resilience (lines 162-173)
3. Implement circuit breaker (lines 93-115)

**Interface**:
```python
def get_molecule_by_chembl_id(chembl_id, use_cache=True):
    # Fetch molecule with caching and rate limiting
    pass
    
def search_molecules(query, limit=10):
    # Search with caching and rate limiting
    pass
```

**Verification**:
```python
client = ResilientChEMBLClient()
aspirin = client.get_molecule_by_chembl_id("CHEMBL25")
assert aspirin.get("pref_name") == "ASPIRIN"
```

## Task 3: Data Transformation

**File**: `ChEMBL_Integrated_Import.py`  
**Priority**: 3  
**Reference**: `PubChem_CryoProtectants_Supabase_Enhanced_MCP.py:766-834`

**Implementation**:
1. Update transform_chembl_to_molecule (lines 570-615)
2. Update transform_chembl_to_properties (lines 618-807)

**Interface**:
```python
def transform_chembl_to_molecule(compound):
    # Transform ChEMBL molecule to our schema
    pass
    
def transform_chembl_to_properties(compound, molecule_id, property_type_map):
    # Transform ChEMBL properties to our schema
    pass
```

**Test Data**:
```python
# Use standard reference compound
compound = fetch_compound_by_id("CHEMBL1234")  # Glycerol
transformed = transform_chembl_to_molecule(compound)
assert transformed["name"] == "GLYCEROL"
assert "smiles" in transformed
```

## Task 4: Replace MCP with Direct Connection

**File**: `ChEMBL_Integrated_Import.py`  
**Priority**: 4  
**Reference**: `supabase_direct.py` (from Task 1)

**Changes**:
1. Replace MCP imports with direct connection
   ```python
   # Before:
   from use_mcp_tool import execute_sql, get_project_id
   
   # After:
   from supabase_direct import SupabaseDirectConnection
   ```

2. Replace SQL execution (lines 165-188, 199-220, etc.)
   ```python
   # Before:
   result = execute_sql(sql, get_project_id_for_mcp())
   
   # After:
   db = SupabaseDirectConnection.get_instance()
   result = db.execute_query(sql)
   ```

3. Replace transaction handling (lines 1000-1070)
   ```python
   # Before:
   transaction_queries = ["BEGIN;", ..., "COMMIT;"]
   combined_query = "\n".join(transaction_queries)
   execute_sql_through_mcp(combined_query)
   
   # After:
   db = SupabaseDirectConnection.get_instance()
   db.execute_batch(queries)
   ```

**Verification**:
```python
result = insert_molecule(test_molecule)
assert result is not None
```

## Task 5: Progress Tracking

**File**: `ChEMBL_Integrated_Import.py`  
**Priority**: 5  
**Reference**: `PubChem_CryoProtectants_Supabase_Enhanced_MCP.py:344-534`

**Implementation**:
1. Import required progress tracking components
   ```python
   class ProgressTracker:
       def __init__(self, total_compounds, batch_size, checkpoint_path):
           pass
       
       def start_batch(self, batch_num, batch_size):
           pass
           
       def end_batch(self, processed, imported, skipped, errors):
           pass
           
       def save_checkpoint(self):
           pass
   ```

2. Initialize and use tracker in main function
   ```python
   tracker = ProgressTracker(len(compounds), batch_size, checkpoint_path)
   tracker.start_batch(batch_num, len(batch))
   # ... process batch ...
   tracker.end_batch(processed, imported, skipped, errors)
   ```

**Verification**:
```python
tracker = ProgressTracker(100, 10, "test_checkpoint.json")
tracker.start_batch(0, 10)
tracker.end_batch(10, 8, 1, 1)
assert tracker.total_processed == 10
```

## Task 6: Data Verification

**File**: `verify_chembl_data.py`  
**Priority**: 6  
**Reference**: `verify_database_integrity.py:80-150`

**Implementation**:
```python
def verify_molecule_counts():
    """Verify molecule count matches expected range."""
    db = SupabaseDirectConnection.get_instance()
    result = db.execute_query(
        "SELECT COUNT(*) FROM molecules WHERE data_source LIKE '%ChEMBL%'"
    )
    count = result[0]['count']
    return {
        "count": count,
        "status": "SUCCESS" if count >= 1000 else "WARNING"
    }
    
def verify_property_distribution():
    """Verify property distribution across molecules."""
    # Implementation
    pass
    
def generate_verification_report():
    """Generate comprehensive verification report."""
    # Implementation
    pass
```

**Verification**:
```python
report = generate_verification_report()
assert "molecule_count" in report
assert "property_distribution" in report
```

## Environment Setup

**Required Variables**:
```
SUPABASE_DB_HOST=db.xxxxxxxxxxxx.supabase.co
SUPABASE_DB_PORT=5432
SUPABASE_DB_NAME=postgres
SUPABASE_DB_USER=postgres
SUPABASE_DB_PASSWORD=your-database-password
SUPABASE_DB_MAX_CONNECTIONS=10
```

**Implementation Order**:
1. Task 1: Direct Connection
2. Task 2: ChEMBL Client
3. Task 3: Data Transformation
4. Task 4: MCP Replacement
5. Task 5: Progress Tracking
6. Task 6: Verification