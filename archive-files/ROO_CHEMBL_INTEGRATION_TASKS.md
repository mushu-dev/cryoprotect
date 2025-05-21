# ChEMBL Integration: Task Breakdown for ROO Orchestrator

## Overview

This document provides a structured breakdown of specific implementation tasks for the Master Orchestrator. Each task includes precise file references, implementation details, and clear success criteria.

## Task 1: Direct Supabase Connection Implementation

**File**: `supabase_direct.py` (create new)
**Priority**: Highest (Blocking)
**Estimated Completion**: 2 hours

**Subtasks**:
1. Create connection pool implementation using psycopg2
   - Lines: 1-100
   - Reference: `connection_pool_wrapper.py:10-85`
   
2. Implement query execution functions
   - Lines: 100-200
   - Reference: `PubChem_CryoProtectants_Supabase_Enhanced_MCP.py:766-834`
   
3. Add transaction management
   - Lines: 200-250
   - Reference: `ChEMBL_Integrated_Import.py:257-339`

**Success Criteria**:
- Connection pool maintains stable connections
- Queries execute successfully with parameterized inputs
- Transactions commit or roll back appropriately
- Connection credentials are properly secured

**Testing Commands**:
```python
# Test connection
import supabase_direct
db = supabase_direct.SupabaseDirectConnection.get_instance()
result = db.execute_query("SELECT current_setting('role'), current_user;")
print(result)
```

## Task 2: Enhanced ChEMBL Client

**File**: `chembl/client.py`
**Priority**: High
**Estimated Completion**: 4 hours

**Subtasks**:
1. Update rate limiter for weekday restrictions
   - Lines: 72-85
   - Reference: `chembl/rate_limiter.py:30-65`
   
2. Enhance caching mechanism
   - Lines: 162-173
   - Reference: `chembl/cache.py:15-40`
   
3. Add adaptive request mechanisms
   - Lines: 93-115
   - Reference: `PubChem_CryoProtectants_Supabase_Enhanced_MCP.py:850-900`

**Success Criteria**:
- API requests respect ChEMBL rate limits
- Caching reduces redundant API calls
- Connection degradation is handled gracefully
- Request failures are retried with backoff

**Testing Commands**:
```python
# Test ChEMBL client
from chembl.client import ResilientChEMBLClient
client = ResilientChEMBLClient()
result = client.get_molecule_by_chembl_id("CHEMBL25")  # Aspirin
print(f"Retrieved: {result.get('pref_name')}")
```

## Task 3: ChEMBL Data Transformation

**File**: `ChEMBL_Integrated_Import.py`
**Priority**: High
**Estimated Completion**: 3 hours

**Subtasks**:
1. Update the `transform_chembl_to_molecule` function
   - Lines: 570-615
   - Reference: `PubChem_CryoProtectants_Supabase_Enhanced_MCP.py:766-834`
   
2. Update the `transform_chembl_to_properties` function
   - Lines: 618-807
   - Reference: `PubChem_CryoProtectants_Supabase_Enhanced_MCP.py:766-834`
   
3. Modify the `fetch_cryoprotectant_compounds` function
   - Lines: 363-569
   - Reference: `PubChem_CryoProtectants_Supabase_Enhanced_MCP.py:582-640`

**Success Criteria**:
- ChEMBL molecule data is properly transformed to match schema
- Property mappings are correctly established
- Compound filtering identifies valid cryoprotectants
- Data formats are consistently handled

**Testing Procedure**:
1. Fetch sample compound (e.g., glycerol - CHEMBL1234)
2. Transform to our schema format
3. Verify all required fields are populated correctly

## Task 4: Database Operations Optimization

**File**: `ChEMBL_Integrated_Import.py`
**Priority**: Medium
**Estimated Completion**: 5 hours

**Subtasks**:
1. Replace MCP SQL execution with direct connection
   - Lines: 165-188, 199-220, 234-255, 290-309, 311-359
   - Reference: `supabase_direct.py` (from Task 1)
   
2. Implement batch processing
   - Lines: 810-990
   - Reference: `PubChem_CryoProtectants_Supabase_Enhanced_MCP.py:942-978`
   
3. Add transaction management
   - Lines: 1000-1070
   - Reference: `PubChem_CryoProtectants_Supabase_Enhanced_MCP.py:850-900`

**Success Criteria**:
- Direct SQL execution replaces all MCP calls
- Batch processing improves insertion performance
- Transactions ensure data consistency
- Error handling maintains database integrity

**Testing Procedure**:
1. Compare execution time between MCP and direct connection
2. Verify transaction rollback on simulated error
3. Test insertion of 100+ compounds in batches

## Task 5: Progress Tracking Implementation

**File**: `ChEMBL_Integrated_Import.py`
**Priority**: Medium
**Estimated Completion**: 4 hours

**Subtasks**:
1. Integrate progress tracker
   - Lines: NEW (add after line 60)
   - Reference: `PubChem_CryoProtectants_Supabase_Enhanced_MCP.py:344-534`
   
2. Add checkpoint mechanism
   - Lines: 404-409, 502-512, 558-567
   - Reference: `PubChem_CryoProtectants_Supabase_Enhanced_MCP.py:370-410`
   
3. Implement dashboard server
   - Lines: NEW (add after progress tracker)
   - Reference: `PubChem_CryoProtectants_Supabase_Enhanced_MCP.py:536-581`

**Success Criteria**:
- Progress is tracked and displayed in real-time
- Checkpoints allow resumable operations
- Dashboard provides visual progress feedback
- Import can be paused and resumed without data loss

## Task 6: Verification and Testing

**File**: `verify_chembl_data.py` (create new)
**Priority**: Low
**Estimated Completion**: 3 hours

**Subtasks**:
1. Implement data verification functions
   - Reference: `verify_database_integrity.py:80-150`
   
2. Add statistical analysis
   - Reference: `PubChem_CryoProtectants_Supabase_Enhanced_MCP.py:230-280`
   
3. Create validation report generator
   - Reference: `chembl/logging.py:100-160`

**Success Criteria**:
- Data integrity is verified after import
- Statistical analysis confirms expected patterns
- Validation report provides comprehensive metrics
- Any issues are clearly identified for remediation

## Task 7: Environment Configuration

**File**: `.env` and `config.py`
**Priority**: Highest (Blocking)
**Estimated Completion**: 1 hour

**Subtasks**:
1. Add Supabase direct connection variables to `.env.template`
2. Update `config.py` to load new variables
3. Document environment setup in README

**Required Variables**:
```
SUPABASE_DB_HOST=db.xxxxxxxxxxxx.supabase.co
SUPABASE_DB_PORT=5432
SUPABASE_DB_NAME=postgres
SUPABASE_DB_USER=postgres
SUPABASE_DB_PASSWORD=your-database-password
SUPABASE_DB_MAX_CONNECTIONS=10
```

**Success Criteria**:
- All required environment variables are documented
- Configuration is loaded correctly in application
- Sensitive credentials are properly protected
- Documentation provides clear setup instructions

## Integration Strategy

### Phase 1: Foundation (Day 1)
- Complete Task 7: Environment Configuration
- Complete Task 1: Direct Supabase Connection

### Phase 2: Core Components (Day 2)
- Complete Task 2: Enhanced ChEMBL Client
- Complete Task 3: ChEMBL Data Transformation

### Phase 3: Data Pipeline (Day 3-4)
- Complete Task 4: Database Operations Optimization
- Complete Task 5: Progress Tracking Implementation

### Phase 4: Verification (Day 5)
- Complete Task 6: Verification and Testing
- Perform end-to-end testing and optimization

## Communication Protocol

1. Task Assignment:
   - Master Orchestrator assigns tasks to specialized agents
   - Each assignment includes file paths, line numbers, and success criteria

2. Status Reporting:
   - Agents report completion percentage every 30 minutes
   - Blockers are reported immediately
   - Completed tasks include tests validating success criteria

3. Integration Points:
   - Task 1 → Task 4 (Database Operations depend on Direct Connection)
   - Task 2 → Task 3 (Data Transformation depends on ChEMBL Client)
   - Task 3 → Task 4 (Database Operations use transformed data)
   - All Tasks → Task 6 (Verification uses all components)

## Reference Files

- `PubChem_CryoProtectants_Supabase_Enhanced_MCP.py`: Complete PubChem integration
- `connection_pool_wrapper.py`: Connection pooling implementation
- `ChEMBL_Integrated_Import.py`: Current ChEMBL integration (MCP-based)
- `chembl/client.py`: ChEMBL API client
- `chembl/rate_limiter.py`: Rate limiting implementation
- `chembl/cache.py`: Caching mechanism
- `verify_database_integrity.py`: Database verification utilities