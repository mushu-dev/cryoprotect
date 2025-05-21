# ChEMBL Integration: Task Delegation

Master Orchestrator,

Below are the specific implementation tasks for ChEMBL integration with direct Supabase connection. Assign these tasks to appropriate ROO agents, focusing solely on implementation with minimal token usage.

## Task 1: Direct Supabase Connection

**IMPLEMENT**: PostgreSQL connection pool with direct Supabase access
**FILE**: Create `supabase_direct.py`
**REFERENCE**: `connection_pool_wrapper.py:10-85` and `PubChem_CryoProtectants_Supabase_Enhanced_MCP.py:43-46`
**VERIFY**: Connection successfully executes `SELECT 1` query

Key requirements:
- Read credentials from environment variables
- Implement connection pooling
- Add query execution with parameterization
- Include transaction management
- Follow singleton pattern from reference

## Task 2: ChEMBL Client Enhancement

**IMPLEMENT**: Enhanced ChEMBL client with Monday-specific rate limiting
**FILE**: `chembl/client.py`
**LINES**: 72-85, 162-173
**REFERENCE**: `PubChem_CryoProtectants_Supabase_Enhanced_MCP.py:65-71`
**VERIFY**: Client successfully retrieves reference compound data with rate limiting

Key requirements:
- Implement weekday-specific rate limiting
- Add resilient caching mechanism
- Include circuit breaker pattern for API failures
- Ensure Monday performance does not trigger API limits

## Task 3: Data Transformation Implementation

**IMPLEMENT**: ChEMBL data transformation to CryoProtect schema
**FILE**: `ChEMBL_Integrated_Import.py`
**LINES**: 570-615, 618-807
**REFERENCE**: `PubChem_CryoProtectants_Supabase_Enhanced_MCP.py:766-834`
**VERIFY**: Successfully transforms standard reference compounds (glycerol, DMSO)

Key requirements:
- Map ChEMBL molecule data to our schema
- Transform ChEMBL properties to our property schema
- Handle missing or non-standard data gracefully
- Include data source attribution

## Task 4: MCP Replacement with Direct Connection

**IMPLEMENT**: Replace MCP calls with direct Supabase connection
**FILE**: `ChEMBL_Integrated_Import.py`
**LINES**: 165-188, 199-220, 234-255, 290-309, 311-359
**REFERENCE**: Direct connection from Task 1
**VERIFY**: Database operations complete successfully with direct connection

Key requirements:
- Replace all MCP imports and function calls
- Use direct connection for all SQL operations
- Implement proper error handling
- Ensure transaction integrity

## Task 5: Batch Processing Implementation

**IMPLEMENT**: Efficient batch processing for ChEMBL data
**FILE**: `ChEMBL_Integrated_Import.py`
**LINES**: 810-990
**REFERENCE**: `PubChem_CryoProtectants_Supabase_Enhanced_MCP.py:942-978`
**VERIFY**: Successfully processes 100+ compounds in batches

Key requirements:
- Process compounds in configurable batches
- Implement transaction-based batch insertion
- Add error handling with partial success tracking
- Optimize batch size for performance

## Task 6: Progress Tracking

**IMPLEMENT**: Reuse progress tracking from PubChem integration
**FILE**: `ChEMBL_Integrated_Import.py`
**REFERENCE**: `PubChem_CryoProtectants_Supabase_Enhanced_MCP.py:344-534`
**VERIFY**: Displays accurate progress during import with resumability

Key requirements:
- Adapt existing tracking system for ChEMBL
- Add checkpoint-based resumable operations
- Include dashboard for visual progress monitoring
- Track success/failure/skipped statistics

## Environment Setup

The following environment variables must be set for direct connection:

```
SUPABASE_DB_HOST=db.xxxxxxxxxxxx.supabase.co
SUPABASE_DB_PORT=5432
SUPABASE_DB_NAME=postgres
SUPABASE_DB_USER=postgres
SUPABASE_DB_PASSWORD=your-database-password
SUPABASE_DB_MAX_CONNECTIONS=10
```

## Implementation Guidelines

1. Assign tasks sequentially as dependencies require
2. Focus agents solely on implementation, not explanation
3. Reference existing code rather than creating from scratch
4. Verify each component before proceeding to dependent tasks
5. Optimize for token efficiency in all communications

Report task completion status only, without explanations unless there are specific implementation issues requiring guidance.

Proceed with implementation immediately.