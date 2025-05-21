# ROO Agent Guide: ChEMBL Integration with Direct Supabase Connection

## Overview

This document provides implementation guidance for the Master Orchestrator to coordinate the integration of ChEMBL data using a direct Supabase connection for improved performance.

## Implementation Approach

### 1. Direct Connection Module

**Task**: Implement the `supabase_direct.py` module
**Agent**: Database Connection Specialist
**Reference Files**: 
- `PubChem_CryoProtectants_Supabase_Enhanced_MCP.py` (lines 43-46)
- `connection_pool_wrapper.py`

**Implementation Guidelines**:
- Use psycopg2 for direct PostgreSQL connection
- Implement connection pooling for efficiency
- Extract credentials from environment variables
- Build singleton connection manager pattern

**Key Functions**:
- `get_connection()`: Acquire connection from pool
- `execute_query(query, params)`: Execute parameterized SQL
- `execute_batch(queries)`: Run multiple queries in transaction

**Expected Output**:
```python
# Core interface example
def execute_sql(query, params=None):
    """Execute SQL with direct connection instead of MCP."""
    conn = get_connection()
    try:
        with conn.cursor() as cursor:
            cursor.execute(query, params)
            result = cursor.fetchall() if cursor.description else None
            conn.commit()
            return result
    except Exception as e:
        conn.rollback()
        raise
    finally:
        release_connection(conn)
```

### 2. ChEMBL Client Adaptation

**Task**: Enhance the ChEMBL client module
**Agent**: API Integration Specialist
**Reference Files**:
- `chembl/client.py` (lines 141-202)
- `chembl/rate_limiter.py`

**Implementation Guidelines**:
- Add Monday-specific rate limiting (ChEMBL has stricter limits on weekdays)
- Enhance caching mechanisms for resilience
- Implement circuit breaker pattern
- Add specific data transformations for cryoprotectant properties

**Key Functions**:
- `get_molecule_by_chembl_id()`: Fetch complete molecule data
- `search_cryoprotectant_compounds()`: Find relevant molecules
- `normalize_chembl_data()`: Transform to CryoProtect schema

### 3. Database Operations

**Task**: Implement efficient database operations
**Agent**: Database Performance Specialist
**Reference Files**:
- `ChEMBL_Integrated_Import.py` (lines 257-339)
- `PubChem_CryoProtectants_Supabase_Enhanced_MCP.py` (lines 766-836)

**Implementation Guidelines**:
- Use bulk inserts instead of individual statements
- Implement proper transaction management
- Add retry logic for transient failures
- Enforce proper escaping for SQL injection prevention

**Key Functions**:
- `batch_insert_molecules(molecules_batch)`: Insert multiple molecules efficiently
- `insert_property_types(property_types)`: Insert property definitions
- `link_properties_to_molecules(properties_map)`: Create property relationships

### 4. Progress Tracking

**Task**: Implement real-time progress tracking
**Agent**: UI/Dashboard Specialist
**Reference Files**:
- `PubChem_CryoProtectants_Supabase_Enhanced_MCP.py` (lines 344-534)

**Implementation Guidelines**:
- Reuse the existing dashboard framework
- Add ChEMBL-specific metrics
- Implement checkpointing for resumable operations
- Create detailed logging for troubleshooting

**Key Functions**:
- `update_progress(current, total)`: Update progress metrics
- `log_chembl_import_status()`: Record detailed import statistics
- `generate_summary_report()`: Create final import report

## Integration Strategy

### Phased Execution

1. **Foundation (Day 1)**
   - Implement direct Supabase connection
   - Set up environment configuration
   - Test connection stability

2. **Data Acquisition (Day 1-2)**
   - Implement enhanced ChEMBL client
   - Build query strategies for cryoprotectants
   - Test with small sample set

3. **Data Processing (Day 2-3)**
   - Implement data transformation
   - Develop property mappings
   - Validate scientific accuracy

4. **Database Population (Day 3-4)**
   - Implement batch insertion
   - Add data verification
   - Test with full dataset

5. **Performance Optimization (Day 4-5)**
   - Add indexing for query performance
   - Optimize batch sizes
   - Fine-tune connection pooling

### Success Criteria

- **Volume**: 1,000+ molecules imported
- **Efficiency**: 90%+ reduction in import time vs. MCP
- **Reliability**: 99%+ success rate for operations
- **Consistency**: Data matches ChEMBL source exactly

## Code Structure Guidelines

### Modularity

- Each component should have single responsibility
- Use dependency injection for testability
- Separate core logic from database operations

### Error Handling

- Use specific exception types
- Implement retry logic for transient failures
- Log detailed error information
- Add graceful degradation

### Security

- Never hardcode credentials
- Use parameterized queries exclusively
- Escape all user inputs
- Implement proper connection closure

## Example Task Delegation

```
Task: Implement execute_batch function in supabase_direct.py
Agent: Database Connection Specialist
File: supabase_direct.py
Lines: 120-158

This function should:
1. Accept a list of SQL queries
2. Execute them in a single transaction
3. Return results for queries that return data
4. Handle errors with proper rollback
5. Use connection pooling efficiently

Reference implementation in connection_pool_wrapper.py:75-110
```

## Environment Setup

Add these variables to `.env`:

```
SUPABASE_DB_HOST=db.xxxxxxxxxxxx.supabase.co
SUPABASE_DB_PORT=5432
SUPABASE_DB_NAME=postgres
SUPABASE_DB_USER=postgres
SUPABASE_DB_PASSWORD=your-database-password
SUPABASE_DB_MAX_CONNECTIONS=10
```

## Performance Optimization Guidelines

1. **Batch Operations**
   - Insert 50-100 molecules per batch
   - Use explicit transactions
   - Minimize roundtrips to database

2. **Connection Management**
   - Reuse connections for sequential operations
   - Implement proper connection pooling
   - Close connections when idle

3. **Query Optimization**
   - Use prepared statements
   - Add proper indexes
   - Minimize data transferred

## Testing Strategy

1. **Connection Testing**:
   ```python
   def test_connection():
       """Test direct connection to Supabase."""
       result = execute_sql("SELECT current_setting('role'), current_user;")
       assert result[0]['current_user'] == 'postgres'
   ```

2. **Performance Testing**:
   ```python
   def test_performance():
       """Compare MCP vs direct connection performance."""
       start = time.time()
       # Execute test operations
       end = time.time()
       assert (end - start) < 1.0  # Should be under 1 second
   ```

## References

- Official ChEMBL Python Client: https://github.com/chembl/chembl_webresource_client
- Psycopg2 Documentation: https://www.psycopg.org/docs/
- Supabase PostgreSQL Documentation: https://supabase.com/docs/guides/database