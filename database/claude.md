# Database Module Reference

## Overview

The database module centralizes all database-related functionality for CryoProtect, featuring a flexible adapter-based architecture that supports multiple database backends, connection pooling, migrations, and comprehensive verification tools.

## Architecture

### Adapter Pattern Implementation
- **Core Interface**: `DatabaseAdapter` in `adapter.py`
- **Implementations**:
  - `SupabaseAdapter`: For Supabase cloud database
  - `MCPAdapter`: For MCP-enabled direct database access
  - `LocalAdapter`: For local PostgreSQL development

### Connection Management
- **Connection Pooling**: Managed through `connection_manager.py`
- **Retry Logic**: Automatic retry with exponential backoff
- **Error Handling**: Standardized error translation
- **Circuit Breaking**: Prevents cascading failures

## Key Components

### Migrations
- **Path**: `database/migrations/`
- **Version Control**: Sequential version system (001, 002, etc.)
- **CLI**: Command-line interface for applying/rolling back
- **Idempotency**: Migrations track applied state to prevent duplication

### Models
- **Path**: `database/models/`
- **Table Definitions**: Database schema representation
- **Validation Rules**: Data integrity constraints
- **Relationships**: Cross-table relationship definitions

### Population
- **Path**: `database/population/`
- **Data Sources**: PubChem, ChEMBL, reference compounds
- **Checkpointing**: Resume-capable data imports
- **Validation**: Pre and post-import validation

### Verification
- **Path**: `database/verification/`
- **Schema Checks**: Validates database schema matches expectations
- **Data Integrity**: Validates relationships and constraints
- **Performance Tests**: Query performance validation

### Utilities
- **Path**: `database/utils/`
- **Health Checks**: Database health monitoring
- **Query Builders**: SQL generation helpers
- **Performance Analysis**: Query plan analysis

## RLS Policies

The database implements Row Level Security (RLS) for fine-grained access control:

1. **Authentication-Based**: Records are associated with specific users
2. **Role-Based**: Different access levels based on user roles
3. **Organization-Based**: Data segregation by organization
4. **Service Role Bypass**: Administrative access for system operations

## Optimization Techniques

1. **Materialized Views**: For expensive calculations
2. **Performance Indexes**: Strategic indexing for common queries
3. **Query Optimization**: Carefully crafted queries to minimize execution time
4. **Connection Pooling**: Reuse connections to reduce overhead
5. **Security Definer Functions**: Optimized RLS policy execution

## Best Practices

1. **Always Use Transactions**: Especially for multi-step operations
2. **Close Connections**: Return connections to the pool promptly
3. **Parameterized Queries**: Never use string concatenation for queries
4. **Use the Adapter Pattern**: Don't bypass the database adapter
5. **Test Migrations**: Always verify migrations before production
6. **Validate Data**: Use schema validation before insertion

## Common Pitfalls

1. **Connection Leaks**: Failing to close connections
2. **RLS Policy Conflicts**: Contradictory policies causing lockout
3. **Transaction Timeouts**: Long-running transactions
4. **Index Bloat**: Too many indexes slowing writes
5. **Timestamp Issues**: Timezone inconsistencies
6. **Cursor Management**: Not properly handling large result sets