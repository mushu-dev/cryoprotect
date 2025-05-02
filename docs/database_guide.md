# Database Guide

This document describes the database architecture and operations for CryoProtect v2.

## Database Architecture

CryoProtect v2 uses a flexible database approach that supports:

1. **Local PostgreSQL**: For development and testing
2. **Supabase PostgreSQL**: For staging and production
3. **MCP API**: As a fallback mechanism

### Schema Overview

The core tables in the database are:

- `molecules`: Basic molecular information
- `molecular_properties`: Properties of molecules in a normalized structure
- `property_types`: Types of molecular properties
- `mixtures`: Cryoprotectant mixtures
- `mixture_components`: Components of mixtures with concentrations

## Database Adapters

The application uses a database adapter pattern to abstract connection details:

1. **DatabaseAdapter**: Abstract interface for all adapters
2. **LocalPostgreSQLAdapter**: Connects to local PostgreSQL
3. **SupabaseDirectAdapter**: Connects directly to Supabase PostgreSQL
4. **MCPAdapter**: Connects via MCP API

## Connection Management

The `ConnectionManager` handles database connections with features:

1. **Connection Pooling**: Efficiently manages database connections
2. **Fallback Mechanism**: Automatically tries alternative connection methods
3. **Transaction Support**: Manages database transactions
4. **Retry Logic**: Handles transient errors

## Local Database Setup

To set up the local database:

1. Install PostgreSQL
2. Configure environment variables:
   ```
   LOCAL_DB_HOST=localhost
   LOCAL_DB_PORT=5432
   LOCAL_DB_NAME=cryoprotect
   LOCAL_DB_USER=postgres
   LOCAL_DB_PASSWORD=your-password
   ```
3. Run the initialization script:
   ```bash
   python database/init_local_db.py
   ```

## Database Migrations

Database schema changes are managed through migrations:

1. Migration files are stored in `migrations/`
2. The format is `NNN_description.sql`
3. Migrations are applied in numerical order

To apply migrations:

```bash
python database/migrations/runner.py
```

## Data Population

To populate the database with scientific data:

1. Configure database connection
2. Run the population scripts:
   ```bash
   # For reference compounds
   python import_reference_compounds.py
   
   # For ChEMBL data
   python import_full_chembl.py
   
   # For PubChem data
   python PubChem_CryoProtectants_Supabase.py
   
   # For property enhancement
   python enhance_pubchem_properties.py
   ```

## Verification

To verify the database population:

```bash
python verify_imported_data.py
```

This will check:
- Molecule counts (>=500 expected)
- Property completeness
- Reference compound presence
- Cross-referencing
- Query performance

## Backup and Restore

### Creating a Backup

```bash
# Full backup
python database/utils/backup.py --full

# Partial backup (molecules only)
python database/utils/backup.py --tables molecules,molecular_properties
```

### Restoring from Backup

```bash
python database/utils/restore.py --file backup_20250430_120000.sql
```

## Performance Optimization

1. **Indexes**: Added to frequently queried columns
2. **Connection Pooling**: Reduces connection overhead
3. **Query Optimization**: Complex queries are optimized

## Troubleshooting

### Common Issues

1. **Connection Failures**
   - Check credentials in `.env`
   - Verify network connectivity
   - Check if PostgreSQL is running

2. **DNS Resolution Issues**
   - Use IP address fallback: Set `SUPABASE_DB_IP_ADDRESS`

3. **Performance Issues**
   - Check query plans: `EXPLAIN ANALYZE SELECT...`
   - Verify indexes

## Database Utilities

The `database/utils.py` module provides helper functions:

- `execute_query`: Execute SQL query
- `execute_batch`: Execute multiple queries
- `get_molecule_by_id`: Get molecule by ID
- `get_molecule_properties`: Get properties for a molecule
- `insert_molecule`: Insert new molecule
- `set_molecule_property`: Set molecule property
- `test_database_connection`: Test connection

## Best Practices

1. Use connection pool for multiple operations
2. Use transactions for related operations
3. Handle connection errors gracefully
4. Use parameterized queries to prevent SQL injection
5. Close connections when done