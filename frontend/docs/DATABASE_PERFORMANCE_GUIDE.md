# Database Performance Optimization Guide

This guide explains the database performance optimization implementation for the CryoProtect application, with a focus on database indexing strategies and query optimization.

## Overview

The database performance optimization ensures efficient data access by implementing:

1. Strategic indexes on frequently queried columns
2. Specialized indexes for text search and pattern matching
3. Optimized join performance through foreign key indexes
4. Automated analysis of query patterns
5. Performance monitoring and reporting

## Key Optimizations

### Added Indexes

The implementation adds several types of indexes to improve performance:

1. **B-tree Indexes**: Standard indexes for equality and range queries
   - Primary lookup fields (IDs, foreign keys)
   - Numeric fields used in range queries
   - Date fields used for sorting

2. **Text Search Indexes**:
   - Case-insensitive search with lowercase transformations
   - Pattern matching with trigram indexes (GIN)
   - Full-text search capabilities

3. **Foreign Key Indexes**:
   - Optimized join performance
   - Efficient referential integrity checking

### Target Tables

Indexes have been added to the following key tables:

1. **molecules**: For efficient molecule lookups by ID, name, formula, type, etc.
2. **molecular_properties**: For property lookups by molecule, name, value, etc.
3. **consolidated_molecules**: For duplicate management lookups
4. **property_types**: For property metadata lookups
5. **units**: For unit conversion queries
6. **mixtures** and **mixture_components**: For mixture-related queries
7. **experiments**: For experiment data queries
8. **User-related tables**: For access control queries

## Index Creation Process

The indexing process involves:

1. **Analysis Phase**:
   - Identify frequently queried columns
   - Examine query patterns and execution plans
   - Check existing indexes to avoid duplication
   - Determine optimal index types for each column

2. **Implementation Phase**:
   - Create indexes through SQL migration
   - Verify index creation success
   - Check query plans to confirm index usage
   - Monitor performance improvements

## Usage

### Running the Index Creation

The indexing process can be run using the provided script:

```bash
# Check existing indexes (dry-run mode)
./run_add_performance_indexes.sh --dry-run

# Add indexes and simulate query improvements
./run_add_performance_indexes.sh --simulate

# Add indexes without simulation
./run_add_performance_indexes.sh
```

### Monitoring Performance

After adding indexes, you can monitor performance:

```bash
# Generate performance analysis report
python3 analyze_query_performance.py

# Run simulation tests
./run_add_performance_indexes.sh --simulate
```

## Database Maintenance

Regular maintenance ensures indexes remain efficient:

1. **Index Analysis**:
   ```sql
   -- Find unused indexes
   SELECT * FROM pg_stat_user_indexes 
   WHERE idx_scan = 0 
   AND idx_tup_read = 0 
   AND idx_tup_fetch = 0;
   
   -- Find bloated indexes
   SELECT * FROM pg_stat_user_indexes 
   ORDER BY idx_tup_read + idx_tup_fetch DESC;
   ```

2. **Index Rebuilding**:
   ```sql
   -- Rebuild an index to remove bloat
   REINDEX INDEX index_name;
   ```

3. **Table Vacuuming**:
   ```sql
   -- Vacuum a table to reclaim space
   VACUUM ANALYZE table_name;
   ```

## Query Optimization Tips

Beyond indexing, optimize queries with these approaches:

1. **Use appropriate WHERE clauses**:
   - Place the most selective conditions first
   - Use indexed columns in WHERE clauses when possible

2. **Limit result sets**:
   - Use LIMIT to restrict large result sets
   - Use paging for large data retrieval

3. **Avoid expensive operations**:
   - Use EXISTS instead of COUNT(*) when checking existence
   - Use JOINs instead of subqueries when possible
   - Avoid functions on indexed columns in WHERE clauses

4. **Monitor query plans**:
   - Use EXPLAIN ANALYZE to understand query execution
   - Look for sequential scans that should be using indexes
   - Check for unexpected nested loops or hash joins

## Troubleshooting

If performance issues persist:

1. **Verify index usage**:
   - Check query plans with EXPLAIN ANALYZE
   - Ensure indexes are being used as expected
   - Verify queries don't contain functions on indexed columns

2. **Check index health**:
   - Look for bloated indexes that need rebuilding
   - Verify tables are being properly vacuumed
   - Check for fragmentation issues

3. **Connection pooling**:
   - Verify connection pool settings are appropriate
   - Check for connection leaks or excessive connections
   - Monitor connection wait times

4. **Resource allocation**:
   - Review PostgreSQL memory settings
   - Check shared_buffers and work_mem settings
   - Monitor disk I/O for bottlenecks