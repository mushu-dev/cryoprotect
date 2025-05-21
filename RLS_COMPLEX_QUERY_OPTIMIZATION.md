# RLS Complex Query Optimization Guide

This guide explains how to optimize complex queries that use Row-Level Security (RLS) policies in the CryoProtect database. The optimization techniques focus on improving performance while maintaining security.

## Overview

Complex queries with RLS policies can be slow because:

1. RLS policy evaluation happens for every row
2. Complex conditions in policies add significant overhead
3. Joins across tables with RLS policies compound the performance impact
4. Group by and aggregate operations become more expensive

The optimization utilities in this repository address these issues through:

1. Security definer functions for common access patterns
2. Performance-tuned indexes for RLS policy expressions
3. Materialized views for frequently accessed data patterns
4. Query rewriting to use more efficient access patterns
5. Cache layer for temporary query results

## Optimization Script

The main optimization script is `optimize_complex_rls_queries.py`, which applies various RLS optimizations to the database. It can be run directly or using the convenient wrapper script `run_optimize_complex_rls_queries.sh`.

### Prerequisites

- Python 3.6+
- For direct database connection: psycopg2
- For Supabase: supabase-py
- For MCP: MCP-provided Supabase tools

### Usage

```bash
# Using the wrapper script (recommended)
./run_optimize_complex_rls_queries.sh --supabase --test-performance

# Direct usage
python optimize_complex_rls_queries.py --direct --db-host localhost --db-name cryoprotect_db --db-user postgres --db-password your_password --test-performance
```

### Connection Options

The script supports three connection modes:

1. **Direct**: Connect directly to PostgreSQL
2. **Supabase**: Use Supabase client
3. **MCP**: Use Supabase MCP for Claude Code

### Optimization Options

- `--verify-only`: Only verify optimizations without applying them
- `--test-performance`: Test query performance after optimization
- `--skip-functions`: Skip creating security definer functions
- `--skip-indexes`: Skip creating performance indexes
- `--skip-views`: Skip creating materialized views
- `--skip-policies`: Skip creating RLS policies
- `--skip-complex-queries`: Skip complex query optimizations
- `--dry-run`: Show what would be done without making changes
- `--force`: Apply optimizations even if they already exist

## Testing Performance

A separate script, `test_rls_complex_queries.py`, allows you to test the performance of complex queries before and after optimization.

```bash
# Test query performance
python test_rls_complex_queries.py --supabase --iterations 5 --output-report performance_report.md
```

This script will:

1. Run a set of complex queries with the original SQL
2. Run the same queries using optimized SQL
3. Compare execution times and generate a detailed report

## Optimization Techniques

### 1. Security Definer Functions

Security definer functions run with elevated privileges but enforce security checks internally. This reduces the overhead of repetitive RLS policy evaluations.

Example:
```sql
CREATE OR REPLACE FUNCTION has_molecule_access(molecule_id uuid)
RETURNS boolean AS $$
BEGIN
  RETURN EXISTS (
    SELECT 1 FROM molecules m
    WHERE m.id = molecule_id AND (
      m.is_public = true OR
      m.created_by = auth.uid() OR
      EXISTS (
        SELECT 1 FROM project_molecules pm
        JOIN team_projects tp ON pm.project_id = tp.project_id
        JOIN user_profile up ON tp.team_id = up.team_id
        WHERE pm.molecule_id = molecule_id AND up.auth_user_id = auth.uid()
      )
    )
  );
END;
$$ LANGUAGE plpgsql SECURITY DEFINER;
```

### 2. Performance Indexes

Specialized indexes support efficient execution of RLS policies and complex queries:

- Indexes on foreign keys referenced in joins
- Indexes for boolean flags like `is_public`
- Indexes for user IDs and ownership
- Partial indexes for common cases like public data
- Text search indexes for name and description fields
- Covering indexes for common query patterns

### 3. Materialized Views

Materialized views pre-compute complex joins and aggregations:

```sql
CREATE MATERIALIZED VIEW public_molecules_summary AS
SELECT 
    m.id, 
    m.name, 
    m.molecular_formula, 
    m.smiles, 
    m.cid, 
    m.pubchem_link,
    m.is_public,
    COUNT(mp.id) AS property_count
FROM 
    molecules m
LEFT JOIN 
    molecular_properties mp ON m.id = mp.molecule_id
WHERE 
    m.is_public = true
GROUP BY 
    m.id;
```

### 4. Query Rewriting

Complex queries are rewritten to use security definer functions:

Before:
```sql
SELECT 
    m.id,
    m.name,
    m.smiles,
    m.molecular_formula,
    COUNT(mp.id) AS property_count
FROM 
    molecules m
LEFT JOIN
    molecular_properties mp ON m.id = mp.molecule_id
WHERE 
    m.is_public = true 
    OR m.created_by = auth.uid()
    OR EXISTS (
        SELECT 1 FROM project_molecules pm
        JOIN team_projects tp ON pm.project_id = tp.project_id
        JOIN user_profile up ON tp.team_id = up.team_id
        WHERE pm.molecule_id = m.id AND up.auth_user_id = auth.uid()
    )
GROUP BY
    m.id, m.name, m.smiles, m.molecular_formula
```

After:
```sql
SELECT * FROM get_molecules_with_properties();
```

### 5. Query Result Caching

A caching layer for temporary query results:

```sql
-- Create a cache table
CREATE TABLE query_cache (
    cache_key TEXT PRIMARY KEY,
    result JSONB NOT NULL,
    created_at TIMESTAMP WITH TIME ZONE DEFAULT CURRENT_TIMESTAMP,
    expires_at TIMESTAMP WITH TIME ZONE NOT NULL
);

-- Cache utility function
CREATE OR REPLACE FUNCTION cache_query_result(
    query_id TEXT,
    result JSONB,
    ttl_minutes INTEGER DEFAULT 10
) RETURNS void AS $$
BEGIN
    DELETE FROM query_cache WHERE cache_key = query_id;
    
    INSERT INTO query_cache (cache_key, result, expires_at)
    VALUES (query_id, result, CURRENT_TIMESTAMP + (ttl_minutes || ' minutes')::INTERVAL);
END;
$$ LANGUAGE plpgsql;
```

## Performance Expectations

Optimizations typically provide:

- 50-95% reduction in query execution time for complex queries
- 10-30x speedup for materialized view access
- Consistent sub-100ms response times for common queries
- Significantly reduced database CPU utilization
- Better scalability with growing dataset size

## Best Practices

1. **Use security definer functions**: Replace complex RLS policy logic with function calls
2. **Use materialized views for public data**: Pre-compute commonly accessed public data
3. **Refresh materialized views regularly**: Set up a schedule to keep views updated
4. **Add indexes for RLS policy columns**: Ensure all columns used in policy conditions are indexed
5. **Monitor query performance**: Regularly check for slow queries and optimize them
6. **Use query result caching**: Cache results for expensive queries that don't change frequently
7. **Optimize joins across tables with RLS**: Use specialized functions for cross-table queries

## Troubleshooting

Common issues:

1. **Slow materialized view refresh**: Update scheduled refresh timing or use concurrent refresh
2. **Missing indexes**: Check explain plans to identify missing indexes for RLS policies
3. **Stale cache results**: Reduce cache TTL or add targeted cache invalidation
4. **Security definer function errors**: Ensure functions properly validate user access
5. **Query optimizer not using indexes**: Update database statistics or rewrite queries

## Additional Resources

- [PostgreSQL RLS Documentation](https://www.postgresql.org/docs/current/ddl-rowsecurity.html)
- [PostgreSQL Function Performance](https://www.postgresql.org/docs/current/xfunc-optimization.html)
- [PostgreSQL Index Types](https://www.postgresql.org/docs/current/indexes-types.html)
- [PostgreSQL Materialized Views](https://www.postgresql.org/docs/current/rules-materializedviews.html)