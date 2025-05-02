# CryoProtect v2 Database Performance Optimization

This document provides information about the database performance optimization implemented for CryoProtect v2, including the indexes added and how to apply the migration.

## Overview

Based on performance testing results, several database indexes have been added to improve query performance. These optimizations target the most frequently executed queries and address the performance bottlenecks identified during testing.

## Implemented Optimizations

### 1. Composite Index on Predictions Table

```sql
CREATE INDEX idx_predictions_mixture_property ON public.predictions(mixture_id, property_type_id);
```

**Purpose:** This composite index significantly improves performance for queries that filter on both `mixture_id` and `property_type_id`, which are common in the application. This is particularly beneficial for the "compare_predictions_experiments" query, which joins predictions and experiments on these columns.

**Expected Improvement:** 30-50% reduction in query times for prediction-related queries.

### 2. Text Search Index for Molecule Names

```sql
CREATE EXTENSION IF NOT EXISTS pg_trgm;
CREATE INDEX idx_molecule_name_trgm ON public.molecules USING gin (name gin_trgm_ops);
```

**Purpose:** This GIN index with the pg_trgm extension optimizes text search operations on molecule names, particularly for ILIKE queries. This improves performance when users search for molecules by name.

**Expected Improvement:** 70-90% reduction in response times for text search operations.

### 3. Text Search Index for Mixture Names

```sql
CREATE INDEX idx_mixture_name_trgm ON public.mixtures USING gin (name gin_trgm_ops);
```

**Purpose:** Similar to the molecule name index, this optimizes text search operations on mixture names.

**Expected Improvement:** 70-90% reduction in response times for mixture name searches.

## Applying the Migration

The performance indexes can be applied using the provided scripts:

### Windows

```
apply_performance_indexes.bat
```

### Linux/macOS

```
chmod +x apply_performance_indexes.sh
./apply_performance_indexes.sh
```

These scripts will:
1. Check if Node.js is installed
2. Run the migration script
3. Apply the indexes to your Supabase database

### Manual Application

If you prefer to apply the migration manually:

1. Go to the Supabase dashboard
2. Select your project
3. Go to SQL Editor
4. Copy the contents of `migrations/010_performance_indexes.sql`
5. Paste into the SQL Editor and run the query

## Verification

After applying the indexes, you can verify they were created successfully by:

1. Going to the Supabase dashboard
2. Selecting your project
3. Going to Database > Tables
4. Selecting a table (e.g., predictions)
5. Checking the "Indexes" tab

You should see the new indexes listed.

## Performance Impact

The added indexes are expected to improve performance in the following areas:

1. **Prediction Queries:** Faster retrieval of predictions for specific mixtures and property types
2. **Text Search:** Significantly faster searches for molecules and mixtures by name
3. **Join Operations:** Improved performance for queries that join predictions with other tables

## Monitoring

After applying these indexes, it's recommended to:

1. Run the performance tests again to measure the improvement
2. Monitor database performance in production
3. Check for any unexpected side effects (e.g., slower write operations)

## Additional Considerations

- **Index Maintenance:** These indexes will be automatically maintained by PostgreSQL
- **Storage Impact:** The indexes will increase the database size slightly
- **Write Performance:** There might be a small impact on write performance, but it should be negligible compared to the read performance benefits

## Future Optimizations

Consider implementing the following additional optimizations if needed:

1. **Query Caching:** Add caching for frequently accessed data
2. **Connection Pooling:** Configure proper connection pooling to handle concurrent users
3. **Database Configuration:** Adjust PostgreSQL configuration parameters for better performance