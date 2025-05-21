# RLS Query Optimization for Convex

## Overview

This document describes the implementation of optimized access control for complex queries in the CryoProtect application with Convex. The implementation provides better performance, maintainability, and reusability for access control logic.

## Key Improvements

1. **Centralized Access Control**
   - Consolidated access checking functions in one place
   - Improved maintainability and consistency
   - Better separation of concerns

2. **Batch Access Checking**
   - Optimized filtering of record IDs for access control
   - Reduced database roundtrips
   - Improved performance for complex queries

3. **Access Caching**
   - Cached access control decisions
   - Reduced redundant permission checks
   - Configurable TTL for cache entries

4. **Optimized Indexes**
   - Additional indexes for common access patterns
   - Composite indexes for efficient joins
   - Improved query performance

## Implementation Details

### Access Control Module

The `access_control.ts` module in `convex/utils` provides a comprehensive set of functions to handle access control for different entity types:

- **userAccess** - Core user permission functions
- **projectAccess** - Project-level access control
- **moleculeAccess** - Molecule-specific permissions
- **mixtureAccess** - Mixture-specific permissions
- **experimentAccess** - Experiment-specific permissions
- **modelAccess** - Scientific model permissions

Each module provides three key methods:
- `hasAccess` - Check if a user has read access to a resource
- `canModify` - Check if a user has write access to a resource
- `filterAccessible` - Filter a list of IDs to only those accessible by the user

### Optimized Indexes

The `optimized_indexes.ts` module provides additional indexes for common query patterns:

- Composite indexes for team memberships (`by_project_user`)
- Status and visibility indexes (`by_status_public`)
- Property name and value indexes (`by_property_name_numeric`)
- Parameter and value indexes (`by_parameter_numeric`)
- Relationship indexes (`by_molecule_mixture`)

### Access Caching

The `accessCache` module provides caching of access control decisions:

- TTL-based cache with automatic expiration
- Efficient key generation based on function arguments
- Automatic cleanup of expired entries
- Pre-wrapped functions for common access patterns

## Usage Examples

### Basic Access Checking

```typescript
import { userAccess } from "../utils/access_control";

export const getMolecule = query({
  args: { id: v.id("molecules") },
  handler: async (ctx, args) => {
    // Check if the current user has access to this molecule
    const userId = ctx.auth.userId;
    if (!userId) throw new Error("Authentication required");

    const hasAccess = await userAccess.hasAccess(
      ctx.db, 
      "molecules", 
      args.id, 
      userId
    );
    
    if (!hasAccess) {
      throw new Error("Access denied");
    }
    
    return await ctx.db.get(args.id);
  }
});
```

### Batch Access Filtering

```typescript
import { moleculeAccess } from "../utils/access_control";

export const getMoleculesByIds = query({
  args: { ids: v.array(v.id("molecules")) },
  handler: async (ctx, args) => {
    const userId = ctx.auth.userId;
    if (!userId) throw new Error("Authentication required");
    
    // Filter to only the molecules the user can access
    const accessibleIds = await moleculeAccess.filterAccessible(
      ctx.db,
      args.ids,
      userId
    );
    
    // Get the accessible molecules
    return await Promise.all(
      accessibleIds.map(id => ctx.db.get(id))
    );
  }
});
```

### Using Cached Access

```typescript
import { cachedAccess } from "../utils/access_control";

export const getExperiment = query({
  args: { id: v.id("experiments") },
  handler: async (ctx, args) => {
    const userId = ctx.auth.userId;
    if (!userId) throw new Error("Authentication required");

    // Use cached access check for better performance
    const hasAccess = await cachedAccess.experiment.hasAccess(
      ctx.db, 
      args.id, 
      userId
    );
    
    if (!hasAccess) {
      throw new Error("Access denied");
    }
    
    return await ctx.db.get(args.id);
  }
});
```

## Performance Benefits

The optimized access control implementation provides significant performance improvements:

1. **Reduced Database Roundtrips**
   - Batch access checking reduces the number of database operations
   - Caching avoids redundant checks for the same resources

2. **More Efficient Queries**
   - Optimized indexes improve query performance
   - Composite indexes enable more efficient filtering

3. **Improved Caching**
   - In-memory cache for frequent access patterns
   - TTL-based expiration to ensure freshness while maintaining performance

4. **Better Resource Utilization**
   - Reduced CPU usage for access checks
   - Lower memory usage for complex queries
   - Faster response times for end users

## Testing

The implementation includes a comprehensive test suite in `convex/utils/__tests__/access_control.test.ts` that verifies:

- Basic access control for different entity types
- Batch access filtering
- Cache functionality
- Access inheritance (project -> resources)
- Role-based permissions

## Best Practices for Query Optimization

1. **Use Batch Operations**
   - Always use `filterAccessible` for multiple documents
   - Avoid checking access for each document individually

2. **Leverage Caching**
   - Use `cachedAccess` for frequently accessed resources
   - Consider the caching TTL when data changes frequently

3. **Optimize Indexes**
   - Use the optimized indexes for query patterns
   - Create additional indexes for specific query patterns if needed

4. **Minimize Query Scope**
   - Filter as early as possible in the query chain
   - Only fetch the fields you need
   - Use pagination for large result sets

5. **Consider Access Patterns**
   - Group resources by project for efficient access checking
   - Use public/private flags for resources that don't need complex access control
   - Store ownership information directly on resources

## Conclusion

The optimized access control implementation for Convex provides significantly better performance for complex queries while maintaining strong security guarantees. By centralizing access control logic, implementing caching, and optimizing indexes, we've reduced query latency and improved overall application responsiveness.