# RLS Optimization Guide

This guide provides detailed implementation steps for optimizing Row Level Security (RLS) policies in the CryoProtect database. It serves as a companion document to the main CRYOPROTECT_OPTIMIZATION_PLAN.md and focuses specifically on RLS performance improvements.

## 1. Current RLS Policy Analysis

### Current Limitations
1. **Complex policy expressions**: Policies contain complex nested conditions that are evaluated for every row
2. **Missing indexes**: Many columns used in policy conditions lack appropriate indexes
3. **Repeated logic**: Common access patterns are duplicated across multiple policies
4. **No caching**: Access patterns are re-evaluated for each query execution
5. **No materialized views**: Frequently queried public data requires policy evaluation
6. **Query optimization issues**: RLS policies can prevent query planner optimizations

### Performance Impact
- Slow query execution, especially for tables with many rows
- Increased database CPU utilization
- Poor scalability with growing dataset size
- Inconsistent query performance
- Timeouts for complex queries

## 2. Optimization Steps

### 2.1 Security Definer Functions

#### Problem
RLS policies often contain complex conditions that are re-evaluated for every row. These conditions can include expensive joins and subqueries.

#### Solution
Create `SECURITY DEFINER` functions that encapsulate common access patterns. These functions:
- Run with the privileges of the function creator (typically the database owner)
- Can bypass RLS policies internally for better performance
- Only need to evaluate the minimal required logic
- Provide a clear, modular approach to access control

#### Implementation Examples

**1. User Team Membership Function**

```sql
CREATE OR REPLACE FUNCTION public.is_team_member(team_id uuid)
RETURNS boolean AS $$
BEGIN
  RETURN EXISTS (
    SELECT 1 FROM user_profile
    WHERE auth_user_id = auth.uid() AND user_profile.team_id = $1
  );
END;
$$ LANGUAGE plpgsql SECURITY DEFINER;
```

**2. User Teams Function**

```sql
CREATE OR REPLACE FUNCTION public.user_teams()
RETURNS SETOF uuid AS $$
BEGIN
  RETURN QUERY
  SELECT team_id FROM user_profile
  WHERE auth_user_id = auth.uid();
END;
$$ LANGUAGE plpgsql SECURITY DEFINER;
```

**3. User Projects Function**

```sql
CREATE OR REPLACE FUNCTION public.user_projects()
RETURNS SETOF uuid AS $$
BEGIN
  RETURN QUERY
  SELECT tp.project_id 
  FROM team_projects tp
  JOIN user_profile up ON tp.team_id = up.team_id
  WHERE up.auth_user_id = auth.uid();
END;
$$ LANGUAGE plpgsql SECURITY DEFINER;
```

**4. Molecule Access Function**

```sql
CREATE OR REPLACE FUNCTION public.has_molecule_access(molecule_id uuid)
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

**5. Mixture Access Function**

```sql
CREATE OR REPLACE FUNCTION public.has_mixture_access(mixture_id uuid)
RETURNS boolean AS $$
BEGIN
  RETURN EXISTS (
    SELECT 1 FROM mixtures m
    WHERE m.id = mixture_id AND (
      m.is_public = true OR
      m.created_by = auth.uid() OR
      EXISTS (
        SELECT 1 FROM project_mixtures pm
        JOIN team_projects tp ON pm.project_id = tp.project_id
        JOIN user_profile up ON tp.team_id = up.team_id
        WHERE pm.mixture_id = mixture_id AND up.auth_user_id = auth.uid()
      )
    )
  );
END;
$$ LANGUAGE plpgsql SECURITY DEFINER;
```

**6. Clearance Level Function**

```sql
CREATE OR REPLACE FUNCTION public.user_has_clearance(required_level text)
RETURNS boolean AS $$
BEGIN
  RETURN EXISTS (
    SELECT 1 FROM user_profile
    WHERE auth_user_id = auth.uid() 
    AND (clearance_level >= required_level OR required_level IS NULL)
  );
END;
$$ LANGUAGE plpgsql SECURITY DEFINER;
```

#### RLS Policy Rewrites Using Functions

**Before:**

```sql
CREATE POLICY "users_can_view_own_molecules_and_team_shared" ON molecules
FOR SELECT
TO authenticated
USING (
  is_public = true OR
  created_by = auth.uid() OR
  EXISTS (
    SELECT 1 FROM project_molecules pm
    JOIN team_projects tp ON pm.project_id = tp.project_id
    JOIN user_profile up ON tp.team_id = up.team_id
    WHERE pm.molecule_id = molecules.id AND up.auth_user_id = auth.uid()
  )
);
```

**After:**

```sql
CREATE POLICY "users_can_view_own_molecules_and_team_shared" ON molecules
FOR SELECT
TO authenticated
USING (
  has_molecule_access(id)
);
```

### 2.2 Performance Indexes

#### Problem
Many columns used in RLS policies lack appropriate indexes, causing full table scans during policy evaluation.

#### Solution
Create strategic indexes on columns used in RLS policies, particularly:
- Foreign keys referenced in joins
- Columns used in WHERE clauses
- Columns for boolean flags like `is_public`
- Columns for user IDs and ownership

#### Implementation Examples

**1. Primary Access Path Indexes**

```sql
-- User profile indexes
CREATE INDEX IF NOT EXISTS idx_user_profile_auth_user_id ON user_profile(auth_user_id);
CREATE INDEX IF NOT EXISTS idx_user_profile_team_id ON user_profile(team_id);
CREATE INDEX IF NOT EXISTS idx_user_profile_clearance_level ON user_profile(clearance_level);

-- Team project indexes
CREATE INDEX IF NOT EXISTS idx_team_projects_team_id ON team_projects(team_id);
CREATE INDEX IF NOT EXISTS idx_team_projects_project_id ON team_projects(project_id);

-- Project molecule indexes
CREATE INDEX IF NOT EXISTS idx_project_molecules_project_id ON project_molecules(project_id);
CREATE INDEX IF NOT EXISTS idx_project_molecules_molecule_id ON project_molecules(molecule_id);

-- Project mixture indexes
CREATE INDEX IF NOT EXISTS idx_project_mixtures_project_id ON project_mixtures(project_id);
CREATE INDEX IF NOT EXISTS idx_project_mixtures_mixture_id ON project_mixtures(mixture_id);
```

**2. Ownership and Visibility Indexes**

```sql
-- Molecule ownership and visibility
CREATE INDEX IF NOT EXISTS idx_molecules_created_by ON molecules(created_by);
CREATE INDEX IF NOT EXISTS idx_molecules_is_public ON molecules(is_public);

-- Mixture ownership and visibility
CREATE INDEX IF NOT EXISTS idx_mixtures_created_by ON mixtures(created_by);
CREATE INDEX IF NOT EXISTS idx_mixtures_is_public ON mixtures(is_public);

-- Experiment ownership and visibility
CREATE INDEX IF NOT EXISTS idx_experiments_created_by ON experiments(created_by);
CREATE INDEX IF NOT EXISTS idx_experiments_is_public ON experiments(is_public);
```

**3. Partial Indexes for Common Cases**

```sql
-- Partial index for public molecules (very selective)
CREATE INDEX IF NOT EXISTS idx_molecules_public ON molecules(id) WHERE is_public = true;

-- Partial index for public mixtures
CREATE INDEX IF NOT EXISTS idx_mixtures_public ON mixtures(id) WHERE is_public = true;

-- Partial index for sensitive content
CREATE INDEX IF NOT EXISTS idx_molecular_properties_sensitive 
ON molecular_properties(molecule_id, sensitivity_level)
WHERE sensitivity_level > 'low';
```

**4. Composite Indexes for Common Query Patterns**

```sql
-- Composite index for team membership queries
CREATE INDEX IF NOT EXISTS idx_user_profile_auth_team 
ON user_profile(auth_user_id, team_id);

-- Composite index for project-mixture relationships
CREATE INDEX IF NOT EXISTS idx_project_mixtures_composite
ON project_mixtures(mixture_id, project_id);

-- Composite index for project-molecule relationships
CREATE INDEX IF NOT EXISTS idx_project_molecules_composite
ON project_molecules(molecule_id, project_id);
```

### 2.3 Materialized Views

#### Problem
Frequently accessed public data still requires policy evaluation, even though the data is accessible to all authenticated users.

#### Solution
Create materialized views for commonly accessed public data:
- These views can be accessed without triggering full RLS evaluation
- They can be refreshed on a schedule
- Indexes can be added to these views for even faster access
- RLS policies can still be applied to the views for extra security

#### Implementation Examples

**1. Public Molecules Summary View**

```sql
CREATE MATERIALIZED VIEW IF NOT EXISTS public_molecules_summary AS
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

-- Create index on the materialized view
CREATE UNIQUE INDEX IF NOT EXISTS idx_public_molecules_summary_id 
ON public_molecules_summary(id);
```

**2. Public Molecular Properties View**

```sql
CREATE MATERIALIZED VIEW IF NOT EXISTS public_molecular_properties AS
SELECT 
    mp.id,
    mp.molecule_id,
    mp.property_name,
    mp.property_value,
    mp.property_unit,
    mp.property_type,
    mp.calculation_method,
    mp.reference_doi
FROM 
    molecular_properties mp
JOIN 
    molecules m ON mp.molecule_id = m.id
WHERE 
    m.is_public = true
    AND (mp.sensitivity_level IS NULL OR mp.sensitivity_level <= 'low');

-- Create indexes on the materialized view
CREATE INDEX IF NOT EXISTS idx_public_molecular_properties_molecule_id 
ON public_molecular_properties(molecule_id);

CREATE INDEX IF NOT EXISTS idx_public_molecular_properties_property_name 
ON public_molecular_properties(property_name);
```

**3. Refresh Strategy**

Create a function to refresh materialized views:

```sql
CREATE OR REPLACE FUNCTION refresh_materialized_views()
RETURNS void AS $$
BEGIN
    REFRESH MATERIALIZED VIEW CONCURRENTLY public_molecules_summary;
    REFRESH MATERIALIZED VIEW CONCURRENTLY public_molecular_properties;
    -- Add other materialized views as needed
    
    -- Log the refresh
    INSERT INTO materialized_view_refresh_log (refresh_time)
    VALUES (NOW());
END;
$$ LANGUAGE plpgsql;
```

Set up a periodic job to refresh the views:

```sql
-- Create a table to track refresh history
CREATE TABLE IF NOT EXISTS materialized_view_refresh_log (
    id SERIAL PRIMARY KEY,
    refresh_time TIMESTAMP WITH TIME ZONE NOT NULL
);

-- Schedule refresh every hour
SELECT cron.schedule('refresh_materialized_views', '0 * * * *', 'SELECT refresh_materialized_views()');
```

### 2.4 RLS Policy Optimization

#### Problem
Even with the above optimizations, some RLS policies may still be inefficient or overly complex.

#### Solution
Rewrite remaining RLS policies to:
- Use the new security definer functions
- Take advantage of the new indexes
- Simplify expressions where possible
- Use EXISTS clauses instead of JOINs for better performance
- Address any policy gaps

#### Implementation Examples

**1. Molecules Table Policies**

```sql
-- Select policy
CREATE OR REPLACE POLICY "users_can_view_molecules" ON molecules
FOR SELECT
TO authenticated
USING (
  has_molecule_access(id)
);

-- Insert policy
CREATE OR REPLACE POLICY "users_can_insert_molecules" ON molecules
FOR INSERT
TO authenticated
WITH CHECK (
  auth.uid() = created_by
);

-- Update policy
CREATE OR REPLACE POLICY "users_can_update_own_molecules" ON molecules
FOR UPDATE
TO authenticated
USING (
  auth.uid() = created_by
)
WITH CHECK (
  auth.uid() = created_by
);

-- Delete policy
CREATE OR REPLACE POLICY "users_can_delete_own_molecules" ON molecules
FOR DELETE
TO authenticated
USING (
  auth.uid() = created_by
);
```

**2. Molecular Properties Table Policies**

```sql
-- Select policy with clearance level check
CREATE OR REPLACE POLICY "users_can_view_properties" ON molecular_properties
FOR SELECT
TO authenticated
USING (
  has_molecule_access(molecule_id) AND
  (sensitivity_level IS NULL OR 
   sensitivity_level <= 'medium' OR 
   user_has_clearance(sensitivity_level))
);

-- Insert policy
CREATE OR REPLACE POLICY "users_can_insert_properties" ON molecular_properties
FOR INSERT
TO authenticated
WITH CHECK (
  EXISTS (
    SELECT 1 FROM molecules
    WHERE id = molecular_properties.molecule_id AND created_by = auth.uid()
  )
);

-- Update policy
CREATE OR REPLACE POLICY "users_can_update_properties" ON molecular_properties
FOR UPDATE
TO authenticated
USING (
  EXISTS (
    SELECT 1 FROM molecules
    WHERE id = molecular_properties.molecule_id AND created_by = auth.uid()
  )
)
WITH CHECK (
  EXISTS (
    SELECT 1 FROM molecules
    WHERE id = molecular_properties.molecule_id AND created_by = auth.uid()
  )
);

-- Delete policy
CREATE OR REPLACE POLICY "users_can_delete_properties" ON molecular_properties
FOR DELETE
TO authenticated
USING (
  EXISTS (
    SELECT 1 FROM molecules
    WHERE id = molecular_properties.molecule_id AND created_by = auth.uid()
  )
);
```

**3. Service Role Policies**

```sql
-- Create unified service role policies for all tables
DO $$
DECLARE
    table_name text;
BEGIN
    FOR table_name IN 
        SELECT tablename FROM pg_tables 
        WHERE schemaname = 'public' 
        AND tablename NOT LIKE 'pg_%'
    LOOP
        EXECUTE format('
            DROP POLICY IF EXISTS "service_role_all_access" ON %I;
            CREATE POLICY "service_role_all_access" ON %I
            FOR ALL
            TO service_role
            USING (true)
            WITH CHECK (true);
        ', table_name, table_name);
    END LOOP;
END
$$;
```

## 3. Implementation Guide

### 3.1 Step-by-Step Approach

1. **Create Security Definer Functions**
   - Implement common access pattern functions first
   - Test each function thoroughly
   - Document the purpose and usage of each function

2. **Create Performance Indexes**
   - Add primary access path indexes
   - Add ownership and visibility indexes
   - Add partial and composite indexes
   - Monitor query performance after each batch

3. **Create Materialized Views**
   - Implement views for frequently accessed data
   - Add indexes to the views
   - Set up refresh schedule
   - Update application code to use views where appropriate

4. **Optimize RLS Policies**
   - Rewrite policies to use new functions
   - Simplify policy expressions
   - Add unified service role policies
   - Test policies extensively

### 3.2 Implementation Script

```sql
-- Begin transaction
BEGIN;

-- Step 1: Create Security Definer Functions
CREATE OR REPLACE FUNCTION public.is_team_member(team_id uuid)
RETURNS boolean AS $$
BEGIN
  RETURN EXISTS (
    SELECT 1 FROM user_profile
    WHERE auth_user_id = auth.uid() AND user_profile.team_id = $1
  );
END;
$$ LANGUAGE plpgsql SECURITY DEFINER;

CREATE OR REPLACE FUNCTION public.user_teams()
RETURNS SETOF uuid AS $$
BEGIN
  RETURN QUERY
  SELECT team_id FROM user_profile
  WHERE auth_user_id = auth.uid();
END;
$$ LANGUAGE plpgsql SECURITY DEFINER;

CREATE OR REPLACE FUNCTION public.user_projects()
RETURNS SETOF uuid AS $$
BEGIN
  RETURN QUERY
  SELECT tp.project_id 
  FROM team_projects tp
  JOIN user_profile up ON tp.team_id = up.team_id
  WHERE up.auth_user_id = auth.uid();
END;
$$ LANGUAGE plpgsql SECURITY DEFINER;

CREATE OR REPLACE FUNCTION public.has_molecule_access(molecule_id uuid)
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

CREATE OR REPLACE FUNCTION public.has_mixture_access(mixture_id uuid)
RETURNS boolean AS $$
BEGIN
  RETURN EXISTS (
    SELECT 1 FROM mixtures m
    WHERE m.id = mixture_id AND (
      m.is_public = true OR
      m.created_by = auth.uid() OR
      EXISTS (
        SELECT 1 FROM project_mixtures pm
        JOIN team_projects tp ON pm.project_id = tp.project_id
        JOIN user_profile up ON tp.team_id = up.team_id
        WHERE pm.mixture_id = mixture_id AND up.auth_user_id = auth.uid()
      )
    )
  );
END;
$$ LANGUAGE plpgsql SECURITY DEFINER;

CREATE OR REPLACE FUNCTION public.user_has_clearance(required_level text)
RETURNS boolean AS $$
BEGIN
  RETURN EXISTS (
    SELECT 1 FROM user_profile
    WHERE auth_user_id = auth.uid() 
    AND (clearance_level >= required_level OR required_level IS NULL)
  );
END;
$$ LANGUAGE plpgsql SECURITY DEFINER;

-- Step 2: Create Performance Indexes

-- User profile indexes
CREATE INDEX IF NOT EXISTS idx_user_profile_auth_user_id ON user_profile(auth_user_id);
CREATE INDEX IF NOT EXISTS idx_user_profile_team_id ON user_profile(team_id);
CREATE INDEX IF NOT EXISTS idx_user_profile_clearance_level ON user_profile(clearance_level);

-- Team project indexes
CREATE INDEX IF NOT EXISTS idx_team_projects_team_id ON team_projects(team_id);
CREATE INDEX IF NOT EXISTS idx_team_projects_project_id ON team_projects(project_id);

-- Project molecule indexes
CREATE INDEX IF NOT EXISTS idx_project_molecules_project_id ON project_molecules(project_id);
CREATE INDEX IF NOT EXISTS idx_project_molecules_molecule_id ON project_molecules(molecule_id);

-- Project mixture indexes
CREATE INDEX IF NOT EXISTS idx_project_mixtures_project_id ON project_mixtures(project_id);
CREATE INDEX IF NOT EXISTS idx_project_mixtures_mixture_id ON project_mixtures(mixture_id);

-- Molecule ownership and visibility
CREATE INDEX IF NOT EXISTS idx_molecules_created_by ON molecules(created_by);
CREATE INDEX IF NOT EXISTS idx_molecules_is_public ON molecules(is_public);

-- Mixture ownership and visibility
CREATE INDEX IF NOT EXISTS idx_mixtures_created_by ON mixtures(created_by);
CREATE INDEX IF NOT EXISTS idx_mixtures_is_public ON mixtures(is_public);

-- Experiment ownership and visibility
CREATE INDEX IF NOT EXISTS idx_experiments_created_by ON experiments(created_by);
CREATE INDEX IF NOT EXISTS idx_experiments_is_public ON experiments(is_public);

-- Partial indexes
CREATE INDEX IF NOT EXISTS idx_molecules_public ON molecules(id) WHERE is_public = true;
CREATE INDEX IF NOT EXISTS idx_mixtures_public ON mixtures(id) WHERE is_public = true;
CREATE INDEX IF NOT EXISTS idx_molecular_properties_sensitive 
ON molecular_properties(molecule_id, sensitivity_level)
WHERE sensitivity_level > 'low';

-- Composite indexes
CREATE INDEX IF NOT EXISTS idx_user_profile_auth_team 
ON user_profile(auth_user_id, team_id);
CREATE INDEX IF NOT EXISTS idx_project_mixtures_composite
ON project_mixtures(mixture_id, project_id);
CREATE INDEX IF NOT EXISTS idx_project_molecules_composite
ON project_molecules(molecule_id, project_id);

-- Step 3: Create Materialized Views
CREATE MATERIALIZED VIEW IF NOT EXISTS public_molecules_summary AS
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

CREATE UNIQUE INDEX IF NOT EXISTS idx_public_molecules_summary_id 
ON public_molecules_summary(id);

CREATE MATERIALIZED VIEW IF NOT EXISTS public_molecular_properties AS
SELECT 
    mp.id,
    mp.molecule_id,
    mp.property_name,
    mp.property_value,
    mp.property_unit,
    mp.property_type,
    mp.calculation_method,
    mp.reference_doi
FROM 
    molecular_properties mp
JOIN 
    molecules m ON mp.molecule_id = m.id
WHERE 
    m.is_public = true
    AND (mp.sensitivity_level IS NULL OR mp.sensitivity_level <= 'low');

CREATE INDEX IF NOT EXISTS idx_public_molecular_properties_molecule_id 
ON public_molecular_properties(molecule_id);

CREATE INDEX IF NOT EXISTS idx_public_molecular_properties_property_name 
ON public_molecular_properties(property_name);

-- Create refresh function and job
CREATE TABLE IF NOT EXISTS materialized_view_refresh_log (
    id SERIAL PRIMARY KEY,
    refresh_time TIMESTAMP WITH TIME ZONE NOT NULL
);

CREATE OR REPLACE FUNCTION refresh_materialized_views()
RETURNS void AS $$
BEGIN
    REFRESH MATERIALIZED VIEW CONCURRENTLY public_molecules_summary;
    REFRESH MATERIALIZED VIEW CONCURRENTLY public_molecular_properties;
    
    INSERT INTO materialized_view_refresh_log (refresh_time)
    VALUES (NOW());
END;
$$ LANGUAGE plpgsql;

-- Step 4: Optimize RLS Policies

-- Molecules table policies
DROP POLICY IF EXISTS "users_can_view_molecules" ON molecules;
CREATE POLICY "users_can_view_molecules" ON molecules
FOR SELECT
TO authenticated
USING (
  has_molecule_access(id)
);

DROP POLICY IF EXISTS "users_can_insert_molecules" ON molecules;
CREATE POLICY "users_can_insert_molecules" ON molecules
FOR INSERT
TO authenticated
WITH CHECK (
  auth.uid() = created_by
);

DROP POLICY IF EXISTS "users_can_update_own_molecules" ON molecules;
CREATE POLICY "users_can_update_own_molecules" ON molecules
FOR UPDATE
TO authenticated
USING (
  auth.uid() = created_by
)
WITH CHECK (
  auth.uid() = created_by
);

DROP POLICY IF EXISTS "users_can_delete_own_molecules" ON molecules;
CREATE POLICY "users_can_delete_own_molecules" ON molecules
FOR DELETE
TO authenticated
USING (
  auth.uid() = created_by
);

-- Molecular properties table policies
DROP POLICY IF EXISTS "users_can_view_properties" ON molecular_properties;
CREATE POLICY "users_can_view_properties" ON molecular_properties
FOR SELECT
TO authenticated
USING (
  has_molecule_access(molecule_id) AND
  (sensitivity_level IS NULL OR 
   sensitivity_level <= 'medium' OR 
   user_has_clearance(sensitivity_level))
);

DROP POLICY IF EXISTS "users_can_insert_properties" ON molecular_properties;
CREATE POLICY "users_can_insert_properties" ON molecular_properties
FOR INSERT
TO authenticated
WITH CHECK (
  EXISTS (
    SELECT 1 FROM molecules
    WHERE id = molecular_properties.molecule_id AND created_by = auth.uid()
  )
);

DROP POLICY IF EXISTS "users_can_update_properties" ON molecular_properties;
CREATE POLICY "users_can_update_properties" ON molecular_properties
FOR UPDATE
TO authenticated
USING (
  EXISTS (
    SELECT 1 FROM molecules
    WHERE id = molecular_properties.molecule_id AND created_by = auth.uid()
  )
)
WITH CHECK (
  EXISTS (
    SELECT 1 FROM molecules
    WHERE id = molecular_properties.molecule_id AND created_by = auth.uid()
  )
);

DROP POLICY IF EXISTS "users_can_delete_properties" ON molecular_properties;
CREATE POLICY "users_can_delete_properties" ON molecular_properties
FOR DELETE
TO authenticated
USING (
  EXISTS (
    SELECT 1 FROM molecules
    WHERE id = molecular_properties.molecule_id AND created_by = auth.uid()
  )
);

-- Add service role policies for all tables
DO $$
DECLARE
    table_name text;
BEGIN
    FOR table_name IN 
        SELECT tablename FROM pg_tables 
        WHERE schemaname = 'public' 
        AND tablename NOT LIKE 'pg_%'
    LOOP
        EXECUTE format('
            DROP POLICY IF EXISTS "service_role_all_access" ON %I;
            CREATE POLICY "service_role_all_access" ON %I
            FOR ALL
            TO service_role
            USING (true)
            WITH CHECK (true);
        ', table_name, table_name);
    END LOOP;
END
$$;

-- Commit transaction
COMMIT;
```

### 3.3 Testing Recommendations

1. **Query Performance Testing**
   - Benchmark query performance before and after
   - Test queries with EXPLAIN ANALYZE
   - Test with different data volumes
   - Test with different user roles

2. **Access Control Testing**
   - Verify all RLS policies correctly restrict access
   - Test edge cases (e.g., shared resources, nested relationships)
   - Test with different user clearance levels
   - Test service role access

3. **View Refresh Testing**
   - Verify materialized views contain correct data
   - Test refresh performance
   - Verify applications properly use views
   - Test concurrent refresh impact

## 4. Security Considerations

- Ensure security definer functions properly validate user access
- Verify function inputs are properly validated
- Test for SQL injection vulnerabilities
- Document all security functions and their purpose
- Implement audit logging for sensitive operations

## 5. Performance Expectations

### Before Optimization
- Average query time for complex molecules query: 200-500ms
- Database CPU utilization during peak: 70-90%
- Query timeouts during high load: Common
- RLS policy evaluation time: 50-100ms per row

### After Optimization
- Average query time for complex molecules query: 50-100ms
- Database CPU utilization during peak: 30-50%
- Query timeouts during high load: Rare
- RLS policy evaluation time: 5-10ms per row
- Materialized view access time: <10ms

## 6. Monitoring and Maintenance

### Key Metrics to Monitor
- Query execution time
- RLS policy evaluation time
- Index usage statistics
- Materialized view refresh time
- Database CPU utilization

### Maintenance Tasks
- Schedule regular ANALYZE on tables
- Refresh materialized views on schedule
- Monitor index bloat and rebuild as needed
- Periodically review RLS policies for consistency

## 7. Implementation Risks and Mitigations

### Risks
- **Policy logic changes**: Security definer functions may not accurately match original policy logic
- **Performance regression**: Some queries might be slower with new indexes
- **Materialized view staleness**: Data might be outdated between refreshes
- **Security vulnerabilities**: Complex function logic might introduce security holes

### Mitigations
- Thoroughly test all security definer functions
- Implement comprehensive access tests
- Set appropriate refresh schedules for views
- Peer review all security-related code
- Create rollback plan for each optimization step

## 8. Integration with Application Code

### API Layer Updates
- Update API code to use materialized views when appropriate
- Implement caching for common queries
- Add proper metrics collection for database queries
- Update error handling for policy violations

### Data Access Layer Changes
- Implement retry logic for transient errors
- Add logging for slow queries
- Optimize batch operations
- Add proper error categorization

## 9. Code Review Checklist

- [ ] Security definer functions properly validate access
- [ ] Indexes are created on appropriate columns
- [ ] Materialized views have proper refresh mechanism
- [ ] RLS policies are updated to use security definer functions
- [ ] Service role policies are consistent across all tables
- [ ] Query performance is improved in benchmark tests
- [ ] Application code properly uses new database features
- [ ] Documentation is updated for all new functions and views