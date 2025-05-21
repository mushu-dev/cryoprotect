# RLS Practical Implementation for CryoProtect

This document provides a pragmatic approach to implementing RLS optimizations for the CryoProtect project. Based on the analysis of existing RLS policies in `006_rls_policies.sql` and the optimization recommendations in `RLS_OPTIMIZATION_GUIDE.md`, this plan focuses on practical solutions that provide real benefits without unnecessary complexity.

## 1. Current Implementation Analysis

The current RLS implementation has the following characteristics:

- Policies are defined for major tables (molecule, mixture, experiment, etc.)
- Each table has four policies (SELECT, INSERT, UPDATE, DELETE)
- Policies use similar patterns for project membership verification
- Many policies use complex nested EXISTS statements and multiple joins
- No security definer functions are used for reusable access patterns
- No specific indexes are created for optimizing RLS policy evaluation
- No materialized views for frequently accessed public data

## 2. Practical Implementation Steps

### 2.1 Security Definer Functions

**Problem**: We have repeated access patterns across many policies that result in duplicated logic and inefficient queries.

**Solution**: Create a minimal set of security definer functions that encapsulate common access patterns:

1. **Project Membership Function**

```sql
CREATE OR REPLACE FUNCTION public.is_project_member(project_id uuid)
RETURNS boolean AS $$
BEGIN
  RETURN EXISTS (
    SELECT 1 FROM user_profile
    WHERE user_profile.project_id = is_project_member.project_id
      AND user_profile.user_id = auth.uid()
  );
END;
$$ LANGUAGE plpgsql SECURITY DEFINER;
```

2. **User Projects Function**

```sql
CREATE OR REPLACE FUNCTION public.user_projects()
RETURNS SETOF uuid AS $$
BEGIN
  RETURN QUERY
  SELECT project_id FROM user_profile
  WHERE user_profile.user_id = auth.uid();
END;
$$ LANGUAGE plpgsql SECURITY DEFINER;
```

3. **Team Membership Function**

```sql
CREATE OR REPLACE FUNCTION public.is_team_member(team_id uuid)
RETURNS boolean AS $$
BEGIN
  RETURN EXISTS (
    SELECT 1 FROM user_profile
    WHERE user_profile.team_id = is_team_member.team_id
      AND user_profile.user_id = auth.uid()
  );
END;
$$ LANGUAGE plpgsql SECURITY DEFINER;
```

4. **Is Project Owner Function**

```sql
CREATE OR REPLACE FUNCTION public.is_project_owner(project_id uuid)
RETURNS boolean AS $$
BEGIN
  RETURN EXISTS (
    SELECT 1 FROM user_profile
    WHERE user_profile.project_id = is_project_owner.project_id
      AND user_profile.user_id = auth.uid()
      AND user_profile.role = 'owner'
  );
END;
$$ LANGUAGE plpgsql SECURITY DEFINER;
```

5. **Molecule Project Membership Function**

```sql
CREATE OR REPLACE FUNCTION public.molecule_in_user_project(molecule_id uuid)
RETURNS boolean AS $$
BEGIN
  RETURN EXISTS (
    SELECT 1 FROM molecule
    JOIN user_profile ON user_profile.project_id = molecule.project_id
    WHERE molecule.id = molecule_in_user_project.molecule_id
      AND user_profile.user_id = auth.uid()
  );
END;
$$ LANGUAGE plpgsql SECURITY DEFINER;
```

6. **Mixture Project Membership Function**

```sql
CREATE OR REPLACE FUNCTION public.mixture_in_user_project(mixture_id uuid)
RETURNS boolean AS $$
BEGIN
  RETURN EXISTS (
    SELECT 1 FROM mixture
    JOIN user_profile ON user_profile.project_id = mixture.project_id
    WHERE mixture.id = mixture_in_user_project.mixture_id
      AND user_profile.user_id = auth.uid()
  );
END;
$$ LANGUAGE plpgsql SECURITY DEFINER;
```

### 2.2 Performance Indexes

**Problem**: RLS policy evaluation requires multiple joins on columns that may not be indexed.

**Solution**: Create a focused set of indexes for columns used in RLS policies:

```sql
-- User profile indexes
CREATE INDEX IF NOT EXISTS idx_user_profile_user_id ON user_profile(user_id);
CREATE INDEX IF NOT EXISTS idx_user_profile_project_id ON user_profile(project_id);
CREATE INDEX IF NOT EXISTS idx_user_profile_team_id ON user_profile(team_id);
CREATE INDEX IF NOT EXISTS idx_user_profile_role ON user_profile(role);

-- Molecule indexes
CREATE INDEX IF NOT EXISTS idx_molecule_project_id ON molecule(project_id);

-- Mixture indexes
CREATE INDEX IF NOT EXISTS idx_mixture_project_id ON mixture(project_id);

-- Mixture component indexes
CREATE INDEX IF NOT EXISTS idx_mixture_component_mixture_id ON mixture_component(mixture_id);

-- Molecular property indexes
CREATE INDEX IF NOT EXISTS idx_molecular_property_molecule_id ON molecular_property(molecule_id);

-- Experiment indexes
CREATE INDEX IF NOT EXISTS idx_experiment_project_id ON experiment(project_id);

-- Prediction indexes
CREATE INDEX IF NOT EXISTS idx_prediction_molecule_id ON prediction(molecule_id);

-- Experiment property indexes
CREATE INDEX IF NOT EXISTS idx_experiment_property_experiment_id ON experiment_property(experiment_id);
```

### 2.3 RLS Policy Optimization

**Problem**: Current policies use complex expressions with multiple EXISTS clauses and joins.

**Solution**: Simplify policies by using the security definer functions:

**1. Molecule Table Policies**

```sql
-- Select policy
DROP POLICY IF EXISTS "Select molecules for project members" ON molecule;
CREATE POLICY "Select molecules for project members"
  ON molecule
  FOR SELECT
  USING (is_project_member(project_id));

-- Insert policy
DROP POLICY IF EXISTS "Insert molecules for project members" ON molecule;
CREATE POLICY "Insert molecules for project members"
  ON molecule
  FOR INSERT
  WITH CHECK (is_project_member(project_id));

-- Update policy
DROP POLICY IF EXISTS "Update molecules for project members" ON molecule;
CREATE POLICY "Update molecules for project members"
  ON molecule
  FOR UPDATE
  USING (is_project_member(project_id));

-- Delete policy
DROP POLICY IF EXISTS "Delete molecules for project members" ON molecule;
CREATE POLICY "Delete molecules for project members"
  ON molecule
  FOR DELETE
  USING (is_project_member(project_id));
```

**2. Mixture Table Policies**

```sql
-- Select policy
DROP POLICY IF EXISTS "Select mixtures for project members" ON mixture;
CREATE POLICY "Select mixtures for project members"
  ON mixture
  FOR SELECT
  USING (is_project_member(project_id));

-- Insert policy
DROP POLICY IF EXISTS "Insert mixtures for project members" ON mixture;
CREATE POLICY "Insert mixtures for project members"
  ON mixture
  FOR INSERT
  WITH CHECK (is_project_member(project_id));

-- Update policy
DROP POLICY IF EXISTS "Update mixtures for project members" ON mixture;
CREATE POLICY "Update mixtures for project members"
  ON mixture
  FOR UPDATE
  USING (is_project_member(project_id));

-- Delete policy
DROP POLICY IF EXISTS "Delete mixtures for project members" ON mixture;
CREATE POLICY "Delete mixtures for project members"
  ON mixture
  FOR DELETE
  USING (is_project_member(project_id));
```

**3. Mixture Component Table Policies**

```sql
-- Select policy
DROP POLICY IF EXISTS "Select mixture_components for project members" ON mixture_component;
CREATE POLICY "Select mixture_components for project members"
  ON mixture_component
  FOR SELECT
  USING (mixture_in_user_project(mixture_id));

-- Insert policy
DROP POLICY IF EXISTS "Insert mixture_components for project members" ON mixture_component;
CREATE POLICY "Insert mixture_components for project members"
  ON mixture_component
  FOR INSERT
  WITH CHECK (mixture_in_user_project(mixture_id));

-- Update policy
DROP POLICY IF EXISTS "Update mixture_components for project members" ON mixture_component;
CREATE POLICY "Update mixture_components for project members"
  ON mixture_component
  FOR UPDATE
  USING (mixture_in_user_project(mixture_id));

-- Delete policy
DROP POLICY IF EXISTS "Delete mixture_components for project members" ON mixture_component;
CREATE POLICY "Delete mixture_components for project members"
  ON mixture_component
  FOR DELETE
  USING (mixture_in_user_project(mixture_id));
```

**4. Other Table Policies**

Apply the same pattern to the remaining tables:
- experiment
- molecular_property
- prediction
- experiment_property
- calculation_method
- project
- team
- user_profile

### 2.4 Service Role Policies

**Problem**: Service role access requires explicit policies for each table.

**Solution**: Create a unified approach for service role access:

```sql
-- Function to create service role policies for all tables
CREATE OR REPLACE FUNCTION create_service_role_policies()
RETURNS void AS $$
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
$$ LANGUAGE plpgsql;

-- Execute the function
SELECT create_service_role_policies();
```

### 2.5 Materialized Views for Public Data (Optional)

If the application has a significant amount of public reference data that is frequently accessed, create materialized views:

```sql
-- Example for public molecules
CREATE MATERIALIZED VIEW IF NOT EXISTS public_molecules_summary AS
SELECT 
    m.id, 
    m.name, 
    m.molecular_formula, 
    m.smiles, 
    m.cid, 
    COUNT(mp.id) AS property_count
FROM 
    molecule m
LEFT JOIN 
    molecular_property mp ON m.id = mp.molecule_id
WHERE 
    m.is_public = true
GROUP BY 
    m.id;

-- Create index on the materialized view
CREATE UNIQUE INDEX IF NOT EXISTS idx_public_molecules_summary_id 
ON public_molecules_summary(id);
```

## 3. Implementation Migration Script

Create a complete migration script that can be applied to update the RLS policies:

```sql
-- Migration: Optimized RLS Policies with Security Definer Functions
-- Purpose: Improve performance of RLS policies with security definer functions

BEGIN;

-- Step 1: Create Security Definer Functions
CREATE OR REPLACE FUNCTION public.is_project_member(project_id uuid)
RETURNS boolean AS $$
BEGIN
  RETURN EXISTS (
    SELECT 1 FROM user_profile
    WHERE user_profile.project_id = is_project_member.project_id
      AND user_profile.user_id = auth.uid()
  );
END;
$$ LANGUAGE plpgsql SECURITY DEFINER;

CREATE OR REPLACE FUNCTION public.user_projects()
RETURNS SETOF uuid AS $$
BEGIN
  RETURN QUERY
  SELECT project_id FROM user_profile
  WHERE user_profile.user_id = auth.uid();
END;
$$ LANGUAGE plpgsql SECURITY DEFINER;

CREATE OR REPLACE FUNCTION public.is_team_member(team_id uuid)
RETURNS boolean AS $$
BEGIN
  RETURN EXISTS (
    SELECT 1 FROM user_profile
    WHERE user_profile.team_id = is_team_member.team_id
      AND user_profile.user_id = auth.uid()
  );
END;
$$ LANGUAGE plpgsql SECURITY DEFINER;

CREATE OR REPLACE FUNCTION public.is_project_owner(project_id uuid)
RETURNS boolean AS $$
BEGIN
  RETURN EXISTS (
    SELECT 1 FROM user_profile
    WHERE user_profile.project_id = is_project_owner.project_id
      AND user_profile.user_id = auth.uid()
      AND user_profile.role = 'owner'
  );
END;
$$ LANGUAGE plpgsql SECURITY DEFINER;

CREATE OR REPLACE FUNCTION public.molecule_in_user_project(molecule_id uuid)
RETURNS boolean AS $$
BEGIN
  RETURN EXISTS (
    SELECT 1 FROM molecule
    JOIN user_profile ON user_profile.project_id = molecule.project_id
    WHERE molecule.id = molecule_in_user_project.molecule_id
      AND user_profile.user_id = auth.uid()
  );
END;
$$ LANGUAGE plpgsql SECURITY DEFINER;

CREATE OR REPLACE FUNCTION public.mixture_in_user_project(mixture_id uuid)
RETURNS boolean AS $$
BEGIN
  RETURN EXISTS (
    SELECT 1 FROM mixture
    JOIN user_profile ON user_profile.project_id = mixture.project_id
    WHERE mixture.id = mixture_in_user_project.mixture_id
      AND user_profile.user_id = auth.uid()
  );
END;
$$ LANGUAGE plpgsql SECURITY DEFINER;

-- Step 2: Create Performance Indexes
CREATE INDEX IF NOT EXISTS idx_user_profile_user_id ON user_profile(user_id);
CREATE INDEX IF NOT EXISTS idx_user_profile_project_id ON user_profile(project_id);
CREATE INDEX IF NOT EXISTS idx_user_profile_team_id ON user_profile(team_id);
CREATE INDEX IF NOT EXISTS idx_user_profile_role ON user_profile(role);

CREATE INDEX IF NOT EXISTS idx_molecule_project_id ON molecule(project_id);
CREATE INDEX IF NOT EXISTS idx_mixture_project_id ON mixture(project_id);
CREATE INDEX IF NOT EXISTS idx_mixture_component_mixture_id ON mixture_component(mixture_id);
CREATE INDEX IF NOT EXISTS idx_molecular_property_molecule_id ON molecular_property(molecule_id);
CREATE INDEX IF NOT EXISTS idx_experiment_project_id ON experiment(project_id);
CREATE INDEX IF NOT EXISTS idx_prediction_molecule_id ON prediction(molecule_id);
CREATE INDEX IF NOT EXISTS idx_experiment_property_experiment_id ON experiment_property(experiment_id);

-- Step 3: Replace RLS Policies with Optimized Versions

-- Molecule Table Policies
DROP POLICY IF EXISTS "Select molecules for project members" ON molecule;
CREATE POLICY "Select molecules for project members"
  ON molecule
  FOR SELECT
  USING (is_project_member(project_id));

DROP POLICY IF EXISTS "Insert molecules for project members" ON molecule;
CREATE POLICY "Insert molecules for project members"
  ON molecule
  FOR INSERT
  WITH CHECK (is_project_member(project_id));

DROP POLICY IF EXISTS "Update molecules for project members" ON molecule;
CREATE POLICY "Update molecules for project members"
  ON molecule
  FOR UPDATE
  USING (is_project_member(project_id));

DROP POLICY IF EXISTS "Delete molecules for project members" ON molecule;
CREATE POLICY "Delete molecules for project members"
  ON molecule
  FOR DELETE
  USING (is_project_member(project_id));

-- Mixture Table Policies
DROP POLICY IF EXISTS "Select mixtures for project members" ON mixture;
CREATE POLICY "Select mixtures for project members"
  ON mixture
  FOR SELECT
  USING (is_project_member(project_id));

DROP POLICY IF EXISTS "Insert mixtures for project members" ON mixture;
CREATE POLICY "Insert mixtures for project members"
  ON mixture
  FOR INSERT
  WITH CHECK (is_project_member(project_id));

DROP POLICY IF EXISTS "Update mixtures for project members" ON mixture;
CREATE POLICY "Update mixtures for project members"
  ON mixture
  FOR UPDATE
  USING (is_project_member(project_id));

DROP POLICY IF EXISTS "Delete mixtures for project members" ON mixture;
CREATE POLICY "Delete mixtures for project members"
  ON mixture
  FOR DELETE
  USING (is_project_member(project_id));

-- Mixture Component Table Policies
DROP POLICY IF EXISTS "Select mixture_components for project members" ON mixture_component;
CREATE POLICY "Select mixture_components for project members"
  ON mixture_component
  FOR SELECT
  USING (mixture_in_user_project(mixture_id));

DROP POLICY IF EXISTS "Insert mixture_components for project members" ON mixture_component;
CREATE POLICY "Insert mixture_components for project members"
  ON mixture_component
  FOR INSERT
  WITH CHECK (mixture_in_user_project(mixture_id));

DROP POLICY IF EXISTS "Update mixture_components for project members" ON mixture_component;
CREATE POLICY "Update mixture_components for project members"
  ON mixture_component
  FOR UPDATE
  USING (mixture_in_user_project(mixture_id));

DROP POLICY IF EXISTS "Delete mixture_components for project members" ON mixture_component;
CREATE POLICY "Delete mixture_components for project members"
  ON mixture_component
  FOR DELETE
  USING (mixture_in_user_project(mixture_id));

-- Add other table policies here...

-- Step 4: Add Service Role Policies
CREATE OR REPLACE FUNCTION create_service_role_policies()
RETURNS void AS $$
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
$$ LANGUAGE plpgsql;

SELECT create_service_role_policies();

-- Log completion
DO $$
BEGIN
    RAISE NOTICE 'RLS optimization migration completed successfully';
END $$;

COMMIT;
```

## 4. Testing Strategy

Creating a testing script is essential to verify that the RLS policies work correctly after optimization:

```sql
-- Begin transaction - allows rollback if tests fail
BEGIN;

-- Create test user identities
INSERT INTO auth.users (id, email) VALUES 
('11111111-1111-1111-1111-111111111111', 'test_user1@example.com'),
('22222222-2222-2222-2222-222222222222', 'test_user2@example.com');

-- Create test projects
INSERT INTO project (id, name) VALUES 
('aaaaaaaa-aaaa-aaaa-aaaa-aaaaaaaaaaaa', 'Test Project 1'),
('bbbbbbbb-bbbb-bbbb-bbbb-bbbbbbbbbbbb', 'Test Project 2');

-- Create test user profiles
INSERT INTO user_profile (user_id, project_id, role) VALUES 
('11111111-1111-1111-1111-111111111111', 'aaaaaaaa-aaaa-aaaa-aaaa-aaaaaaaaaaaa', 'owner'),
('22222222-2222-2222-2222-222222222222', 'bbbbbbbb-bbbb-bbbb-bbbb-bbbbbbbbbbbb', 'owner');

-- Create test molecules
INSERT INTO molecule (id, project_id, name, smiles) VALUES 
('cccccccc-cccc-cccc-cccc-cccccccccccc', 'aaaaaaaa-aaaa-aaaa-aaaa-aaaaaaaaaaaa', 'Test Molecule 1', 'C'),
('dddddddd-dddd-dddd-dddd-dddddddddddd', 'bbbbbbbb-bbbb-bbbb-bbbb-bbbbbbbbbbbb', 'Test Molecule 2', 'CC');

-- Test RLS as user 1
SET LOCAL ROLE authenticated;
SET LOCAL auth.uid = '11111111-1111-1111-1111-111111111111';

-- Test project access
SELECT * FROM project;
-- Should see only Test Project 1

-- Test molecule access
SELECT * FROM molecule;
-- Should see only Test Molecule 1

-- Switch to user 2
SET LOCAL auth.uid = '22222222-2222-2222-2222-222222222222';

-- Test project access
SELECT * FROM project;
-- Should see only Test Project 2

-- Test molecule access
SELECT * FROM molecule;
-- Should see only Test Molecule 2

-- Switch to service role
SET ROLE service_role;

-- Test all access
SELECT * FROM project;
-- Should see all projects

SELECT * FROM molecule;
-- Should see all molecules

-- Clean up test data
ROLLBACK;
```

## 5. Performance Monitoring

To ensure the RLS optimizations actually improve performance, implement a simple benchmark script:

```sql
-- Function to benchmark RLS performance
CREATE OR REPLACE FUNCTION benchmark_rls_performance()
RETURNS TABLE (query_type text, execution_time_ms numeric) AS $$
DECLARE
    start_time timestamptz;
    end_time timestamptz;
    user_id uuid := auth.uid();
    project_count int;
    molecule_count int;
BEGIN
    -- Count number of projects a user can access
    start_time := clock_timestamp();
    SELECT COUNT(*) INTO project_count FROM project;
    end_time := clock_timestamp();
    query_type := 'project_access';
    execution_time_ms := extract(epoch from (end_time - start_time)) * 1000;
    RETURN NEXT;
    
    -- Count number of molecules a user can access
    start_time := clock_timestamp();
    SELECT COUNT(*) INTO molecule_count FROM molecule;
    end_time := clock_timestamp();
    query_type := 'molecule_access';
    execution_time_ms := extract(epoch from (end_time - start_time)) * 1000;
    RETURN NEXT;
    
    -- Join query performance
    start_time := clock_timestamp();
    PERFORM COUNT(*) FROM molecule m 
    JOIN molecular_property mp ON m.id = mp.molecule_id;
    end_time := clock_timestamp();
    query_type := 'molecule_property_join';
    execution_time_ms := extract(epoch from (end_time - start_time)) * 1000;
    RETURN NEXT;
END;
$$ LANGUAGE plpgsql;

-- Run benchmark as a specific user
SET LOCAL ROLE authenticated;
SET LOCAL auth.uid = '11111111-1111-1111-1111-111111111111';
SELECT * FROM benchmark_rls_performance();
```

## 6. Implementation Timeline

1. **Day 1**: Create security definer functions and performance indexes
2. **Day 2**: Update RLS policies for core tables
3. **Day 3**: Create service role policies and test the implementation
4. **Day 4**: Monitor performance and make adjustments

## 7. Conclusion

This practical approach to RLS optimization focuses on real-world benefits:

1. **Reduced query complexity**: Security definer functions simplify policy expressions
2. **Better indexing**: Strategic indexes improve policy evaluation performance
3. **Service role consistency**: Uniform access for service roles across all tables
4. **Maintainability**: Clear, modular approach to access control
5. **Testability**: Comprehensive testing strategy ensures security is maintained

By implementing these optimizations, the CryoProtect database should see significant performance improvements while maintaining proper security.