# Database Verification Plan

This document outlines a comprehensive verification plan for the CryoProtect database implementation, with a focus on verifying the RLS optimization, connection pooling, and overall database architecture. This plan ensures we have a solid foundation before moving to the next phase of development.

## 1. RLS Policy Verification

### 1.1 Security Definer Functions

| Test Case | Description | Expected Result |
|-----------|-------------|-----------------|
| TC-RLS-001 | Verify `is_project_member` function with valid project | Returns true for valid project membership |
| TC-RLS-002 | Verify `is_project_member` function with invalid project | Returns false for invalid project membership |
| TC-RLS-003 | Verify `user_projects` function | Returns all projects for the current user |
| TC-RLS-004 | Verify `is_team_member` function | Returns correct team membership status |
| TC-RLS-005 | Verify `is_project_owner` function | Returns true only for project owners |
| TC-RLS-006 | Verify `molecule_in_user_project` function | Returns true for molecules in user's projects |
| TC-RLS-007 | Verify `mixture_in_user_project` function | Returns true for mixtures in user's projects |

```sql
-- Example test script for TC-RLS-001
DO $$
DECLARE
    result boolean;
BEGIN
    -- Set up test user and project
    SET LOCAL ROLE authenticated;
    SET LOCAL auth.uid = '00000000-0000-0000-0000-000000000001';
    
    -- Test with a project the user is a member of
    SELECT is_project_member('11111111-1111-1111-1111-111111111111') INTO result;
    
    -- Verify result
    ASSERT result = true, 'Expected is_project_member to return true for a valid project membership';
    
    -- Test with a project the user is not a member of
    SELECT is_project_member('22222222-2222-2222-2222-222222222222') INTO result;
    
    -- Verify result
    ASSERT result = false, 'Expected is_project_member to return false for an invalid project membership';
END $$;
```

### 1.2 Table Access Policies

| Test Case | Description | Expected Result |
|-----------|-------------|-----------------|
| TC-RLS-010 | Test SELECT access to `molecule` table as project member | Can view molecules in own projects |
| TC-RLS-011 | Test SELECT access to `molecule` table as non-member | Cannot view molecules in other projects |
| TC-RLS-012 | Test INSERT into `molecule` table as project member | Can insert molecules into own projects |
| TC-RLS-013 | Test INSERT into `molecule` table as non-member | Cannot insert molecules into other projects |
| TC-RLS-014 | Test UPDATE on `molecule` table as project member | Can update molecules in own projects |
| TC-RLS-015 | Test UPDATE on `molecule` table as non-member | Cannot update molecules in other projects |
| TC-RLS-016 | Test DELETE from `molecule` table as project member | Can delete molecules in own projects |
| TC-RLS-017 | Test DELETE from `molecule` table as non-member | Cannot delete molecules in other projects |

```sql
-- Example test script for TC-RLS-010 and TC-RLS-011
DO $$
DECLARE
    molecule_count integer;
BEGIN
    -- Set up test user 1
    SET LOCAL ROLE authenticated;
    SET LOCAL auth.uid = '00000000-0000-0000-0000-000000000001';
    
    -- Count molecules for user 1 (should see their project's molecules)
    SELECT COUNT(*) INTO molecule_count FROM molecule;
    
    -- Verify result
    ASSERT molecule_count > 0, 'Expected user 1 to see molecules in their project';
    
    -- Set up test user 2
    SET LOCAL ROLE authenticated;
    SET LOCAL auth.uid = '00000000-0000-0000-0000-000000000002';
    
    -- Count molecules for user 2 (should only see their project's molecules)
    SELECT COUNT(*) INTO molecule_count FROM molecule 
    WHERE project_id = '11111111-1111-1111-1111-111111111111';
    
    -- Verify result
    ASSERT molecule_count = 0, 'Expected user 2 to NOT see molecules from user 1''s project';
END $$;
```

### 1.3 Relationship Policies

| Test Case | Description | Expected Result |
|-----------|-------------|-----------------|
| TC-RLS-020 | Test access to `molecular_property` for molecule in user's project | Can access properties |
| TC-RLS-021 | Test access to `molecular_property` for molecule not in user's project | Cannot access properties |
| TC-RLS-022 | Test access to `mixture_component` for mixture in user's project | Can access components |
| TC-RLS-023 | Test access to `mixture_component` for mixture not in user's project | Cannot access components |
| TC-RLS-024 | Test access to `experiment_property` for experiment in user's project | Can access properties |
| TC-RLS-025 | Test access to `experiment_property` for experiment not in user's project | Cannot access properties |

### 1.4 Service Role Access

| Test Case | Description | Expected Result |
|-----------|-------------|-----------------|
| TC-RLS-030 | Test SELECT access as service role to all tables | Can view all data |
| TC-RLS-031 | Test INSERT access as service role to all tables | Can insert data anywhere |
| TC-RLS-032 | Test UPDATE access as service role to all tables | Can update all data |
| TC-RLS-033 | Test DELETE access as service role to all tables | Can delete all data |

```sql
-- Example test script for TC-RLS-030
DO $$
DECLARE
    table_name text;
    row_count integer;
BEGIN
    -- Set service role
    SET ROLE service_role;
    
    -- Loop through main tables
    FOR table_name IN 
        SELECT tablename FROM pg_tables 
        WHERE schemaname = 'public' 
        AND tablename IN ('molecule', 'mixture', 'experiment', 'molecular_property', 'prediction')
    LOOP
        -- Get row count for table
        EXECUTE format('SELECT COUNT(*) FROM %I', table_name) INTO row_count;
        
        -- Log result
        RAISE NOTICE 'Service role can access % rows in table %', row_count, table_name;
    END LOOP;
END $$;
```

## 2. Connection Pool Verification

### 2.1 Basic Functionality

| Test Case | Description | Expected Result |
|-----------|-------------|-----------------|
| TC-CP-001 | Initialize connection pool | Pool created successfully |
| TC-CP-002 | Get and return connection | Connection successfully obtained and returned |
| TC-CP-003 | Execute simple query | Query executes successfully |
| TC-CP-004 | Execute transaction | Transaction commits successfully |
| TC-CP-005 | Handle connection errors | Errors handled gracefully |

### 2.2 Performance Testing

| Test Case | Description | Expected Result |
|-----------|-------------|-----------------|
| TC-CP-010 | Measure query time with direct connection | Baseline performance established |
| TC-CP-011 | Measure query time with connection pool | Improved or similar performance |
| TC-CP-012 | Measure concurrent query performance | Handles concurrent requests efficiently |
| TC-CP-013 | Measure connection acquisition time | Fast connection acquisition |
| TC-CP-014 | Measure pool scaling under load | Proper scaling behavior |

```python
# Example test script for TC-CP-011 and TC-CP-012
import time
import concurrent.futures
from db_pool import execute_query

def test_connection_pool_performance():
    # Test query
    query = "SELECT COUNT(*) FROM molecule"
    
    # Single query performance
    start_time = time.time()
    result = execute_query(query)
    single_query_time = time.time() - start_time
    print(f"Single query time: {single_query_time:.6f}s")
    
    # Concurrent query performance
    num_concurrent = 10
    start_time = time.time()
    
    with concurrent.futures.ThreadPoolExecutor(max_workers=num_concurrent) as executor:
        futures = [executor.submit(execute_query, query) for _ in range(num_concurrent)]
        results = [future.result() for future in concurrent.futures.as_completed(futures)]
    
    concurrent_total_time = time.time() - start_time
    print(f"Concurrent queries ({num_concurrent}): {concurrent_total_time:.6f}s")
    print(f"Average time per query: {concurrent_total_time/num_concurrent:.6f}s")
    
    # Compare with theoretical sequential time
    sequential_time = single_query_time * num_concurrent
    print(f"Sequential equivalent: {sequential_time:.6f}s")
    print(f"Speedup factor: {sequential_time/concurrent_total_time:.2f}x")

if __name__ == "__main__":
    test_connection_pool_performance()
```

### 2.3 Resilience Testing

| Test Case | Description | Expected Result |
|-----------|-------------|-----------------|
| TC-CP-020 | Test connection retry logic | Retries failed connections |
| TC-CP-021 | Test connection timeout | Handles timeouts gracefully |
| TC-CP-022 | Test connection validation | Validates and replaces invalid connections |
| TC-CP-023 | Test connection lifecycle | Manages connection lifecycle properly |
| TC-CP-024 | Test circuit breaker | Activates and resets circuit breaker |

## 3. Database Schema Verification

### 3.1 Table Structure

| Test Case | Description | Expected Result |
|-----------|-------------|-----------------|
| TC-DB-001 | Verify molecule table schema | Matches expected schema |
| TC-DB-002 | Verify mixture table schema | Matches expected schema |
| TC-DB-003 | Verify experiment table schema | Matches expected schema |
| TC-DB-004 | Verify molecular_property table schema | Matches expected schema |
| TC-DB-005 | Verify user_profile table schema | Matches expected schema |

### 3.2 Foreign Key Relationships

| Test Case | Description | Expected Result |
|-----------|-------------|-----------------|
| TC-DB-010 | Verify molecule → project relationship | Foreign key enforced |
| TC-DB-011 | Verify mixture → project relationship | Foreign key enforced |
| TC-DB-012 | Verify mixture_component → mixture relationship | Foreign key enforced |
| TC-DB-013 | Verify molecular_property → molecule relationship | Foreign key enforced |
| TC-DB-014 | Verify experiment → project relationship | Foreign key enforced |
| TC-DB-015 | Verify experiment_property → experiment relationship | Foreign key enforced |

```sql
-- Example test script for TC-DB-010
DO $$
DECLARE
    result boolean;
BEGIN
    -- Check if foreign key exists
    SELECT EXISTS (
        SELECT 1
        FROM information_schema.table_constraints tc
        JOIN information_schema.constraint_column_usage ccu ON ccu.constraint_name = tc.constraint_name
        WHERE tc.constraint_type = 'FOREIGN KEY'
        AND tc.table_name = 'molecule'
        AND ccu.table_name = 'project'
        AND ccu.column_name = 'id'
    ) INTO result;
    
    -- Verify result
    ASSERT result = true, 'Foreign key from molecule to project not found';
    
    -- Test foreign key enforcement
    BEGIN
        -- Try to insert a molecule with non-existent project_id
        INSERT INTO molecule (id, project_id, name, smiles)
        VALUES (gen_random_uuid(), 'ffffffff-ffff-ffff-ffff-ffffffffffff', 'Test Molecule', 'C');
        
        -- Should not reach here
        ASSERT false, 'Foreign key constraint not enforced';
    EXCEPTION
        WHEN foreign_key_violation THEN
            -- Expected behavior
            NULL;
    END;
END $$;
```

### 3.3 Index Coverage

| Test Case | Description | Expected Result |
|-----------|-------------|-----------------|
| TC-DB-020 | Verify indexes for RLS policy columns | Required indexes exist |
| TC-DB-021 | Verify indexes for common query patterns | Required indexes exist |
| TC-DB-022 | Verify index usage in queries | Indexes are used in execution plans |

```sql
-- Example test script for TC-DB-020
DO $$
DECLARE
    index_exists boolean;
BEGIN
    -- Check for user_profile(user_id) index
    SELECT EXISTS (
        SELECT 1
        FROM pg_indexes
        WHERE schemaname = 'public'
        AND tablename = 'user_profile'
        AND indexname LIKE '%user_id%'
    ) INTO index_exists;
    
    -- Verify result
    ASSERT index_exists = true, 'Index for user_profile(user_id) not found';
    
    -- Check for molecule(project_id) index
    SELECT EXISTS (
        SELECT 1
        FROM pg_indexes
        WHERE schemaname = 'public'
        AND tablename = 'molecule'
        AND indexname LIKE '%project_id%'
    ) INTO index_exists;
    
    -- Verify result
    ASSERT index_exists = true, 'Index for molecule(project_id) not found';
END $$;
```

## 4. Migration Framework Verification

### 4.1 Migration Application

| Test Case | Description | Expected Result |
|-----------|-------------|-----------------|
| TC-MG-001 | Apply simple migration | Migration applies successfully |
| TC-MG-002 | Apply migration with rollback | Rollback functions correctly |
| TC-MG-003 | Apply migration with data transformation | Data transformed correctly |
| TC-MG-004 | Apply migration with constraints | Constraints enforced correctly |

### 4.2 Migration Tracking

| Test Case | Description | Expected Result |
|-----------|-------------|-----------------|
| TC-MG-010 | Track migration versions | Versions tracked correctly |
| TC-MG-011 | Detect and prevent duplicate migrations | Duplicates prevented |
| TC-MG-012 | Handle migration dependencies | Dependencies resolved correctly |
| TC-MG-013 | Generate migration report | Report accurately reflects status |

## 5. Implementation Plan

1. Create a test database environment
2. Implement test scripts for RLS policy verification
3. Implement test scripts for connection pool verification
4. Implement test scripts for database schema verification
5. Implement test scripts for migration framework verification
6. Run tests and document results
7. Fix any issues found during testing
8. Create performance baseline measurements
9. Document the test results and verification process

## 6. Verification Metrics

- **Coverage**: Percentage of test cases executed and passed
- **Performance Improvement**: Query execution time compared to baseline
- **Concurrent Load Handling**: Number of concurrent connections handled
- **Error Rates**: Percentage of operations that result in errors
- **Recovery Time**: Time to recover from simulated failures

This verification plan ensures that we have a solid database foundation in place before proceeding to the next phase of the project. By thoroughly testing the RLS optimization, connection pooling, and overall database architecture, we can be confident that our system will perform well and maintain data security.

## 7. Tools and Utilities

### RLS Test Helper

```python
# rls_test_helper.py
import psycopg2
from psycopg2.extras import RealDictCursor
import uuid
import json

class RLSTestHelper:
    def __init__(self, conn_params):
        self.conn_params = conn_params
        
    def create_test_data(self):
        """Create test data for RLS testing"""
        conn = None
        try:
            # Connect to database
            conn = psycopg2.connect(**self.conn_params)
            conn.autocommit = False
            cursor = conn.cursor()
            
            # Create test users
            user1_id = uuid.uuid4()
            user2_id = uuid.uuid4()
            
            # Create test projects
            project1_id = uuid.uuid4()
            project2_id = uuid.uuid4()
            
            # Create test user profiles
            cursor.execute(f"""
                INSERT INTO user_profile (user_id, project_id, role)
                VALUES 
                  ('{user1_id}', '{project1_id}', 'owner'),
                  ('{user2_id}', '{project2_id}', 'owner')
            """)
            
            # Create test molecules
            cursor.execute(f"""
                INSERT INTO molecule (id, project_id, name, smiles)
                VALUES 
                  ('{uuid.uuid4()}', '{project1_id}', 'Test Molecule 1', 'C'),
                  ('{uuid.uuid4()}', '{project2_id}', 'Test Molecule 2', 'CC')
            """)
            
            # Commit changes
            conn.commit()
            
            # Return test data IDs
            return {
                'users': [str(user1_id), str(user2_id)],
                'projects': [str(project1_id), str(project2_id)]
            }
        except Exception as e:
            # Rollback on error
            if conn:
                conn.rollback()
            raise e
        finally:
            # Close connection
            if conn:
                conn.close()
    
    def test_user_access(self, user_id, table_name, condition=None):
        """Test user access to a table"""
        conn = None
        try:
            # Connect to database
            conn = psycopg2.connect(**self.conn_params)
            cursor = conn.cursor(cursor_factory=RealDictCursor)
            
            # Set auth context
            cursor.execute(f"SET LOCAL ROLE authenticated")
            cursor.execute(f"SET LOCAL auth.uid = '{user_id}'")
            
            # Execute query
            query = f"SELECT * FROM {table_name}"
            if condition:
                query += f" WHERE {condition}"
            cursor.execute(query)
            
            # Get results
            results = cursor.fetchall()
            
            # Return result count and data
            return {
                'count': len(results),
                'data': results
            }
        except Exception as e:
            # Handle error
            return {
                'error': str(e)
            }
        finally:
            # Close connection
            if conn:
                conn.close()
    
    def test_service_role_access(self, table_name):
        """Test service role access to a table"""
        conn = None
        try:
            # Connect to database
            conn = psycopg2.connect(**self.conn_params)
            cursor = conn.cursor(cursor_factory=RealDictCursor)
            
            # Set auth context
            cursor.execute(f"SET ROLE service_role")
            
            # Execute query
            cursor.execute(f"SELECT * FROM {table_name}")
            
            # Get results
            results = cursor.fetchall()
            
            # Return result count and data
            return {
                'count': len(results),
                'data': results
            }
        except Exception as e:
            # Handle error
            return {
                'error': str(e)
            }
        finally:
            # Close connection
            if conn:
                conn.close()
    
    def cleanup_test_data(self, test_data):
        """Clean up test data"""
        conn = None
        try:
            # Connect to database
            conn = psycopg2.connect(**self.conn_params)
            conn.autocommit = False
            cursor = conn.cursor()
            
            # Clean up data in reverse order of dependencies
            for project_id in test_data['projects']:
                # Delete molecules
                cursor.execute(f"DELETE FROM molecule WHERE project_id = '{project_id}'")
                
                # Delete user profiles
                cursor.execute(f"DELETE FROM user_profile WHERE project_id = '{project_id}'")
                
                # Delete projects
                cursor.execute(f"DELETE FROM project WHERE id = '{project_id}'")
            
            # Commit changes
            conn.commit()
            
            return True
        except Exception as e:
            # Rollback on error
            if conn:
                conn.rollback()
            raise e
        finally:
            # Close connection
            if conn:
                conn.close()
```

### Connection Pool Test Helper

```python
# connection_pool_test.py
import time
import concurrent.futures
import statistics
from db_pool import ConnectionPool, execute_query

class ConnectionPoolTester:
    def __init__(self, config):
        self.config = config
    
    def test_basic_functionality(self):
        """Test basic connection pool functionality"""
        results = {}
        
        # Initialize pool
        pool = ConnectionPool(self.config)
        results['initialization'] = "Success"
        
        # Get connection
        try:
            conn, conn_id = pool.get_connection()
            results['get_connection'] = "Success"
            
            # Return connection
            pool.return_connection(conn, conn_id)
            results['return_connection'] = "Success"
        except Exception as e:
            results['error'] = str(e)
        
        # Execute simple query
        try:
            result = execute_query("SELECT 1 as test")
            if result and result[0]['test'] == 1:
                results['execute_query'] = "Success"
            else:
                results['execute_query'] = "Failed"
        except Exception as e:
            results['execute_query_error'] = str(e)
        
        return results
    
    def test_performance(self, query, iterations=5, concurrent=10):
        """Test connection pool performance"""
        results = {
            'single_query_times': [],
            'concurrent_query_times': []
        }
        
        # Single query performance
        for i in range(iterations):
            start_time = time.time()
            execute_query(query)
            query_time = time.time() - start_time
            results['single_query_times'].append(query_time)
        
        # Calculate single query statistics
        results['single_query_avg'] = statistics.mean(results['single_query_times'])
        results['single_query_median'] = statistics.median(results['single_query_times'])
        results['single_query_max'] = max(results['single_query_times'])
        results['single_query_min'] = min(results['single_query_times'])
        
        # Concurrent query performance
        for i in range(iterations):
            start_time = time.time()
            
            with concurrent.futures.ThreadPoolExecutor(max_workers=concurrent) as executor:
                futures = [executor.submit(execute_query, query) for _ in range(concurrent)]
                results_list = [future.result() for future in concurrent.futures.as_completed(futures)]
            
            total_time = time.time() - start_time
            results['concurrent_query_times'].append(total_time)
        
        # Calculate concurrent query statistics
        results['concurrent_query_avg'] = statistics.mean(results['concurrent_query_times'])
        results['concurrent_query_median'] = statistics.median(results['concurrent_query_times'])
        results['concurrent_query_max'] = max(results['concurrent_query_times'])
        results['concurrent_query_min'] = min(results['concurrent_query_times'])
        
        # Calculate speedup
        theoretical_sequential = results['single_query_avg'] * concurrent
        actual_concurrent = results['concurrent_query_avg']
        results['speedup_factor'] = theoretical_sequential / actual_concurrent
        
        return results
```

## 8. Verification Checklist

- [ ] RLS security definer functions correctly implemented
- [ ] RLS policies correctly applied to all tables
- [ ] Performance indexes created and verified
- [ ] Service role policies implemented and tested
- [ ] Connection pool basic functionality verified
- [ ] Connection pool performance metrics established
- [ ] Connection pool resilience verified
- [ ] Database schema fully verified
- [ ] Foreign key relationships verified
- [ ] Migration framework verified
- [ ] All test cases executed and documented
- [ ] Performance baseline established
- [ ] Documentation completed