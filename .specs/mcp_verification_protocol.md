# MCP Verification Protocol for ChEMBL Integration

## Purpose
To verify that the Model Context Protocol (MCP) server, as configured in `.roo/mcp.json`, is fully capable of supporting all database operations required for ChEMBL integration into CryoProtect v2. This includes executing SQL, modifying Row Level Security (RLS) policies, and performing schema/permission management as outlined in `CHEMBL_INTEGRATION_ROADMAP.md`.

## Prerequisites
- `.roo/mcp.json` is present and correctly configured for the Supabase MCP server.
- MCP server is running and accessible.
- Sufficient permissions are granted to the MCP service account to perform all required operations.

## Verification Steps

### 1. SQL Execution Capability
- **Action:** Execute a simple SQL query to verify connectivity and role.
- **SQL:**
  ```sql
  SELECT current_setting('role'), current_user;
  ```
- **Expected Result:** Returns the current database role and user.

### 2. RLS Policy Inspection
- **Action:** List all RLS policies affecting the `molecules` table.
- **SQL:**
  ```sql
  SELECT tablename, policyname, permissive, roles, cmd, qual, with_check 
  FROM pg_policies 
  WHERE tablename = 'molecules';
  ```
- **Expected Result:** Returns all RLS policies for the `molecules` table.

### 3. Schema Validation
- **Action:** Verify the structure of the `molecules` table.
- **SQL:**
  ```sql
  SELECT column_name, data_type, is_nullable
  FROM information_schema.columns
  WHERE table_name = 'molecules'
  ORDER BY ordinal_position;
  ```
- **Expected Result:** Returns the schema of the `molecules` table.

### 4. RLS Policy Modification
- **Action:** Test the ability to modify RLS policies.
- **SQL Examples:**
  - Disable RLS:
    ```sql
    ALTER TABLE molecules DISABLE ROW LEVEL SECURITY;
    ```
  - Create a permissive policy:
    ```sql
    CREATE POLICY "Service role has full access to molecules"
    ON molecules
    FOR ALL
    USING (auth.jwt() ->> 'role' = 'service_role')
    WITH CHECK (auth.jwt() ->> 'role' = 'service_role');
    ```
  - Re-enable RLS:
    ```sql
    ALTER TABLE molecules ENABLE ROW LEVEL SECURITY;
    ```
  - Drop temporary policy:
    ```sql
    DROP POLICY IF EXISTS "Service role has full access to molecules" ON molecules;
    ```
- **Expected Result:** Each command executes successfully and the policy state changes as expected.

### 5. Permission Management
- **Action:** Test granting privileges to roles.
- **SQL:**
  ```sql
  GRANT ALL PRIVILEGES ON TABLE molecules TO authenticated;
  GRANT ALL PRIVILEGES ON TABLE molecular_properties TO authenticated;
  GRANT USAGE, SELECT ON SEQUENCE molecules_id_seq TO authenticated;
  ```
- **Expected Result:** Privileges are granted without error.

### 6. Batch DML and Transaction Handling
- **Action:** Test batch inserts and transaction management.
- **SQL Example:**
  ```sql
  BEGIN;
  INSERT INTO molecules (name, smiles, inchi, inchikey, formula, molecular_weight, created_by, data_source, version)
  VALUES ('TestMol', 'C1=CC=CC=C1', 'InChI=1S/C6H6/c1-2-4-6-5-3-1/h1-6H', 'UHOVQNZJYSORNB-UHFFFAOYSA-N', 'C6H6', 78.11, 'test_user', 'TestSource', 1)
  ON CONFLICT (inchikey) DO NOTHING
  RETURNING id;
  COMMIT;
  ```
- **Expected Result:** Insert succeeds and transaction is committed.

## Documentation & Issue Reporting

- For each verification step, record:
  - The SQL executed
  - The result (success/failure, output)
  - Any errors or issues encountered
- If any step fails, document the error and attempt to diagnose (e.g., permission denied, unsupported operation).
- Summarize findings and recommendations for remediation if issues are found.

## References

- `.roo/mcp.json` (MCP configuration)
- `CHEMBL_INTEGRATION_ROADMAP.md` (integration requirements)
- Supabase/Postgres documentation for RLS and permissions