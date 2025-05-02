# MCP DDL Execution Guide

This document provides guidance on executing DDL (Data Definition Language) statements via the MCP tool suite in the CryoProtect system.

## Problem Statement

The MCP `execute_sql` tool only supports SQL queries that return tuples (e.g., SELECT statements). It cannot directly execute DDL statements like CREATE TABLE, ALTER TABLE, or CREATE INDEX, resulting in the error:

```
{'code': '42601', 'message': 'ALTER TABLE query does not return tuples'}
```

This limitation can block critical database schema changes required for system enhancements.

## Solution Patterns

### 1. SQL Function Wrapper Pattern

Create a PostgreSQL function that wraps DDL commands in a tuple-returning function:

```sql
-- Create the helper function (only needs to be done once)
CREATE OR REPLACE FUNCTION execute_ddl(ddl_command text) 
RETURNS TABLE(result text) AS $$
BEGIN
    EXECUTE ddl_command;
    RETURN QUERY SELECT 'Success: ' || ddl_command AS result;
EXCEPTION WHEN OTHERS THEN
    RETURN QUERY SELECT 'Error: ' || SQLERRM AS result;
END;
$$ LANGUAGE plpgsql SECURITY DEFINER;

-- Use the function to execute DDL
SELECT execute_ddl('ALTER TABLE molecules ADD COLUMN chembl_id VARCHAR;') as result;
```

### 2. RETURNING Clause for DML Statements

For data manipulation operations, add a RETURNING clause to make them return tuples:

```sql
-- Instead of:
UPDATE molecules SET some_field = 'value' WHERE condition;

-- Use:
UPDATE molecules SET some_field = 'value' WHERE condition RETURNING id, some_field;
```

### 3. Split Multi-Statement Migrations

Split migration files into individual statements that can be executed separately:

```python
# Split the migration SQL into separate statements
sql_statements = migration_sql.split(';')
sql_statements = [stmt.strip() for stmt in sql_statements if stmt.strip()]

# Execute each SQL statement
for i, sql in enumerate(sql_statements):
    if sql.strip().upper().startswith(('SELECT', 'UPDATE', 'INSERT', 'DELETE')):
        # Execute the SQL statement using MCP
        result = execute_sql(sql, project_id)
    else:
        # Skip non-tuple-returning statements or handle differently
        logger.warning(f"Skipping non-tuple-returning statement")
```

### 4. Direct Database Connection Fallback

For critical DDL operations, consider bypassing MCP and using a direct database connection:

```python
def execute_direct_ddl(sql):
    """Execute DDL using direct psycopg2 connection."""
    import psycopg2
    from config import Config
    
    config = Config()
    
    try:
        conn = psycopg2.connect(
            host=config.SUPABASE_DB_HOST,
            port=config.SUPABASE_DB_PORT,
            database=config.SUPABASE_DB_NAME,
            user=config.SUPABASE_DB_USER,
            password=config.SUPABASE_DB_PASSWORD
        )
        conn.autocommit = True
        
        cursor = conn.cursor()
        cursor.execute(sql)
        cursor.close()
        conn.close()
        
        return True
    except Exception as e:
        logger.error(f"Error executing DDL: {str(e)}")
        return False
```

## Example: ChEMBL Remediation Application

The ChEMBL remediation process required schema changes that were blocked by this limitation. We developed a workaround using these patterns:

1. Created a helper function in the database to wrap DDL in SELECT statements
2. Modified the migration logic to split statements and process them appropriately
3. Added verification steps to ensure schema changes were correctly applied

For the complete implementation, see the `apply_migration` function in `chembl_remediation_main_fixed.py`.

## Best Practices

1. **Always Check First**: Use `information_schema` to check if objects exist before trying to create them
2. **Transaction Safety**: Ensure DDL changes are wrapped in appropriate transactions
3. **Verification**: Always verify that schema changes were applied correctly
4. **Documentation**: Document any workarounds used for DDL execution
5. **Error Handling**: Implement robust error handling for DDL operations

## Long-Term Recommendations

For a more permanent solution:

1. Enhance the MCP tool suite to handle DDL statements directly
2. Create a dedicated MCP function specifically for executing DDL (e.g., `execute_ddl`)
3. Maintain a library of SQL functions that can wrap DDL in tuple-returning operations
4. Document this limitation clearly in the MCP documentation

---

*Last Updated: April 27, 2025*