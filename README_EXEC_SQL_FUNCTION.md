# Creating the exec_sql Function for Database Verification

The database verification script requires a PostgreSQL function called `exec_sql` to execute SQL queries directly from Python. This document explains how to create this function in your Supabase database.

## Why is this function needed?

The verification script needs to run complex SQL queries that join multiple tables and perform aggregations. The standard Supabase client doesn't support these complex queries directly, so we use a PostgreSQL function that can execute arbitrary SQL.

## Security Considerations

The `exec_sql` function executes SQL queries directly, which could potentially be a security risk if misused. However, in our implementation:

1. The function is only accessible to authenticated users
2. RLS policies still apply to the results
3. We're only using it for read-only verification queries

## Creating the Function

You can create the function using the SQL Editor in the Supabase dashboard or by running the `create_exec_sql_function.sql` script.

### Option 1: Using the Supabase Dashboard

1. Log in to your Supabase dashboard
2. Go to the SQL Editor
3. Create a new query
4. Paste the following SQL:

```sql
-- Create a function to execute SQL queries
CREATE OR REPLACE FUNCTION exec_sql(query text)
RETURNS SETOF json AS $$
BEGIN
  RETURN QUERY EXECUTE query;
END;
$$ LANGUAGE plpgsql SECURITY DEFINER;

-- Grant execute permission to authenticated users
GRANT EXECUTE ON FUNCTION exec_sql TO authenticated;
```

5. Run the query

### Option 2: Using the Provided Script

Run the `create_exec_sql_function.bat` (Windows) or `create_exec_sql_function.sh` (Unix) script:

```bash
# Windows
create_exec_sql_function.bat

# Unix
./create_exec_sql_function.sh
```

## Verifying the Function

To verify that the function was created successfully, you can run:

```sql
SELECT exec_sql('SELECT 1 as test');
```

This should return:

```json
[{"test": 1}]
```

## Troubleshooting

If you encounter errors when creating or using the function:

1. **Permission denied**: Make sure you're using an account with sufficient privileges (service role or admin)
2. **Function already exists**: You can drop the existing function with `DROP FUNCTION IF EXISTS exec_sql(text);` and then recreate it
3. **RLS errors**: The function respects RLS policies, so make sure your user has access to the tables being queried