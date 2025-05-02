-- Create a function to execute SQL queries for verification purposes
-- This function allows executing arbitrary SQL queries from Python
-- It is used by the verification script to run complex queries

-- Drop the function if it already exists
DROP FUNCTION IF EXISTS exec_sql(text);

-- Create the function
CREATE OR REPLACE FUNCTION exec_sql(query text)
RETURNS SETOF json AS $$
BEGIN
  RETURN QUERY EXECUTE query;
END;
$$ LANGUAGE plpgsql SECURITY DEFINER;

-- Grant execute permission to authenticated users
GRANT EXECUTE ON FUNCTION exec_sql TO authenticated;

-- Grant execute permission to anon users (if needed for public access)
-- GRANT EXECUTE ON FUNCTION exec_sql TO anon;

-- Output confirmation
DO $$
BEGIN
  RAISE NOTICE 'exec_sql function created successfully';
END $$;