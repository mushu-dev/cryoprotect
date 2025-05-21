-- Function to check if the current user has the service role
CREATE OR REPLACE FUNCTION auth.is_service_role() 
RETURNS BOOLEAN AS 
$$
BEGIN
  RETURN (SELECT (current_setting('request.jwt.claims', true)::jsonb->'role')::text = '"service_role"');
EXCEPTION
  WHEN OTHERS THEN
    RETURN FALSE;
END;
$$
LANGUAGE plpgsql SECURITY DEFINER;