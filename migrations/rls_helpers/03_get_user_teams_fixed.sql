CREATE OR REPLACE FUNCTION auth.get_user_teams() 
RETURNS SETOF UUID AS 
$$
BEGIN
  RETURN QUERY
  SELECT team_id FROM team_members 
  WHERE user_id = auth.uid();
END;
$$ 
LANGUAGE plpgsql SECURITY DEFINER;