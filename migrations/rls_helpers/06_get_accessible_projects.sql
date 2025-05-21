CREATE OR REPLACE FUNCTION auth.get_accessible_projects() 
RETURNS SETOF UUID AS 
$$
BEGIN
  RETURN QUERY
  SELECT p.id
  FROM projects p
  JOIN teams t ON p.team_id = t.id
  JOIN user_teams ut ON t.id = ut.team_id
  WHERE ut.user_id = auth.uid();
END;
$$ 
LANGUAGE plpgsql SECURITY DEFINER;