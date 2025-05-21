CREATE OR REPLACE FUNCTION auth.get_accessible_projects() 
RETURNS SETOF UUID AS 
$$
BEGIN
  RETURN QUERY
  SELECT p.id
  FROM projects p
  JOIN teams t ON p.team_id = t.id
  JOIN team_members tm ON t.id = tm.team_id
  WHERE tm.user_id = auth.uid();
END;
$$ 
LANGUAGE plpgsql SECURITY DEFINER;