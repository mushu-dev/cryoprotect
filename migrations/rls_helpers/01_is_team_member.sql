CREATE OR REPLACE FUNCTION auth.is_team_member(team_id UUID) 
RETURNS BOOLEAN AS 
$$
BEGIN
  RETURN EXISTS (
    SELECT 1 FROM user_teams 
    WHERE team_id = $1 
    AND user_id = auth.uid()
  );
END;
$$ 
LANGUAGE plpgsql SECURITY DEFINER;