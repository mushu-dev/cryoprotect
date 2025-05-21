CREATE OR REPLACE FUNCTION auth.is_team_admin(team_id UUID) 
RETURNS BOOLEAN AS 
$$
BEGIN
  RETURN EXISTS (
    SELECT 1 FROM team_members 
    WHERE team_id = $1 
    AND user_id = auth.uid()
    AND role = 'admin'
  );
END;
$$ 
LANGUAGE plpgsql SECURITY DEFINER;