CREATE OR REPLACE FUNCTION auth.can_access_molecule(molecule_id UUID) 
RETURNS BOOLEAN AS 
$$
BEGIN
  RETURN EXISTS (
    SELECT 1
    FROM experiments e
    JOIN projects p ON e.project_id = p.id
    JOIN teams t ON p.team_id = t.id
    JOIN user_teams ut ON t.id = ut.team_id
    WHERE e.molecule_id = molecule_id
    AND ut.user_id = auth.uid()
  );
END;
$$ 
LANGUAGE plpgsql SECURITY DEFINER;