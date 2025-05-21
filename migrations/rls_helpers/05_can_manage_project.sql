CREATE OR REPLACE FUNCTION auth.can_manage_project(project_id UUID) 
RETURNS BOOLEAN AS 
$$
DECLARE
  project_team_id UUID;
BEGIN
  -- Get the team_id for the project
  SELECT team_id INTO project_team_id
  FROM projects
  WHERE id = project_id;
  
  -- Check if user is an admin of this team
  RETURN auth.is_team_admin(project_team_id);
END;
$$ 
LANGUAGE plpgsql SECURITY DEFINER;