-- Projects table policies
DROP POLICY IF EXISTS projects_select_policy ON projects;
CREATE POLICY projects_select_policy ON projects
    FOR SELECT USING (auth.can_access_project(id));

DROP POLICY IF EXISTS projects_insert_policy ON projects;
CREATE POLICY projects_insert_policy ON projects
    FOR INSERT WITH CHECK (auth.is_team_member(team_id));

DROP POLICY IF EXISTS projects_update_policy ON projects;
CREATE POLICY projects_update_policy ON projects
    FOR UPDATE USING (auth.can_manage_project(id)) WITH CHECK (auth.can_manage_project(id));

DROP POLICY IF EXISTS projects_delete_policy ON projects;
CREATE POLICY projects_delete_policy ON projects
    FOR DELETE USING (auth.can_manage_project(id));