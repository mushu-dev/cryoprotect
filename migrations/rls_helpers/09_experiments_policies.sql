-- Experiments table policies
DROP POLICY IF EXISTS experiments_select_policy ON experiments;
CREATE POLICY experiments_select_policy ON experiments
    FOR SELECT USING (auth.can_access_project(project_id));

DROP POLICY IF EXISTS experiments_insert_policy ON experiments;
CREATE POLICY experiments_insert_policy ON experiments
    FOR INSERT WITH CHECK (auth.can_access_project(project_id));

DROP POLICY IF EXISTS experiments_update_policy ON experiments;
CREATE POLICY experiments_update_policy ON experiments
    FOR UPDATE USING (auth.can_access_project(project_id)) WITH CHECK (auth.can_access_project(project_id));

DROP POLICY IF EXISTS experiments_delete_policy ON experiments;
CREATE POLICY experiments_delete_policy ON experiments
    FOR DELETE USING (auth.can_access_project(project_id));