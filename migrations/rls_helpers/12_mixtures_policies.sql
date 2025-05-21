-- Mixtures table policies
DROP POLICY IF EXISTS mixtures_select_policy ON mixtures;
CREATE POLICY mixtures_select_policy ON mixtures
    FOR SELECT USING (auth.can_access_project(project_id));

DROP POLICY IF EXISTS mixtures_insert_policy ON mixtures;
CREATE POLICY mixtures_insert_policy ON mixtures
    FOR INSERT WITH CHECK (auth.can_access_project(project_id));

DROP POLICY IF EXISTS mixtures_update_policy ON mixtures;
CREATE POLICY mixtures_update_policy ON mixtures
    FOR UPDATE USING (auth.can_access_project(project_id)) WITH CHECK (auth.can_access_project(project_id));

DROP POLICY IF EXISTS mixtures_delete_policy ON mixtures;
CREATE POLICY mixtures_delete_policy ON mixtures
    FOR DELETE USING (auth.can_access_project(project_id));