-- Create policies for user_profile and experiment_properties
DROP POLICY IF EXISTS user_profile_select_policy ON user_profile;
CREATE POLICY user_profile_select_policy ON user_profile
    FOR SELECT USING (auth_user_id = auth.uid() OR auth.is_service_role());

DROP POLICY IF EXISTS user_profile_update_policy ON user_profile;
CREATE POLICY user_profile_update_policy ON user_profile
    FOR UPDATE USING (auth_user_id = auth.uid()) WITH CHECK (auth_user_id = auth.uid());

DROP POLICY IF EXISTS experiment_properties_select_policy ON experiment_properties;
CREATE POLICY experiment_properties_select_policy ON experiment_properties
    FOR SELECT USING (EXISTS (SELECT 1 FROM experiments e WHERE e.id = experiment_id AND auth.can_access_project(e.project_id)) OR auth.is_service_role());

DROP POLICY IF EXISTS experiment_properties_insert_policy ON experiment_properties;
CREATE POLICY experiment_properties_insert_policy ON experiment_properties
    FOR INSERT WITH CHECK (EXISTS (SELECT 1 FROM experiments e WHERE e.id = experiment_id AND auth.can_access_project(e.project_id)) OR auth.is_service_role());

DROP POLICY IF EXISTS experiment_properties_update_policy ON experiment_properties;
CREATE POLICY experiment_properties_update_policy ON experiment_properties
    FOR UPDATE USING (EXISTS (SELECT 1 FROM experiments e WHERE e.id = experiment_id AND auth.can_access_project(e.project_id)) OR auth.is_service_role())
    WITH CHECK (EXISTS (SELECT 1 FROM experiments e WHERE e.id = experiment_id AND auth.can_access_project(e.project_id)) OR auth.is_service_role());

DROP POLICY IF EXISTS experiment_properties_delete_policy ON experiment_properties;
CREATE POLICY experiment_properties_delete_policy ON experiment_properties
    FOR DELETE USING (EXISTS (SELECT 1 FROM experiments e WHERE e.id = experiment_id AND auth.can_access_project(e.project_id)) OR auth.is_service_role());