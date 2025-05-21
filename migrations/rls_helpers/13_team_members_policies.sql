-- Team members table policies
DROP POLICY IF EXISTS team_members_select_policy ON team_members;
CREATE POLICY team_members_select_policy ON team_members
    FOR SELECT USING (user_id = auth.uid() OR team_id IN (SELECT team_id FROM team_members WHERE user_id = auth.uid() AND role = 'admin'));

DROP POLICY IF EXISTS team_members_insert_policy ON team_members;
CREATE POLICY team_members_insert_policy ON team_members
    FOR INSERT WITH CHECK (team_id IN (SELECT team_id FROM team_members WHERE user_id = auth.uid() AND role = 'admin'));

DROP POLICY IF EXISTS team_members_update_policy ON team_members;
CREATE POLICY team_members_update_policy ON team_members
    FOR UPDATE USING (team_id IN (SELECT team_id FROM team_members WHERE user_id = auth.uid() AND role = 'admin')) 
    WITH CHECK (team_id IN (SELECT team_id FROM team_members WHERE user_id = auth.uid() AND role = 'admin'));

DROP POLICY IF EXISTS team_members_delete_policy ON team_members;
CREATE POLICY team_members_delete_policy ON team_members
    FOR DELETE USING (team_id IN (SELECT team_id FROM team_members WHERE user_id = auth.uid() AND role = 'admin'));