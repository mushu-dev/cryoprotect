-- User teams table policies
DROP POLICY IF EXISTS user_teams_select_policy ON user_teams;
CREATE POLICY user_teams_select_policy ON user_teams
    FOR SELECT USING (user_id = auth.uid() OR team_id IN (SELECT team_id FROM user_teams WHERE user_id = auth.uid() AND role = 'admin'));

DROP POLICY IF EXISTS user_teams_insert_policy ON user_teams;
CREATE POLICY user_teams_insert_policy ON user_teams
    FOR INSERT WITH CHECK (team_id IN (SELECT team_id FROM user_teams WHERE user_id = auth.uid() AND role = 'admin'));

DROP POLICY IF EXISTS user_teams_update_policy ON user_teams;
CREATE POLICY user_teams_update_policy ON user_teams
    FOR UPDATE USING (team_id IN (SELECT team_id FROM user_teams WHERE user_id = auth.uid() AND role = 'admin')) 
    WITH CHECK (team_id IN (SELECT team_id FROM user_teams WHERE user_id = auth.uid() AND role = 'admin'));

DROP POLICY IF EXISTS user_teams_delete_policy ON user_teams;
CREATE POLICY user_teams_delete_policy ON user_teams
    FOR DELETE USING (team_id IN (SELECT team_id FROM user_teams WHERE user_id = auth.uid() AND role = 'admin'));