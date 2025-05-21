-- Service role access policies
DROP POLICY IF EXISTS projects_service_role_policy ON projects;
CREATE POLICY projects_service_role_policy ON projects FOR ALL USING (auth.is_service_role());

DROP POLICY IF EXISTS experiments_service_role_policy ON experiments;
CREATE POLICY experiments_service_role_policy ON experiments FOR ALL USING (auth.is_service_role());

DROP POLICY IF EXISTS molecules_service_role_policy ON molecules;
CREATE POLICY molecules_service_role_policy ON molecules FOR ALL USING (auth.is_service_role());

DROP POLICY IF EXISTS molecular_properties_service_role_policy ON molecular_properties;
CREATE POLICY molecular_properties_service_role_policy ON molecular_properties FOR ALL USING (auth.is_service_role());

DROP POLICY IF EXISTS mixtures_service_role_policy ON mixtures;
CREATE POLICY mixtures_service_role_policy ON mixtures FOR ALL USING (auth.is_service_role());

DROP POLICY IF EXISTS teams_service_role_policy ON teams;
CREATE POLICY teams_service_role_policy ON teams FOR ALL USING (auth.is_service_role());

DROP POLICY IF EXISTS team_members_service_role_policy ON team_members;
CREATE POLICY team_members_service_role_policy ON team_members FOR ALL USING (auth.is_service_role());