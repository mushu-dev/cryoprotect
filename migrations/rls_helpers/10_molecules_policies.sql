-- Molecules table policies
DROP POLICY IF EXISTS molecules_select_policy ON molecules;
CREATE POLICY molecules_select_policy ON molecules
    FOR SELECT USING (EXISTS (SELECT 1 FROM experiments e WHERE e.molecule_id = id AND auth.can_access_project(e.project_id)));

DROP POLICY IF EXISTS molecules_insert_policy ON molecules;
CREATE POLICY molecules_insert_policy ON molecules
    FOR INSERT WITH CHECK (true);  -- Allow insert, association happens in experiments

DROP POLICY IF EXISTS molecules_update_policy ON molecules;
CREATE POLICY molecules_update_policy ON molecules
    FOR UPDATE USING (auth.can_access_molecule(id)) WITH CHECK (auth.can_access_molecule(id));

DROP POLICY IF EXISTS molecules_delete_policy ON molecules;
CREATE POLICY molecules_delete_policy ON molecules
    FOR DELETE USING (auth.can_access_molecule(id));