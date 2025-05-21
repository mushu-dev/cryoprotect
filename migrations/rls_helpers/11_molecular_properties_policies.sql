-- Molecular properties table policies
DROP POLICY IF EXISTS molecular_properties_select_policy ON molecular_properties;
CREATE POLICY molecular_properties_select_policy ON molecular_properties
    FOR SELECT USING (auth.can_access_molecule(molecule_id));

DROP POLICY IF EXISTS molecular_properties_insert_policy ON molecular_properties;
CREATE POLICY molecular_properties_insert_policy ON molecular_properties
    FOR INSERT WITH CHECK (auth.can_access_molecule(molecule_id));

DROP POLICY IF EXISTS molecular_properties_update_policy ON molecular_properties;
CREATE POLICY molecular_properties_update_policy ON molecular_properties
    FOR UPDATE USING (auth.can_access_molecule(molecule_id)) WITH CHECK (auth.can_access_molecule(molecule_id));

DROP POLICY IF EXISTS molecular_properties_delete_policy ON molecular_properties;
CREATE POLICY molecular_properties_delete_policy ON molecular_properties
    FOR DELETE USING (auth.can_access_molecule(molecule_id));