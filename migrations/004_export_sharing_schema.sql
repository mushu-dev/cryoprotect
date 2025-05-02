-- CryoProtect Analyzer - Export and Sharing Schema Migration

-- Create shares table for storing shared items
CREATE TABLE IF NOT EXISTS shares (
    id UUID PRIMARY KEY,
    user_id UUID NOT NULL REFERENCES auth.users(id) ON DELETE CASCADE,
    share_type TEXT NOT NULL CHECK (share_type IN ('link', 'email', 'embed')),
    data_type TEXT NOT NULL CHECK (data_type IN ('molecules', 'mixtures', 'predictions', 'experiments', 'comparisons', 'report', 'visualization')),
    item_id UUID NOT NULL,
    created_at TIMESTAMP WITH TIME ZONE DEFAULT NOW(),
    password_protected BOOLEAN DEFAULT FALSE,
    password_hash TEXT,
    expiration TIMESTAMP WITH TIME ZONE,
    access_count INTEGER DEFAULT 0
);

-- Create report_templates table for storing report templates
CREATE TABLE IF NOT EXISTS report_templates (
    id UUID PRIMARY KEY DEFAULT uuid_generate_v4(),
    name TEXT NOT NULL,
    description TEXT,
    template_data JSONB NOT NULL,
    created_by UUID NOT NULL REFERENCES auth.users(id) ON DELETE CASCADE,
    created_at TIMESTAMP WITH TIME ZONE DEFAULT NOW(),
    updated_at TIMESTAMP WITH TIME ZONE DEFAULT NOW()
);

-- Create visualization_templates table for storing visualization templates
CREATE TABLE IF NOT EXISTS visualization_templates (
    id UUID PRIMARY KEY DEFAULT uuid_generate_v4(),
    name TEXT NOT NULL,
    description TEXT,
    chart_type TEXT NOT NULL,
    template_data JSONB NOT NULL,
    created_by UUID NOT NULL REFERENCES auth.users(id) ON DELETE CASCADE,
    created_at TIMESTAMP WITH TIME ZONE DEFAULT NOW(),
    updated_at TIMESTAMP WITH TIME ZONE DEFAULT NOW()
);

-- Create export_logs table for tracking exports
CREATE TABLE IF NOT EXISTS export_logs (
    id UUID PRIMARY KEY DEFAULT uuid_generate_v4(),
    user_id UUID NOT NULL REFERENCES auth.users(id) ON DELETE CASCADE,
    data_type TEXT NOT NULL,
    item_id UUID,
    format TEXT NOT NULL,
    created_at TIMESTAMP WITH TIME ZONE DEFAULT NOW(),
    success BOOLEAN DEFAULT TRUE,
    error_message TEXT
);

-- Create share_access_logs table for tracking share access
CREATE TABLE IF NOT EXISTS share_access_logs (
    id UUID PRIMARY KEY DEFAULT uuid_generate_v4(),
    share_id UUID NOT NULL REFERENCES shares(id) ON DELETE CASCADE,
    accessed_at TIMESTAMP WITH TIME ZONE DEFAULT NOW(),
    accessed_by UUID REFERENCES auth.users(id) ON DELETE SET NULL,
    ip_address TEXT,
    user_agent TEXT
);

-- Create function to update shares access_count
CREATE OR REPLACE FUNCTION update_share_access_count()
RETURNS TRIGGER AS $$
BEGIN
    UPDATE shares
    SET access_count = access_count + 1
    WHERE id = NEW.share_id;
    RETURN NEW;
END;
$$ LANGUAGE plpgsql;

-- Create trigger to update shares access_count
CREATE TRIGGER update_share_access_count_trigger
AFTER INSERT ON share_access_logs
FOR EACH ROW
EXECUTE FUNCTION update_share_access_count();

-- Create function to check if a share has expired
CREATE OR REPLACE FUNCTION is_share_expired(share_id UUID)
RETURNS BOOLEAN AS $$
DECLARE
    share_expiration TIMESTAMP WITH TIME ZONE;
BEGIN
    SELECT expiration INTO share_expiration
    FROM shares
    WHERE id = share_id;
    
    IF share_expiration IS NULL THEN
        RETURN FALSE;
    ELSE
        RETURN NOW() > share_expiration;
    END IF;
END;
$$ LANGUAGE plpgsql;

-- Create function to get shared item data
CREATE OR REPLACE FUNCTION get_shared_item_data(p_share_id UUID)
RETURNS JSONB AS $$
DECLARE
    v_share shares;
    v_data JSONB;
BEGIN
    -- Get share record
    SELECT * INTO v_share
    FROM shares
    WHERE id = p_share_id;
    
    IF v_share IS NULL THEN
        RAISE EXCEPTION 'Share not found';
    END IF;
    
    -- Check if share has expired
    IF is_share_expired(p_share_id) THEN
        RAISE EXCEPTION 'Share has expired';
    END IF;
    
    -- Get data based on data_type
    CASE v_share.data_type
        WHEN 'molecules' THEN
            SELECT row_to_json(m)::JSONB INTO v_data
            FROM molecule_with_properties m
            WHERE m.id = v_share.item_id;
            
        WHEN 'mixtures' THEN
            SELECT row_to_json(m)::JSONB INTO v_data
            FROM mixture_with_components m
            WHERE m.id = v_share.item_id;
            
        WHEN 'predictions' THEN
            SELECT jsonb_agg(p)::JSONB INTO v_data
            FROM (
                SELECT p.*, pt.name as property_name, cm.name as calculation_method
                FROM predictions p
                JOIN property_types pt ON p.property_type_id = pt.id
                JOIN calculation_methods cm ON p.calculation_method_id = cm.id
                WHERE p.mixture_id = v_share.item_id
            ) p;
            
        WHEN 'experiments' THEN
            SELECT jsonb_agg(e)::JSONB INTO v_data
            FROM (
                SELECT e.*, pt.name as property_name
                FROM experiments e
                JOIN property_types pt ON e.property_type_id = pt.id
                WHERE e.mixture_id = v_share.item_id
            ) e;
            
        WHEN 'comparisons' THEN
            SELECT jsonb_agg(c)::JSONB INTO v_data
            FROM (
                SELECT *
                FROM compare_prediction_with_experiment(v_share.item_id, pt.name)
                FROM property_types pt
            ) c;
            
        WHEN 'report' THEN
            -- For reports, the item_id is the report template id
            SELECT template_data INTO v_data
            FROM report_templates
            WHERE id = v_share.item_id;
            
        WHEN 'visualization' THEN
            -- For visualizations, the item_id is the visualization template id
            SELECT template_data INTO v_data
            FROM visualization_templates
            WHERE id = v_share.item_id;
            
        ELSE
            RAISE EXCEPTION 'Unsupported data type: %', v_share.data_type;
    END CASE;
    
    IF v_data IS NULL THEN
        RAISE EXCEPTION 'Shared item data not found';
    END IF;
    
    -- Log access
    INSERT INTO share_access_logs (share_id)
    VALUES (p_share_id);
    
    RETURN v_data;
END;
$$ LANGUAGE plpgsql;

-- Create RLS policies for shares table
ALTER TABLE shares ENABLE ROW LEVEL SECURITY;

CREATE POLICY shares_select_policy ON shares
    FOR SELECT
    USING (user_id = auth.uid() OR id IN (
        SELECT share_id FROM share_access_logs WHERE accessed_by = auth.uid()
    ));

CREATE POLICY shares_insert_policy ON shares
    FOR INSERT
    WITH CHECK (user_id = auth.uid());

CREATE POLICY shares_update_policy ON shares
    FOR UPDATE
    USING (user_id = auth.uid());

CREATE POLICY shares_delete_policy ON shares
    FOR DELETE
    USING (user_id = auth.uid());

-- Create RLS policies for report_templates table
ALTER TABLE report_templates ENABLE ROW LEVEL SECURITY;

CREATE POLICY report_templates_select_policy ON report_templates
    FOR SELECT
    USING (created_by = auth.uid());

CREATE POLICY report_templates_insert_policy ON report_templates
    FOR INSERT
    WITH CHECK (created_by = auth.uid());

CREATE POLICY report_templates_update_policy ON report_templates
    FOR UPDATE
    USING (created_by = auth.uid());

CREATE POLICY report_templates_delete_policy ON report_templates
    FOR DELETE
    USING (created_by = auth.uid());

-- Create RLS policies for visualization_templates table
ALTER TABLE visualization_templates ENABLE ROW LEVEL SECURITY;

CREATE POLICY visualization_templates_select_policy ON visualization_templates
    FOR SELECT
    USING (created_by = auth.uid());

CREATE POLICY visualization_templates_insert_policy ON visualization_templates
    FOR INSERT
    WITH CHECK (created_by = auth.uid());

CREATE POLICY visualization_templates_update_policy ON visualization_templates
    FOR UPDATE
    USING (created_by = auth.uid());

CREATE POLICY visualization_templates_delete_policy ON visualization_templates
    FOR DELETE
    USING (created_by = auth.uid());

-- Create RLS policies for export_logs table
ALTER TABLE export_logs ENABLE ROW LEVEL SECURITY;

CREATE POLICY export_logs_select_policy ON export_logs
    FOR SELECT
    USING (user_id = auth.uid());

CREATE POLICY export_logs_insert_policy ON export_logs
    FOR INSERT
    WITH CHECK (user_id = auth.uid());

-- Create RLS policies for share_access_logs table
ALTER TABLE share_access_logs ENABLE ROW LEVEL SECURITY;

CREATE POLICY share_access_logs_select_policy ON share_access_logs
    FOR SELECT
    USING (
        share_id IN (SELECT id FROM shares WHERE user_id = auth.uid()) OR
        accessed_by = auth.uid()
    );

CREATE POLICY share_access_logs_insert_policy ON share_access_logs
    FOR INSERT
    WITH CHECK (TRUE);  -- Allow all inserts, will be controlled by application logic