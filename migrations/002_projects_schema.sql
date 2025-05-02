-- Projects and Project Experiments Schema

-- Create projects table
CREATE TABLE IF NOT EXISTS projects (
    id UUID PRIMARY KEY DEFAULT uuid_generate_v4(),
    name TEXT NOT NULL,
    description TEXT,
    is_public BOOLEAN DEFAULT FALSE,
    created_at TIMESTAMP WITH TIME ZONE DEFAULT NOW(),
    updated_at TIMESTAMP WITH TIME ZONE DEFAULT NOW(),
    created_by UUID REFERENCES auth.users(id) ON DELETE SET NULL
);

-- Create project_experiments table for organizing experiments within projects
CREATE TABLE IF NOT EXISTS project_experiments (
    id UUID PRIMARY KEY DEFAULT uuid_generate_v4(),
    project_id UUID NOT NULL REFERENCES projects(id) ON DELETE CASCADE,
    experiment_id UUID NOT NULL REFERENCES experiments(id) ON DELETE CASCADE,
    notes TEXT,
    created_at TIMESTAMP WITH TIME ZONE DEFAULT NOW(),
    created_by UUID REFERENCES auth.users(id) ON DELETE SET NULL,
    UNIQUE(project_id, experiment_id)
);

-- Create function to update updated_at timestamp
CREATE OR REPLACE FUNCTION update_updated_at_column()
RETURNS TRIGGER AS $$
BEGIN
    NEW.updated_at = NOW();
    RETURN NEW;
END;
$$ LANGUAGE plpgsql;

-- Create trigger to update updated_at timestamp on projects table
CREATE TRIGGER update_projects_updated_at
BEFORE UPDATE ON projects
FOR EACH ROW
EXECUTE FUNCTION update_updated_at_column();

-- Create function to get project with experiment count
CREATE OR REPLACE FUNCTION get_project_with_experiment_count(p_project_id UUID)
RETURNS TABLE (
    id UUID,
    name TEXT,
    description TEXT,
    is_public BOOLEAN,
    created_at TIMESTAMP WITH TIME ZONE,
    updated_at TIMESTAMP WITH TIME ZONE,
    created_by UUID,
    experiment_count BIGINT
) AS $$
BEGIN
    RETURN QUERY
    SELECT 
        p.id,
        p.name,
        p.description,
        p.is_public,
        p.created_at,
        p.updated_at,
        p.created_by,
        COUNT(pe.id)::BIGINT AS experiment_count
    FROM 
        projects p
    LEFT JOIN 
        project_experiments pe ON p.id = pe.project_id
    WHERE 
        p.id = p_project_id
    GROUP BY 
        p.id;
END;
$$ LANGUAGE plpgsql;

-- Create function to get user projects
CREATE OR REPLACE FUNCTION get_user_projects(p_user_id UUID, p_limit INT DEFAULT 100, p_offset INT DEFAULT 0)
RETURNS TABLE (
    id UUID,
    name TEXT,
    description TEXT,
    is_public BOOLEAN,
    created_at TIMESTAMP WITH TIME ZONE,
    updated_at TIMESTAMP WITH TIME ZONE,
    created_by UUID,
    experiment_count BIGINT
) AS $$
BEGIN
    RETURN QUERY
    SELECT 
        p.id,
        p.name,
        p.description,
        p.is_public,
        p.created_at,
        p.updated_at,
        p.created_by,
        COUNT(pe.id)::BIGINT AS experiment_count
    FROM 
        projects p
    LEFT JOIN 
        project_experiments pe ON p.id = pe.project_id
    WHERE 
        p.created_by = p_user_id OR p.is_public = TRUE
    GROUP BY 
        p.id
    ORDER BY 
        p.updated_at DESC
    LIMIT p_limit
    OFFSET p_offset;
END;
$$ LANGUAGE plpgsql;

-- Create function to get project experiments
CREATE OR REPLACE FUNCTION get_project_experiments(p_project_id UUID)
RETURNS TABLE (
    id UUID,
    project_id UUID,
    experiment_id UUID,
    mixture_id UUID,
    property_type_id UUID,
    property_name TEXT,
    numeric_value DOUBLE PRECISION,
    text_value TEXT,
    boolean_value BOOLEAN,
    experimental_conditions TEXT,
    date_performed DATE,
    created_at TIMESTAMP WITH TIME ZONE,
    notes TEXT,
    name TEXT
) AS $$
BEGIN
    RETURN QUERY
    SELECT 
        pe.id,
        pe.project_id,
        pe.experiment_id,
        e.mixture_id,
        e.property_type_id,
        pt.name AS property_name,
        e.numeric_value,
        e.text_value,
        e.boolean_value,
        e.experimental_conditions,
        e.date_performed,
        pe.created_at,
        pe.notes,
        m.name
    FROM 
        project_experiments pe
    JOIN 
        experiments e ON pe.experiment_id = e.id
    JOIN 
        property_types pt ON e.property_type_id = pt.id
    LEFT JOIN
        mixtures m ON e.mixture_id = m.id
    WHERE 
        pe.project_id = p_project_id
    ORDER BY 
        e.date_performed DESC;
END;
$$ LANGUAGE plpgsql;

-- Create function to get project activity
CREATE OR REPLACE FUNCTION get_project_activity(p_project_id UUID, p_limit INT DEFAULT 10)
RETURNS TABLE (
    id UUID,
    project_id UUID,
    title TEXT,
    description TEXT,
    timestamp TIMESTAMP WITH TIME ZONE
) AS $$
BEGIN
    -- This is a placeholder function that would normally query an activity log table
    -- For now, we'll return some mock data
    RETURN QUERY
    WITH mock_activity AS (
        SELECT 
            uuid_generate_v4() AS id,
            p_project_id AS project_id,
            'Project Created' AS title,
            'Project was created' AS description,
            (SELECT created_at FROM projects WHERE id = p_project_id) AS timestamp
        UNION ALL
        SELECT 
            uuid_generate_v4() AS id,
            p_project_id AS project_id,
            'Experiment Added' AS title,
            'New experiment was added to the project' AS description,
            (SELECT created_at FROM project_experiments WHERE project_id = p_project_id ORDER BY created_at DESC LIMIT 1) AS timestamp
        WHERE EXISTS (SELECT 1 FROM project_experiments WHERE project_id = p_project_id)
        UNION ALL
        SELECT 
            uuid_generate_v4() AS id,
            p_project_id AS project_id,
            'Project Updated' AS title,
            'Project details were updated' AS description,
            (SELECT updated_at FROM projects WHERE id = p_project_id) AS timestamp
    )
    SELECT * FROM mock_activity
    WHERE timestamp IS NOT NULL
    ORDER BY timestamp DESC
    LIMIT p_limit;
END;
$$ LANGUAGE plpgsql;

-- Add RLS policies for projects table
ALTER TABLE projects ENABLE ROW LEVEL SECURITY;

CREATE POLICY projects_select_policy ON projects
    FOR SELECT
    USING (
        auth.uid() = created_by OR is_public = TRUE
    );

CREATE POLICY projects_insert_policy ON projects
    FOR INSERT
    WITH CHECK (
        auth.uid() = created_by
    );

CREATE POLICY projects_update_policy ON projects
    FOR UPDATE
    USING (
        auth.uid() = created_by
    );

CREATE POLICY projects_delete_policy ON projects
    FOR DELETE
    USING (
        auth.uid() = created_by
    );

-- Add RLS policies for project_experiments table
ALTER TABLE project_experiments ENABLE ROW LEVEL SECURITY;

CREATE POLICY project_experiments_select_policy ON project_experiments
    FOR SELECT
    USING (
        auth.uid() = created_by OR 
        EXISTS (
            SELECT 1 FROM projects p 
            WHERE p.id = project_id AND (p.created_by = auth.uid() OR p.is_public = TRUE)
        )
    );

CREATE POLICY project_experiments_insert_policy ON project_experiments
    FOR INSERT
    WITH CHECK (
        auth.uid() = created_by AND
        EXISTS (
            SELECT 1 FROM projects p 
            WHERE p.id = project_id AND p.created_by = auth.uid()
        )
    );

CREATE POLICY project_experiments_update_policy ON project_experiments
    FOR UPDATE
    USING (
        auth.uid() = created_by
    );

CREATE POLICY project_experiments_delete_policy ON project_experiments
    FOR DELETE
    USING (
        auth.uid() = created_by OR
        EXISTS (
            SELECT 1 FROM projects p 
            WHERE p.id = project_id AND p.created_by = auth.uid()
        )
    );