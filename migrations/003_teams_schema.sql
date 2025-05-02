-- Teams and Collaboration Schema Migration
-- Compatible with Supabase and PostgreSQL

-- Create teams table
CREATE TABLE IF NOT EXISTS public.teams (
    id UUID PRIMARY KEY DEFAULT uuid_generate_v4(),
    name TEXT NOT NULL,
    description TEXT,
    created_at TIMESTAMPTZ NOT NULL DEFAULT NOW(),
    updated_at TIMESTAMPTZ NOT NULL DEFAULT NOW(),
    created_by UUID REFERENCES auth.users(id) ON DELETE SET NULL
);

-- Create team_members table for managing team membership
CREATE TABLE IF NOT EXISTS public.team_members (
    id UUID PRIMARY KEY DEFAULT uuid_generate_v4(),
    team_id UUID NOT NULL REFERENCES public.teams(id) ON DELETE CASCADE,
    user_id UUID NOT NULL REFERENCES auth.users(id) ON DELETE CASCADE,
    role TEXT NOT NULL CHECK (role IN ('admin', 'editor', 'viewer')),
    created_at TIMESTAMPTZ NOT NULL DEFAULT NOW(),
    updated_at TIMESTAMPTZ NOT NULL DEFAULT NOW(),
    created_by UUID REFERENCES auth.users(id) ON DELETE SET NULL,
    UNIQUE(team_id, user_id)
);

-- Create shared_resources table for tracking which resources are shared with which teams
CREATE TABLE IF NOT EXISTS public.shared_resources (
    id UUID PRIMARY KEY DEFAULT uuid_generate_v4(),
    team_id UUID NOT NULL REFERENCES public.teams(id) ON DELETE CASCADE,
    resource_type TEXT NOT NULL CHECK (resource_type IN ('project', 'mixture', 'experiment')),
    resource_id UUID NOT NULL,
    created_at TIMESTAMPTZ NOT NULL DEFAULT NOW(),
    created_by UUID REFERENCES auth.users(id) ON DELETE SET NULL,
    UNIQUE(team_id, resource_type, resource_id)
);

-- Create comments table for discussions
CREATE TABLE IF NOT EXISTS public.comments (
    id UUID PRIMARY KEY DEFAULT uuid_generate_v4(),
    resource_type TEXT NOT NULL CHECK (resource_type IN ('project', 'mixture', 'experiment')),
    resource_id UUID NOT NULL,
    content TEXT NOT NULL,
    parent_id UUID REFERENCES public.comments(id) ON DELETE CASCADE,
    created_at TIMESTAMPTZ NOT NULL DEFAULT NOW(),
    updated_at TIMESTAMPTZ NOT NULL DEFAULT NOW(),
    created_by UUID REFERENCES auth.users(id) ON DELETE SET NULL
);

-- Create notifications table
CREATE TABLE IF NOT EXISTS public.notifications (
    id UUID PRIMARY KEY DEFAULT uuid_generate_v4(),
    user_id UUID NOT NULL REFERENCES auth.users(id) ON DELETE CASCADE,
    team_id UUID REFERENCES public.teams(id) ON DELETE CASCADE,
    title TEXT NOT NULL,
    message TEXT NOT NULL,
    resource_type TEXT CHECK (resource_type IN ('team', 'project', 'mixture', 'experiment', 'comment')),
    resource_id UUID,
    is_read BOOLEAN NOT NULL DEFAULT FALSE,
    created_at TIMESTAMPTZ NOT NULL DEFAULT NOW(),
    created_by UUID REFERENCES auth.users(id) ON DELETE SET NULL
);

-- Create activity_log table for tracking team activities
CREATE TABLE IF NOT EXISTS public.activity_log (
    id UUID PRIMARY KEY DEFAULT uuid_generate_v4(),
    team_id UUID REFERENCES public.teams(id) ON DELETE CASCADE,
    user_id UUID REFERENCES auth.users(id) ON DELETE SET NULL,
    action TEXT NOT NULL,
    resource_type TEXT NOT NULL CHECK (resource_type IN ('team', 'project', 'mixture', 'experiment', 'comment')),
    resource_id UUID,
    details JSONB,
    created_at TIMESTAMPTZ NOT NULL DEFAULT NOW()
);

-- Create triggers for updated_at timestamps
CREATE TRIGGER set_timestamp_teams
BEFORE UPDATE ON public.teams
FOR EACH ROW EXECUTE PROCEDURE public.trigger_set_timestamp();

CREATE TRIGGER set_timestamp_team_members
BEFORE UPDATE ON public.team_members
FOR EACH ROW EXECUTE PROCEDURE public.trigger_set_timestamp();

CREATE TRIGGER set_timestamp_comments
BEFORE UPDATE ON public.comments
FOR EACH ROW EXECUTE PROCEDURE public.trigger_set_timestamp();

-- Enable Row Level Security (RLS)
ALTER TABLE public.teams ENABLE ROW LEVEL SECURITY;
ALTER TABLE public.team_members ENABLE ROW LEVEL SECURITY;
ALTER TABLE public.shared_resources ENABLE ROW LEVEL SECURITY;
ALTER TABLE public.comments ENABLE ROW LEVEL SECURITY;
ALTER TABLE public.notifications ENABLE ROW LEVEL SECURITY;
ALTER TABLE public.activity_log ENABLE ROW LEVEL SECURITY;

-- Create RLS Policies

-- Teams: viewable by members, insertable by authenticated users, updatable/deletable by creator or admin
CREATE POLICY "Teams are viewable by members" 
ON public.teams FOR SELECT USING (
    auth.uid() = created_by OR 
    EXISTS (
        SELECT 1 FROM public.team_members 
        WHERE team_id = id AND user_id = auth.uid()
    )
);

CREATE POLICY "Teams can be inserted by authenticated users" 
ON public.teams FOR INSERT WITH CHECK (auth.role() = 'authenticated');

CREATE POLICY "Teams can be updated by creator or admin" 
ON public.teams FOR UPDATE USING (
    auth.uid() = created_by OR 
    EXISTS (
        SELECT 1 FROM public.team_members 
        WHERE team_id = id AND user_id = auth.uid() AND role = 'admin'
    )
);

CREATE POLICY "Teams can be deleted by creator or admin" 
ON public.teams FOR DELETE USING (
    auth.uid() = created_by OR 
    EXISTS (
        SELECT 1 FROM public.team_members 
        WHERE team_id = id AND user_id = auth.uid() AND role = 'admin'
    )
);

-- Team Members: viewable by team members, insertable/updatable/deletable by team admin
CREATE POLICY "Team members are viewable by team members" 
ON public.team_members FOR SELECT USING (
    user_id = auth.uid() OR 
    EXISTS (
        SELECT 1 FROM public.team_members tm 
        WHERE tm.team_id = team_id AND tm.user_id = auth.uid()
    )
);

CREATE POLICY "Team members can be inserted by team admin" 
ON public.team_members FOR INSERT WITH CHECK (
    EXISTS (
        SELECT 1 FROM public.teams t 
        WHERE t.id = team_id AND t.created_by = auth.uid()
    ) OR
    EXISTS (
        SELECT 1 FROM public.team_members tm 
        WHERE tm.team_id = team_id AND tm.user_id = auth.uid() AND tm.role = 'admin'
    )
);

CREATE POLICY "Team members can be updated by team admin" 
ON public.team_members FOR UPDATE USING (
    EXISTS (
        SELECT 1 FROM public.teams t 
        WHERE t.id = team_id AND t.created_by = auth.uid()
    ) OR
    EXISTS (
        SELECT 1 FROM public.team_members tm 
        WHERE tm.team_id = team_id AND tm.user_id = auth.uid() AND tm.role = 'admin'
    )
);

CREATE POLICY "Team members can be deleted by team admin" 
ON public.team_members FOR DELETE USING (
    EXISTS (
        SELECT 1 FROM public.teams t 
        WHERE t.id = team_id AND t.created_by = auth.uid()
    ) OR
    EXISTS (
        SELECT 1 FROM public.team_members tm 
        WHERE tm.team_id = team_id AND tm.user_id = auth.uid() AND tm.role = 'admin'
    )
);

-- Shared Resources: viewable by team members, insertable/deletable by resource owner or team admin
CREATE POLICY "Shared resources are viewable by team members" 
ON public.shared_resources FOR SELECT USING (
    EXISTS (
        SELECT 1 FROM public.team_members 
        WHERE team_id = team_id AND user_id = auth.uid()
    )
);

CREATE POLICY "Shared resources can be inserted by resource owner or team admin" 
ON public.shared_resources FOR INSERT WITH CHECK (
    auth.uid() = created_by OR
    EXISTS (
        SELECT 1 FROM public.team_members 
        WHERE team_id = team_id AND user_id = auth.uid() AND role = 'admin'
    )
);

CREATE POLICY "Shared resources can be deleted by resource owner or team admin" 
ON public.shared_resources FOR DELETE USING (
    auth.uid() = created_by OR
    EXISTS (
        SELECT 1 FROM public.team_members 
        WHERE team_id = team_id AND user_id = auth.uid() AND role = 'admin'
    )
);

-- Comments: viewable by resource viewers, insertable by authenticated users, updatable/deletable by creator
CREATE POLICY "Comments are viewable by resource viewers" 
ON public.comments FOR SELECT USING (true);

CREATE POLICY "Comments can be inserted by authenticated users" 
ON public.comments FOR INSERT WITH CHECK (auth.role() = 'authenticated');

CREATE POLICY "Comments can be updated by creator" 
ON public.comments FOR UPDATE USING (auth.uid() = created_by);

CREATE POLICY "Comments can be deleted by creator" 
ON public.comments FOR DELETE USING (auth.uid() = created_by);

-- Notifications: viewable/updatable/deletable by recipient
CREATE POLICY "Notifications are viewable by recipient" 
ON public.notifications FOR SELECT USING (user_id = auth.uid());

CREATE POLICY "Notifications can be updated by recipient" 
ON public.notifications FOR UPDATE USING (user_id = auth.uid());

CREATE POLICY "Notifications can be deleted by recipient" 
ON public.notifications FOR DELETE USING (user_id = auth.uid());

-- Activity Log: viewable by team members
CREATE POLICY "Activity log is viewable by team members" 
ON public.activity_log FOR SELECT USING (
    EXISTS (
        SELECT 1 FROM public.team_members 
        WHERE team_id = team_id AND user_id = auth.uid()
    )
);

-- Create indexes for performance
CREATE INDEX idx_team_members_team_id ON public.team_members(team_id);
CREATE INDEX idx_team_members_user_id ON public.team_members(user_id);
CREATE INDEX idx_shared_resources_team_id ON public.shared_resources(team_id);
CREATE INDEX idx_shared_resources_resource_id ON public.shared_resources(resource_id);
CREATE INDEX idx_comments_resource_id ON public.comments(resource_id);
CREATE INDEX idx_comments_parent_id ON public.comments(parent_id);
CREATE INDEX idx_notifications_user_id ON public.notifications(user_id);
CREATE INDEX idx_notifications_team_id ON public.notifications(team_id);
CREATE INDEX idx_activity_log_team_id ON public.activity_log(team_id);
CREATE INDEX idx_activity_log_user_id ON public.activity_log(user_id);

-- Create functions for team operations

-- Function to check if a user is a team member
CREATE OR REPLACE FUNCTION is_team_member(p_team_id UUID, p_user_id UUID)
RETURNS BOOLEAN AS $$
BEGIN
    RETURN EXISTS (
        SELECT 1 FROM public.team_members
        WHERE team_id = p_team_id AND user_id = p_user_id
    );
END;
$$ LANGUAGE plpgsql;

-- Function to check if a user has a specific role in a team
CREATE OR REPLACE FUNCTION has_team_role(p_team_id UUID, p_user_id UUID, p_role TEXT)
RETURNS BOOLEAN AS $$
BEGIN
    RETURN EXISTS (
        SELECT 1 FROM public.team_members
        WHERE team_id = p_team_id AND user_id = p_user_id AND role = p_role
    );
END;
$$ LANGUAGE plpgsql;

-- Function to get team members with user details
CREATE OR REPLACE FUNCTION get_team_members_with_details(p_team_id UUID)
RETURNS TABLE (
    id UUID,
    team_id UUID,
    user_id UUID,
    role TEXT,
    created_at TIMESTAMPTZ,
    email TEXT,
    display_name TEXT
) AS $$
BEGIN
    RETURN QUERY
    SELECT 
        tm.id,
        tm.team_id,
        tm.user_id,
        tm.role,
        tm.created_at,
        u.email,
        u.raw_user_meta_data->>'name' AS display_name
    FROM 
        public.team_members tm
    JOIN 
        auth.users u ON tm.user_id = u.id
    WHERE 
        tm.team_id = p_team_id
    ORDER BY 
        tm.role, u.email;
END;
$$ LANGUAGE plpgsql;

-- Function to get user teams with details
CREATE OR REPLACE FUNCTION get_user_teams(p_user_id UUID)
RETURNS TABLE (
    team_id UUID,
    team_name TEXT,
    team_description TEXT,
    role TEXT,
    member_count BIGINT
) AS $$
BEGIN
    RETURN QUERY
    SELECT 
        t.id AS team_id,
        t.name AS team_name,
        t.description AS team_description,
        tm.role,
        COUNT(tm2.id) AS member_count
    FROM 
        public.teams t
    JOIN 
        public.team_members tm ON t.id = tm.team_id AND tm.user_id = p_user_id
    LEFT JOIN 
        public.team_members tm2 ON t.id = tm2.team_id
    GROUP BY 
        t.id, t.name, t.description, tm.role
    ORDER BY 
        t.name;
END;
$$ LANGUAGE plpgsql;

-- Function to get team activity
CREATE OR REPLACE FUNCTION get_team_activity(p_team_id UUID, p_limit INT DEFAULT 20)
RETURNS TABLE (
    id UUID,
    team_id UUID,
    user_id UUID,
    user_email TEXT,
    action TEXT,
    resource_type TEXT,
    resource_id UUID,
    details JSONB,
    created_at TIMESTAMPTZ
) AS $$
BEGIN
    RETURN QUERY
    SELECT 
        al.id,
        al.team_id,
        al.user_id,
        u.email AS user_email,
        al.action,
        al.resource_type,
        al.resource_id,
        al.details,
        al.created_at
    FROM 
        public.activity_log al
    LEFT JOIN 
        auth.users u ON al.user_id = u.id
    WHERE 
        al.team_id = p_team_id
    ORDER BY 
        al.created_at DESC
    LIMIT p_limit;
END;
$$ LANGUAGE plpgsql;

-- Function to log team activity
CREATE OR REPLACE FUNCTION log_team_activity(
    p_team_id UUID,
    p_user_id UUID,
    p_action TEXT,
    p_resource_type TEXT,
    p_resource_id UUID,
    p_details JSONB DEFAULT NULL
) RETURNS UUID AS $$
DECLARE
    v_activity_id UUID;
BEGIN
    INSERT INTO public.activity_log (
        team_id,
        user_id,
        action,
        resource_type,
        resource_id,
        details
    ) VALUES (
        p_team_id,
        p_user_id,
        p_action,
        p_resource_type,
        p_resource_id,
        p_details
    ) RETURNING id INTO v_activity_id;
    
    RETURN v_activity_id;
END;
$$ LANGUAGE plpgsql;

-- Function to create a notification for team members
CREATE OR REPLACE FUNCTION create_team_notification(
    p_team_id UUID,
    p_title TEXT,
    p_message TEXT,
    p_resource_type TEXT,
    p_resource_id UUID,
    p_created_by UUID
) RETURNS VOID AS $$
BEGIN
    INSERT INTO public.notifications (
        user_id,
        team_id,
        title,
        message,
        resource_type,
        resource_id,
        created_by
    )
    SELECT 
        user_id,
        p_team_id,
        p_title,
        p_message,
        p_resource_type,
        p_resource_id,
        p_created_by
    FROM 
        public.team_members
    WHERE 
        team_id = p_team_id AND user_id != p_created_by;
END;
$$ LANGUAGE plpgsql;

-- Add comments
COMMENT ON TABLE public.teams IS 'Teams for collaborative research';
COMMENT ON TABLE public.team_members IS 'Team membership with role-based permissions';
COMMENT ON TABLE public.shared_resources IS 'Resources shared with teams';
COMMENT ON TABLE public.comments IS 'Comments and discussions on resources';
COMMENT ON TABLE public.notifications IS 'User notifications for team activities';
COMMENT ON TABLE public.activity_log IS 'Log of team activities';