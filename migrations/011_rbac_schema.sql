-- CryoProtect v2 - Role-Based Access Control Schema
-- Migration 011: RBAC Schema

-- Enable RLS on all tables
ALTER TABLE IF EXISTS roles ENABLE ROW LEVEL SECURITY;
ALTER TABLE IF EXISTS permissions ENABLE ROW LEVEL SECURITY;
ALTER TABLE IF EXISTS role_permissions ENABLE ROW LEVEL SECURITY;
ALTER TABLE IF EXISTS user_roles ENABLE ROW LEVEL SECURITY;

-- Create roles table
CREATE TABLE IF NOT EXISTS roles (
    id UUID PRIMARY KEY DEFAULT uuid_generate_v4(),
    name VARCHAR(50) NOT NULL UNIQUE,
    description TEXT,
    is_system_role BOOLEAN DEFAULT FALSE,
    created_at TIMESTAMP WITH TIME ZONE DEFAULT NOW(),
    updated_at TIMESTAMP WITH TIME ZONE DEFAULT NOW(),
    created_by UUID REFERENCES auth.users(id) ON DELETE SET NULL
);

-- Create permissions table
CREATE TABLE IF NOT EXISTS permissions (
    id UUID PRIMARY KEY DEFAULT uuid_generate_v4(),
    name VARCHAR(100) NOT NULL UNIQUE,
    description TEXT,
    resource_type VARCHAR(50) NOT NULL,
    action VARCHAR(50) NOT NULL,
    created_at TIMESTAMP WITH TIME ZONE DEFAULT NOW(),
    updated_at TIMESTAMP WITH TIME ZONE DEFAULT NOW(),
    created_by UUID REFERENCES auth.users(id) ON DELETE SET NULL,
    UNIQUE(resource_type, action)
);

-- Create role_permissions junction table
CREATE TABLE IF NOT EXISTS role_permissions (
    id UUID PRIMARY KEY DEFAULT uuid_generate_v4(),
    role_id UUID NOT NULL REFERENCES roles(id) ON DELETE CASCADE,
    permission_id UUID NOT NULL REFERENCES permissions(id) ON DELETE CASCADE,
    created_at TIMESTAMP WITH TIME ZONE DEFAULT NOW(),
    created_by UUID REFERENCES auth.users(id) ON DELETE SET NULL,
    UNIQUE(role_id, permission_id)
);

-- Create user_roles junction table
CREATE TABLE IF NOT EXISTS user_roles (
    id UUID PRIMARY KEY DEFAULT uuid_generate_v4(),
    user_id UUID NOT NULL REFERENCES auth.users(id) ON DELETE CASCADE,
    role_id UUID NOT NULL REFERENCES roles(id) ON DELETE CASCADE,
    created_at TIMESTAMP WITH TIME ZONE DEFAULT NOW(),
    updated_at TIMESTAMP WITH TIME ZONE DEFAULT NOW(),
    created_by UUID REFERENCES auth.users(id) ON DELETE SET NULL,
    UNIQUE(user_id, role_id)
);

-- Create function to check if a user has a specific role
CREATE OR REPLACE FUNCTION has_role(p_user_id UUID, p_role_name VARCHAR)
RETURNS BOOLEAN AS $$
BEGIN
    RETURN EXISTS (
        SELECT 1
        FROM user_roles ur
        JOIN roles r ON ur.role_id = r.id
        WHERE ur.user_id = p_user_id AND r.name = p_role_name
    );
END;
$$ LANGUAGE plpgsql SECURITY DEFINER;

-- Create function to check if a user has a specific permission
CREATE OR REPLACE FUNCTION has_permission(p_user_id UUID, p_resource_type VARCHAR, p_action VARCHAR)
RETURNS BOOLEAN AS $$
BEGIN
    RETURN EXISTS (
        SELECT 1
        FROM user_roles ur
        JOIN role_permissions rp ON ur.role_id = rp.role_id
        JOIN permissions p ON rp.permission_id = p.id
        WHERE ur.user_id = p_user_id
        AND p.resource_type = p_resource_type
        AND p.action = p_action
    );
END;
$$ LANGUAGE plpgsql SECURITY DEFINER;

-- Create function to get all roles for a user
CREATE OR REPLACE FUNCTION get_user_roles(p_user_id UUID)
RETURNS TABLE (role_name VARCHAR, role_id UUID) AS $$
BEGIN
    RETURN QUERY
    SELECT r.name, r.id
    FROM user_roles ur
    JOIN roles r ON ur.role_id = r.id
    WHERE ur.user_id = p_user_id;
END;
$$ LANGUAGE plpgsql SECURITY DEFINER;

-- Create function to get all permissions for a user
CREATE OR REPLACE FUNCTION get_user_permissions(p_user_id UUID)
RETURNS TABLE (resource_type VARCHAR, action VARCHAR, permission_id UUID) AS $$
BEGIN
    RETURN QUERY
    SELECT DISTINCT p.resource_type, p.action, p.id
    FROM user_roles ur
    JOIN role_permissions rp ON ur.role_id = rp.role_id
    JOIN permissions p ON rp.permission_id = p.id
    WHERE ur.user_id = p_user_id;
END;
$$ LANGUAGE plpgsql SECURITY DEFINER;

-- RLS Policies for roles table
CREATE POLICY "Roles are viewable by authenticated users"
ON roles FOR SELECT
USING (auth.role() = 'authenticated');

CREATE POLICY "Roles can be inserted by admins"
ON roles FOR INSERT
WITH CHECK (has_role(auth.uid(), 'admin'));

CREATE POLICY "Roles can be updated by admins"
ON roles FOR UPDATE
USING (has_role(auth.uid(), 'admin'));

CREATE POLICY "Roles can be deleted by admins"
ON roles FOR DELETE
USING (has_role(auth.uid(), 'admin') AND NOT is_system_role);

-- RLS Policies for permissions table
CREATE POLICY "Permissions are viewable by authenticated users"
ON permissions FOR SELECT
USING (auth.role() = 'authenticated');

CREATE POLICY "Permissions can be inserted by admins"
ON permissions FOR INSERT
WITH CHECK (has_role(auth.uid(), 'admin'));

CREATE POLICY "Permissions can be updated by admins"
ON permissions FOR UPDATE
USING (has_role(auth.uid(), 'admin'));

CREATE POLICY "Permissions can be deleted by admins"
ON permissions FOR DELETE
USING (has_role(auth.uid(), 'admin'));

-- RLS Policies for role_permissions table
CREATE POLICY "Role permissions are viewable by authenticated users"
ON role_permissions FOR SELECT
USING (auth.role() = 'authenticated');

CREATE POLICY "Role permissions can be inserted by admins"
ON role_permissions FOR INSERT
WITH CHECK (has_role(auth.uid(), 'admin'));

CREATE POLICY "Role permissions can be deleted by admins"
ON role_permissions FOR DELETE
USING (has_role(auth.uid(), 'admin'));

-- RLS Policies for user_roles table
CREATE POLICY "User roles are viewable by authenticated users"
ON user_roles FOR SELECT
USING (auth.role() = 'authenticated');

CREATE POLICY "User roles can be inserted by admins"
ON user_roles FOR INSERT
WITH CHECK (has_role(auth.uid(), 'admin'));

CREATE POLICY "User roles can be deleted by admins"
ON user_roles FOR DELETE
USING (has_role(auth.uid(), 'admin'));

-- Insert default roles
INSERT INTO roles (name, description, is_system_role, created_by)
VALUES 
('admin', 'Administrator with full system access', TRUE, auth.uid()),
('user', 'Regular user with standard access', TRUE, auth.uid()),
('curator', 'User with special privileges for curating data', TRUE, auth.uid()),
('viewer', 'User with read-only access to data', TRUE, auth.uid())
ON CONFLICT (name) DO NOTHING;

-- Insert default permissions
-- Admin permissions
INSERT INTO permissions (name, description, resource_type, action, created_by)
VALUES 
('admin:all', 'Full administrative access', 'admin', 'all', auth.uid()),
('users:manage', 'Manage users', 'users', 'manage', auth.uid()),
('roles:manage', 'Manage roles', 'roles', 'manage', auth.uid()),
('permissions:manage', 'Manage permissions', 'permissions', 'manage', auth.uid()),
('system:configure', 'Configure system settings', 'system', 'configure', auth.uid())
ON CONFLICT (resource_type, action) DO NOTHING;

-- Molecule permissions
INSERT INTO permissions (name, description, resource_type, action, created_by)
VALUES 
('molecules:create', 'Create molecules', 'molecules', 'create', auth.uid()),
('molecules:read', 'Read molecules', 'molecules', 'read', auth.uid()),
('molecules:update', 'Update molecules', 'molecules', 'update', auth.uid()),
('molecules:delete', 'Delete molecules', 'molecules', 'delete', auth.uid())
ON CONFLICT (resource_type, action) DO NOTHING;

-- Mixture permissions
INSERT INTO permissions (name, description, resource_type, action, created_by)
VALUES 
('mixtures:create', 'Create mixtures', 'mixtures', 'create', auth.uid()),
('mixtures:read', 'Read mixtures', 'mixtures', 'read', auth.uid()),
('mixtures:update', 'Update mixtures', 'mixtures', 'update', auth.uid()),
('mixtures:delete', 'Delete mixtures', 'mixtures', 'delete', auth.uid())
ON CONFLICT (resource_type, action) DO NOTHING;

-- Experiment permissions
INSERT INTO permissions (name, description, resource_type, action, created_by)
VALUES 
('experiments:create', 'Create experiments', 'experiments', 'create', auth.uid()),
('experiments:read', 'Read experiments', 'experiments', 'read', auth.uid()),
('experiments:update', 'Update experiments', 'experiments', 'update', auth.uid()),
('experiments:delete', 'Delete experiments', 'experiments', 'delete', auth.uid())
ON CONFLICT (resource_type, action) DO NOTHING;

-- Prediction permissions
INSERT INTO permissions (name, description, resource_type, action, created_by)
VALUES 
('predictions:create', 'Create predictions', 'predictions', 'create', auth.uid()),
('predictions:read', 'Read predictions', 'predictions', 'read', auth.uid()),
('predictions:update', 'Update predictions', 'predictions', 'update', auth.uid()),
('predictions:delete', 'Delete predictions', 'predictions', 'delete', auth.uid())
ON CONFLICT (resource_type, action) DO NOTHING;

-- Project permissions
INSERT INTO permissions (name, description, resource_type, action, created_by)
VALUES 
('projects:create', 'Create projects', 'projects', 'create', auth.uid()),
('projects:read', 'Read projects', 'projects', 'read', auth.uid()),
('projects:update', 'Update projects', 'projects', 'update', auth.uid()),
('projects:delete', 'Delete projects', 'projects', 'delete', auth.uid())
ON CONFLICT (resource_type, action) DO NOTHING;

-- Team permissions
INSERT INTO permissions (name, description, resource_type, action, created_by)
VALUES 
('teams:create', 'Create teams', 'teams', 'create', auth.uid()),
('teams:read', 'Read teams', 'teams', 'read', auth.uid()),
('teams:update', 'Update teams', 'teams', 'update', auth.uid()),
('teams:delete', 'Delete teams', 'teams', 'delete', auth.uid()),
('teams:manage_members', 'Manage team members', 'teams', 'manage_members', auth.uid())
ON CONFLICT (resource_type, action) DO NOTHING;

-- Assign permissions to roles
-- Admin role gets all permissions
INSERT INTO role_permissions (role_id, permission_id, created_by)
SELECT r.id, p.id, auth.uid()
FROM roles r, permissions p
WHERE r.name = 'admin'
ON CONFLICT (role_id, permission_id) DO NOTHING;

-- User role gets basic permissions
INSERT INTO role_permissions (role_id, permission_id, created_by)
SELECT r.id, p.id, auth.uid()
FROM roles r, permissions p
WHERE r.name = 'user'
AND p.name IN (
    'molecules:read', 'molecules:create', 
    'mixtures:read', 'mixtures:create',
    'experiments:read', 'experiments:create',
    'predictions:read', 'predictions:create',
    'projects:read', 'projects:create',
    'teams:read'
)
ON CONFLICT (role_id, permission_id) DO NOTHING;

-- Curator role gets extended permissions
INSERT INTO role_permissions (role_id, permission_id, created_by)
SELECT r.id, p.id, auth.uid()
FROM roles r, permissions p
WHERE r.name = 'curator'
AND p.name IN (
    'molecules:read', 'molecules:create', 'molecules:update',
    'mixtures:read', 'mixtures:create', 'mixtures:update',
    'experiments:read', 'experiments:create', 'experiments:update',
    'predictions:read', 'predictions:create', 'predictions:update',
    'projects:read', 'projects:create', 'projects:update',
    'teams:read', 'teams:create'
)
ON CONFLICT (role_id, permission_id) DO NOTHING;

-- Viewer role gets read-only permissions
INSERT INTO role_permissions (role_id, permission_id, created_by)
SELECT r.id, p.id, auth.uid()
FROM roles r, permissions p
WHERE r.name = 'viewer'
AND p.name IN (
    'molecules:read',
    'mixtures:read',
    'experiments:read',
    'predictions:read',
    'projects:read',
    'teams:read'
)
ON CONFLICT (role_id, permission_id) DO NOTHING;

-- Make the first user an admin (if they don't have a role yet)
INSERT INTO user_roles (user_id, role_id, created_by)
SELECT 
    u.id, 
    r.id,
    auth.uid()
FROM 
    auth.users u,
    roles r
WHERE 
    r.name = 'admin'
    AND u.id = (SELECT id FROM auth.users ORDER BY created_at ASC LIMIT 1)
    AND NOT EXISTS (
        SELECT 1 FROM user_roles ur WHERE ur.user_id = u.id
    )
ON CONFLICT (user_id, role_id) DO NOTHING;

-- Update RLS policies on existing tables to use role-based permissions
-- This is just an example for the molecules table, similar policies would be created for other tables
CREATE POLICY "Molecules are viewable by users with read permission"
ON molecules FOR SELECT
USING (has_permission(auth.uid(), 'molecules', 'read'));

CREATE POLICY "Molecules can be inserted by users with create permission"
ON molecules FOR INSERT
WITH CHECK (has_permission(auth.uid(), 'molecules', 'create'));

CREATE POLICY "Molecules can be updated by users with update permission"
ON molecules FOR UPDATE
USING (has_permission(auth.uid(), 'molecules', 'update') OR auth.uid() = created_by);

CREATE POLICY "Molecules can be deleted by users with delete permission"
ON molecules FOR DELETE
USING (has_permission(auth.uid(), 'molecules', 'delete') OR auth.uid() = created_by);