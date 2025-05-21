# CryoProtect v2 - Role-Based Access Control (RBAC)

This document provides comprehensive information on the role-based access control (RBAC) system implemented in CryoProtect v2.

## Overview

The RBAC system provides a flexible and secure way to control access to resources and operations within the CryoProtect v2 application. It is integrated with the existing JWT authentication, session management, and token management systems.

Key features:
- Hierarchical role system
- Fine-grained permission control
- Database-level enforcement through Row Level Security (RLS)
- Application-level enforcement through decorators
- API endpoints for role and permission management

## Role Hierarchy

The RBAC system defines a hierarchical role structure:

| Role | Description | Inherits From |
|------|-------------|---------------|
| `admin` | Administrator with full system access | All roles |
| `curator` | User with special privileges for curating data | `user`, `viewer` |
| `user` | Regular user with standard access | `viewer` |
| `viewer` | User with read-only access to data | None |

This hierarchy means that higher-level roles inherit all permissions from lower-level roles. For example, an `admin` has all the permissions of a `curator`, `user`, and `viewer`.

## Permissions Model

Permissions in the system are defined as a combination of a resource type and an action:

```
permission = resource_type:action
```

### Resource Types

The system defines the following resource types:

- `molecules`: Molecular data
- `mixtures`: Mixture data
- `experiments`: Experiment data
- `predictions`: Prediction data
- `projects`: Project data
- `teams`: Team data
- `users`: User management
- `roles`: Role management
- `permissions`: Permission management
- `system`: System configuration
- `admin`: Administrative functions

### Actions

The system defines the following actions:

- `create`: Create new resources
- `read`: View resources
- `update`: Modify existing resources
- `delete`: Remove resources
- `manage`: Special management operations
- `all`: All possible actions

### Default Permissions by Role

#### Admin
- All permissions on all resource types

#### Curator
- `molecules`: create, read, update
- `mixtures`: create, read, update
- `experiments`: create, read, update
- `predictions`: create, read, update
- `projects`: create, read, update
- `teams`: create, read

#### User
- `molecules`: create, read
- `mixtures`: create, read
- `experiments`: create, read
- `predictions`: create, read
- `projects`: create, read
- `teams`: read

#### Viewer
- `molecules`: read
- `mixtures`: read
- `experiments`: read
- `predictions`: read
- `projects`: read
- `teams`: read

## Database Schema

The RBAC system uses the following database tables:

### roles

Stores the available roles in the system.

```sql
CREATE TABLE roles (
    id UUID PRIMARY KEY DEFAULT uuid_generate_v4(),
    name VARCHAR(50) NOT NULL UNIQUE,
    description TEXT,
    is_system_role BOOLEAN DEFAULT FALSE,
    created_at TIMESTAMP WITH TIME ZONE DEFAULT NOW(),
    updated_at TIMESTAMP WITH TIME ZONE DEFAULT NOW(),
    created_by UUID REFERENCES auth.users(id) ON DELETE SET NULL
);
```

### permissions

Stores the available permissions in the system.

```sql
CREATE TABLE permissions (
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
```

### role_permissions

Maps roles to permissions.

```sql
CREATE TABLE role_permissions (
    id UUID PRIMARY KEY DEFAULT uuid_generate_v4(),
    role_id UUID NOT NULL REFERENCES roles(id) ON DELETE CASCADE,
    permission_id UUID NOT NULL REFERENCES permissions(id) ON DELETE CASCADE,
    created_at TIMESTAMP WITH TIME ZONE DEFAULT NOW(),
    created_by UUID REFERENCES auth.users(id) ON DELETE SET NULL,
    UNIQUE(role_id, permission_id)
);
```

### user_roles

Maps users to roles.

```sql
CREATE TABLE user_roles (
    id UUID PRIMARY KEY DEFAULT uuid_generate_v4(),
    user_id UUID NOT NULL REFERENCES auth.users(id) ON DELETE CASCADE,
    role_id UUID NOT NULL REFERENCES roles(id) ON DELETE CASCADE,
    created_at TIMESTAMP WITH TIME ZONE DEFAULT NOW(),
    updated_at TIMESTAMP WITH TIME ZONE DEFAULT NOW(),
    created_by UUID REFERENCES auth.users(id) ON DELETE SET NULL,
    UNIQUE(user_id, role_id)
);
```

## Database Functions

The RBAC system provides the following database functions:

### has_role(p_user_id UUID, p_role_name VARCHAR)

Checks if a user has a specific role.

```sql
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
```

### has_permission(p_user_id UUID, p_resource_type VARCHAR, p_action VARCHAR)

Checks if a user has a specific permission.

```sql
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
```

### get_user_roles(p_user_id UUID)

Gets all roles for a user.

```sql
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
```

### get_user_permissions(p_user_id UUID)

Gets all permissions for a user.

```sql
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
```

## Row Level Security (RLS) Policies

The RBAC system enforces access control at the database level through RLS policies. Here's an example of how RLS policies are defined for the `molecules` table:

```sql
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
```

Similar policies are defined for other tables in the system.

## Python API

The RBAC system provides a Python API for role and permission management.

### Role Management

```python
from api.rbac import RoleManager

# Get all roles
roles = RoleManager.get_all_roles()

# Get a role by ID
role = RoleManager.get_role(role_id)

# Get a role by name
role = RoleManager.get_role_by_name(role_name)

# Create a new role
role = RoleManager.create_role(name, description, is_system_role)

# Update a role
role = RoleManager.update_role(role_id, data)

# Delete a role
success = RoleManager.delete_role(role_id)

# Get role permissions
permissions = RoleManager.get_role_permissions(role_id)

# Assign permission to role
role_permission = RoleManager.assign_permission_to_role(role_id, permission_id)

# Remove permission from role
success = RoleManager.remove_permission_from_role(role_id, permission_id)
```

### Permission Management

```python
from api.rbac import PermissionManager

# Get all permissions
permissions = PermissionManager.get_all_permissions()

# Get a permission by ID
permission = PermissionManager.get_permission(permission_id)

# Get a permission by resource type and action
permission = PermissionManager.get_permission_by_resource_action(resource_type, action)

# Create a new permission
permission = PermissionManager.create_permission(name, resource_type, action, description)

# Update a permission
permission = PermissionManager.update_permission(permission_id, data)

# Delete a permission
success = PermissionManager.delete_permission(permission_id)
```

### User Role Management

```python
from api.rbac import UserRoleManager

# Get user roles
roles = UserRoleManager.get_user_roles(user_id)

# Get user permissions
permissions = UserRoleManager.get_user_permissions(user_id)

# Assign role to user
user_role = UserRoleManager.assign_role_to_user(user_id, role_id)

# Remove role from user
success = UserRoleManager.remove_role_from_user(user_id, role_id)

# Check if user has role
has_role = UserRoleManager.has_role(user_id, role_name)

# Check if user has permission
has_permission = UserRoleManager.has_permission(user_id, resource_type, action)
```

## Decorators

The RBAC system provides decorators for protecting routes based on roles and permissions.

### role_required

Requires a specific role for API endpoints.

```python
from api.rbac import role_required

@app.route('/admin/dashboard')
@role_required('admin')
def admin_dashboard():
    # Only accessible to users with the 'admin' role
    return jsonify({'message': 'Admin dashboard'})

@app.route('/curator/edit')
@role_required(['curator', 'admin'])
def curator_edit():
    # Accessible to users with either 'curator' or 'admin' role
    return jsonify({'message': 'Curator edit page'})
```

### permission_required

Requires a specific permission for API endpoints.

```python
from api.rbac import permission_required

@app.route('/molecules')
@permission_required('molecules', 'read')
def get_molecules():
    # Only accessible to users with the 'molecules:read' permission
    return jsonify({'message': 'Molecules list'})

@app.route('/molecules', methods=['POST'])
@permission_required('molecules', 'create')
def create_molecule():
    # Only accessible to users with the 'molecules:create' permission
    return jsonify({'message': 'Molecule created'})
```

### admin_required

Shorthand for requiring the admin role.

```python
from api.rbac import admin_required

@app.route('/admin/settings')
@admin_required
def admin_settings():
    # Only accessible to users with the 'admin' role
    return jsonify({'message': 'Admin settings'})
```

## API Endpoints

The RBAC system provides API endpoints for managing roles and permissions.

### Role Management

- `GET /api/v1/rbac/roles`: Get all roles
- `GET /api/v1/rbac/roles/{role_id}`: Get a role by ID
- `POST /api/v1/rbac/roles`: Create a new role
- `PUT /api/v1/rbac/roles/{role_id}`: Update a role
- `DELETE /api/v1/rbac/roles/{role_id}`: Delete a role
- `GET /api/v1/rbac/roles/{role_id}/permissions`: Get all permissions for a role
- `POST /api/v1/rbac/roles/{role_id}/permissions/{permission_id}`: Assign a permission to a role
- `DELETE /api/v1/rbac/roles/{role_id}/permissions/{permission_id}`: Remove a permission from a role

### Permission Management

- `GET /api/v1/rbac/permissions`: Get all permissions
- `GET /api/v1/rbac/permissions/{permission_id}`: Get a permission by ID
- `POST /api/v1/rbac/permissions`: Create a new permission
- `PUT /api/v1/rbac/permissions/{permission_id}`: Update a permission
- `DELETE /api/v1/rbac/permissions/{permission_id}`: Delete a permission

### User Role Management

- `GET /api/v1/rbac/users/{user_id}/roles`: Get all roles for a user
- `GET /api/v1/rbac/users/{user_id}/permissions`: Get all permissions for a user
- `POST /api/v1/rbac/users/{user_id}/roles/{role_id}`: Assign a role to a user
- `DELETE /api/v1/rbac/users/{user_id}/roles/{role_id}`: Remove a role from a user

### Current User

- `GET /api/v1/rbac/me/roles`: Get all roles for the current user
- `GET /api/v1/rbac/me/permissions`: Get all permissions for the current user
- `GET /api/v1/rbac/me/has-role/{role_name}`: Check if the current user has a specific role
- `GET /api/v1/rbac/me/has-permission/{resource_type}/{action}`: Check if the current user has a specific permission

## Integration with Authentication System

The RBAC system is integrated with the existing JWT authentication system. When a user logs in, their roles and permissions are loaded and included in the JWT token. This allows for efficient access control without the need for additional database queries on every request.

### JWT Token Claims

The JWT token includes the following claims related to RBAC:

- `role`: The user's primary role (for backward compatibility)
- `roles`: Array of all roles assigned to the user
- `permissions`: Array of all permissions granted to the user

### Session Management

The RBAC system is also integrated with the session management system. When a user's roles or permissions are changed, their active sessions are updated to reflect the changes.

## Configuration

RBAC configuration is defined in `auth_config.py`:

```python
# RBAC Configuration
RBAC_ENABLED = True
ROLE_HIERARCHY = {
    "admin": ["curator", "user", "viewer"],  # Admin inherits all permissions
    "curator": ["user", "viewer"],           # Curator inherits user and viewer permissions
    "user": ["viewer"],                      # User inherits viewer permissions
    "viewer": []                             # Viewer has no inheritance
}

# Resource types for permissions
RESOURCE_TYPES = [
    "molecules", "mixtures", "experiments", "predictions", 
    "projects", "teams", "users", "roles", "permissions", 
    "system", "admin"
]

# Actions for permissions
ACTIONS = [
    "create", "read", "update", "delete", "manage", "all"
]
```

## Implementation

The RBAC system is implemented in the following files:

- `migrations/011_rbac_schema.sql`: Database schema for RBAC
- `api/rbac.py`: Core RBAC functionality
- `api/rbac_routes.py`: API endpoints for RBAC
- `auth_config.py`: RBAC configuration
- `apply_rbac_schema.py`: Script to apply the RBAC schema migration

## Applying the RBAC Schema

To apply the RBAC schema to your database, run the following command:

```bash
# On Windows
apply_rbac_schema.bat

# On Unix-based systems
./apply_rbac_schema.sh
```

This will create the necessary tables, functions, and policies for the RBAC system.

## Security Considerations

- The RBAC system enforces access control at both the application and database levels.
- RLS policies ensure that users can only access data they are authorized to see.
- Decorators ensure that users can only perform actions they are authorized to perform.
- Role and permission management is restricted to administrators.
- System roles cannot be deleted.
- Users must have at least one role.

## Best Practices

- Always use the `permission_required` decorator for API endpoints that access sensitive data.
- Use the `role_required` decorator for administrative endpoints.
- When creating new tables, define RLS policies that use the `has_permission` function.
- When adding new functionality, define appropriate permissions and assign them to roles.
- Regularly audit user roles and permissions.