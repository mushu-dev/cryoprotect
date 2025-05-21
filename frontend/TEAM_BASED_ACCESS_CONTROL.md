# Team-Based Access Control

This document describes the implementation of comprehensive team-based access control in the CryoProtect database using Row Level Security (RLS).

## Overview

The team-based access control system enables:

1. Organization of users into teams with different permission levels
2. Secure sharing of molecules and mixtures within teams
3. Fine-grained control over who can view and modify resources
4. Invitation-based team membership management
5. Resource ownership transfer between users and teams

## Database Schema

### Core Tables

#### `teams`
```sql
CREATE TABLE public.teams (
    id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
    name TEXT NOT NULL,
    description TEXT,
    created_by UUID NOT NULL REFERENCES auth.users(id),
    created_at TIMESTAMPTZ NOT NULL DEFAULT NOW(),
    updated_at TIMESTAMPTZ NOT NULL DEFAULT NOW()
);
```

#### `team_members`
```sql
CREATE TABLE public.team_members (
    team_id UUID NOT NULL REFERENCES public.teams(id) ON DELETE CASCADE,
    user_id UUID NOT NULL REFERENCES auth.users(id) ON DELETE CASCADE,
    role TEXT NOT NULL DEFAULT 'read',
    joined_at TIMESTAMPTZ NOT NULL DEFAULT NOW(),
    PRIMARY KEY (team_id, user_id),
    CONSTRAINT valid_role CHECK (role IN ('admin', 'write', 'read'))
);
```

#### `team_invitations`
```sql
CREATE TABLE public.team_invitations (
    id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
    team_id UUID NOT NULL REFERENCES teams(id) ON DELETE CASCADE,
    email TEXT NOT NULL,
    role TEXT NOT NULL DEFAULT 'read',
    invited_by UUID NOT NULL REFERENCES auth.users(id),
    created_at TIMESTAMPTZ NOT NULL DEFAULT NOW(),
    expires_at TIMESTAMPTZ NOT NULL DEFAULT NOW() + INTERVAL '7 days',
    token TEXT NOT NULL UNIQUE DEFAULT encode(gen_random_bytes(32), 'hex'),
    CONSTRAINT valid_role CHECK (role IN ('admin', 'write', 'read'))
);
```

### Resource Connection

The `molecules` and `mixtures` tables have a `team_id` column that links resources to teams.

## Team Roles and Permissions

The system implements three role levels:

1. **admin** - Can manage team membership and modify all team resources
2. **write** - Can create and modify resources within the team
3. **read** - Can only view team resources (default)

## Helper Functions

Several helper functions are implemented to facilitate team-based access control:

### `auth.is_team_member(team_id UUID) → BOOLEAN`
Checks if the current user is a member of the specified team.

### `auth.is_team_admin(team_id UUID) → BOOLEAN`
Checks if the current user is an admin of the specified team.

### `auth.get_user_teams() → TABLE(team_id UUID, role TEXT)`
Returns all teams the current user belongs to, along with their role in each team.

### `auth.can_access_team_molecule(molecule_id UUID) → BOOLEAN`
Checks if the current user can access a molecule through team membership.

### `auth.can_access_team_mixture(mixture_id UUID) → BOOLEAN`
Checks if the current user can access a mixture through team membership.

### `auth.can_modify_team_resource(team_id UUID) → BOOLEAN`
Checks if the current user has permission to modify a team resource.

## Team Management Functions

### `auth.create_team(team_name TEXT, team_description TEXT DEFAULT NULL) → UUID`
Creates a new team and adds the current user as an admin.

### `auth.invite_to_team(team_id UUID, email TEXT, role TEXT DEFAULT 'read') → UUID`
Invites a user to a team by email.

### `auth.accept_team_invitation(invitation_token TEXT) → BOOLEAN`
Accepts a team invitation using the invitation token.

### `auth.transfer_resource_to_team(resource_type TEXT, resource_id UUID, team_id UUID) → BOOLEAN`
Transfers ownership of a resource to a team.

## Row Level Security (RLS) Policies

### Teams Table
- **teams_access_policy**: Allows users to see teams they belong to
- **teams_modify_policy**: Allows only team admins and system admins to modify teams

### Team Members Table
- **team_members_access_policy**: Allows members to see who else is in their teams
- **team_members_modify_policy**: Allows only team admins to manage team membership

### Molecules Table
- **molecules_team_access_policy**: Allows team members to access molecules
- **molecules_team_modify_policy**: Allows team members with write or admin roles to modify molecules

### Mixtures Table
- **mixtures_team_access_policy**: Allows team members to access mixtures
- **mixtures_team_modify_policy**: Allows team members with write or admin roles to modify mixtures

### Team Invitations Table
- **team_invitations_access_policy**: Controls who can see invitations
- **team_invitations_modify_policy**: Controls who can create and modify invitations

## Views

### `auth.team_resources`
A view that shows all resources (molecules, mixtures) that belong to teams.

## Triggers

### `check_team_resource_permission`
Ensures that users can only assign resources to teams they have write permissions for.

## Usage Examples

### Creating a New Team
```sql
SELECT auth.create_team('Research Team Alpha', 'Our primary research team');
```

### Inviting a User to a Team
```sql
SELECT auth.invite_to_team('team_id_here', 'colleague@example.com', 'write');
```

### Accepting a Team Invitation
```sql
SELECT auth.accept_team_invitation('invitation_token_here');
```

### Transferring a Resource to a Team
```sql
SELECT auth.transfer_resource_to_team('molecule', 'molecule_id_here', 'team_id_here');
```

### Viewing Team Resources
```sql
SELECT * FROM auth.team_resources WHERE team_id = 'team_id_here';
```

## Testing and Verification

The implementation includes comprehensive tests to verify:

1. The existence and structure of team-related tables
2. The correct implementation of team helper functions
3. The proper application of RLS policies
4. The functionality of the team resources view

To run these tests:
```bash
./run_team_access_implementation.sh
```

## Security Considerations

1. All security-critical functions use `SECURITY DEFINER` to ensure consistent privilege evaluation
2. Resource ownership checks are performed before team permission checks
3. RLS policies are applied consistently across all team-related tables
4. Team invitations expire after 7 days for security

## Frontend Integration

The frontend should integrate with these features through:

1. Team management UI for creating and managing teams
2. Team invitation and acceptance workflow
3. Team-aware resource browsers
4. Permission indicators for team-owned resources

## Next Steps

1. Implement team analytics and usage statistics
2. Add team-based quota management
3. Develop team resource sharing features
4. Create team activity logs for auditing