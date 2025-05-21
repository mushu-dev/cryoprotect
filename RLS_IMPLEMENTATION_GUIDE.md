# RLS Implementation Guide

This document explains the Row Level Security (RLS) implementation in the CryoProtect database. It covers the team-based access model, security definer functions, and policies for each table.

## Team-Based Access Model

The CryoProtect database uses a team-based access model with the following structure:

1. **Teams**: Core organizational units
2. **Team Members**: Users belonging to teams with roles (admin or member)
3. **Projects**: Belong to teams
4. **Experiments**: Belong to projects
5. **Molecules**: Associated with projects through experiments
6. **Mixtures**: Belong to projects

This structure allows for multi-level access control where access to data is determined by team membership.

## Security Definer Functions

The RLS implementation relies on several security definer functions that efficiently implement the access control logic:

### Core Access Functions

- `auth.is_service_role()`: Checks if the current user has the service role
  ```sql
  SELECT (current_setting('request.jwt.claims', true)::jsonb->'role')::text = '"service_role"'
  ```

- `auth.is_team_member(team_id)`: Checks if a user is a member of a team
  ```sql
  SELECT EXISTS (
    SELECT 1 FROM team_members 
    WHERE team_id = $1 
    AND user_id = auth.uid()
  )
  ```

- `auth.is_team_admin(team_id)`: Checks if a user is an admin of a team
  ```sql
  SELECT EXISTS (
    SELECT 1 FROM team_members 
    WHERE team_id = $1 
    AND user_id = auth.uid()
    AND role = 'admin'
  )
  ```

### Project Access Functions

- `auth.can_access_project(project_id)`: Checks if a user can access a project via team membership
  ```sql
  -- Get the team_id for the project then check membership
  SELECT auth.is_team_member((SELECT team_id FROM projects WHERE id = $1))
  ```

- `auth.can_manage_project(project_id)`: Checks if a user can manage a project via team admin role
  ```sql
  -- Get the team_id for the project then check admin role
  SELECT auth.is_team_admin((SELECT team_id FROM projects WHERE id = $1))
  ```

### Data Access Functions

- `auth.can_access_molecule(molecule_id)`: Checks if a user can access a molecule via experiments and projects
  ```sql
  SELECT EXISTS (
    SELECT 1
    FROM experiments e
    JOIN projects p ON e.project_id = p.id
    JOIN teams t ON p.team_id = t.id
    JOIN team_members tm ON t.id = tm.team_id
    WHERE e.molecule_id = molecule_id
    AND tm.user_id = auth.uid()
  )
  ```

- `auth.get_user_teams()`: Returns all teams a user belongs to
  ```sql
  SELECT team_id FROM team_members WHERE user_id = auth.uid()
  ```

- `auth.get_accessible_projects()`: Returns all projects a user can access
  ```sql
  SELECT p.id
  FROM projects p
  JOIN teams t ON p.team_id = t.id
  JOIN team_members tm ON t.id = tm.team_id
  WHERE tm.user_id = auth.uid()
  ```

## Table Policies

Each table in the database has specific RLS policies that govern access:

### Projects Table

- **SELECT**: User can view projects from teams they are members of
- **INSERT**: User can create projects in teams they are members of
- **UPDATE**: Only team admins can update projects
- **DELETE**: Only team admins can delete projects

### Experiments Table

- **SELECT/INSERT/UPDATE/DELETE**: User can access experiments if they can access the associated project

### Molecules Table

- **SELECT**: User can view molecules used in experiments for projects they can access
- **INSERT**: Allow insertion (association happens in experiments)
- **UPDATE/DELETE**: User can modify/delete molecules they can access via experiments

### Molecular Properties Table

- **SELECT/INSERT/UPDATE/DELETE**: User can access molecular properties for molecules they can access

### Mixtures Table

- **SELECT/INSERT/UPDATE/DELETE**: User can access mixtures for projects they can access

### Team Members Table

- **SELECT**: User can view their own team memberships and members of teams they administer
- **INSERT/UPDATE/DELETE**: Only team admins can manage team memberships

## Performance Considerations

The implementation includes performance indexes to optimize RLS policy evaluation:

- `idx_team_members_user_id` on `team_members (user_id)`
- `idx_team_members_team_id` on `team_members (team_id)`
- `idx_team_members_role` on `team_members (role)`
- `idx_projects_team_id` on `projects (team_id)`
- `idx_experiments_project_id` on `experiments (project_id)`
- `idx_experiments_molecule_id` on `experiments (molecule_id)`
- `idx_molecular_properties_molecule_id` on `molecular_properties (molecule_id)`
- `idx_mixtures_project_id` on `mixtures (project_id)`

## Service Role Access

The service role has full access to all tables through dedicated policies:

```sql
CREATE POLICY [table]_service_role_policy ON [table]
    FOR ALL USING (auth.is_service_role());
```

## Troubleshooting

Common issues and their solutions:

1. **Access Denied Errors**: If users cannot access data they should be able to, check:
   - Team membership in `team_members` table
   - Project team assignment in `projects` table
   - Experiment-project association in `experiments` table

2. **Performance Issues**: If queries are slow, check:
   - Existence of necessary indexes
   - Query structure (avoid complex joins in user queries)
   - Function complexity (keep security definer functions simple)

3. **Service Role Issues**: If service role cannot access data, check:
   - Proper JWT configuration
   - Existence of service role policies for all tables
   - Implementation of `auth.is_service_role()` function

## Testing RLS Policies

Use the `run_rls_policy_improvements_test.sh` script to test the RLS policies:

```bash
./run_rls_policy_improvements_test.sh --report
```

This will generate a report in the `reports` directory with detailed test results.