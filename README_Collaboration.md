# CryoProtect Analyzer - Collaboration Features

This document provides an overview of the collaboration features in the CryoProtect Analyzer application, including team management, resource sharing, commenting, and notifications.

## Table of Contents

1. [Overview](#overview)
2. [Team Management](#team-management)
3. [Resource Sharing](#resource-sharing)
4. [Role-Based Permissions](#role-based-permissions)
5. [Comments and Discussions](#comments-and-discussions)
6. [Notifications](#notifications)
7. [Activity Tracking](#activity-tracking)
8. [API Endpoints](#api-endpoints)
9. [Database Schema](#database-schema)

## Overview

The collaboration features in CryoProtect Analyzer enable research teams to work together effectively by:

- Creating and managing teams
- Sharing projects, experiments, and mixtures with team members
- Implementing role-based permissions for access control
- Adding comments and discussions on shared resources
- Providing notifications for team activities
- Tracking team activity for transparency and accountability

## Team Management

### Creating Teams

1. Navigate to the Teams page from the main navigation menu
2. Click the "Create New Team" button
3. Enter a name and optional description for the team
4. Submit the form to create the team

When you create a team, you are automatically added as an admin member.

### Managing Team Members

As a team admin, you can:

1. Add new members by email
2. Assign roles to members (admin, editor, viewer)
3. Update member roles
4. Remove members from the team

## Resource Sharing

You can share the following resources with teams:

- Projects
- Mixtures
- Experiments

To share a resource:

1. Navigate to the Teams page
2. Select a team
3. Go to the "Shared Resources" tab
4. Click "Share Resource"
5. Select the resource type and enter the resource ID
6. Submit the form to share the resource

Shared resources are accessible to all team members according to their roles.

## Role-Based Permissions

The system implements three roles with different permission levels:

### Admin
- Create, edit, and delete teams
- Add, update, and remove team members
- Share and unshare resources
- Edit and delete shared resources
- Add comments and discussions

### Editor
- View team information and members
- View shared resources
- Edit shared resources
- Add comments and discussions

### Viewer
- View team information and members
- View shared resources
- Add comments and discussions

## Comments and Discussions

Users can add comments to shared resources to facilitate discussions:

1. Navigate to a shared resource
2. Scroll to the comments section
3. Enter your comment in the text field
4. Submit the form to add the comment

Comments support:
- Text formatting
- Replies to other comments
- Editing and deleting your own comments

## Notifications

The system provides notifications for various team activities:

- When you are added to a team
- When your role is updated
- When a resource is shared with your team
- When someone replies to your comment
- When a team member makes changes to a shared resource

To view notifications:
1. Click the Notifications icon in the navigation bar
2. View the list of notifications
3. Mark notifications as read or delete them

## Activity Tracking

The system tracks team activities for transparency and accountability:

- Team creation and updates
- Member additions, role updates, and removals
- Resource sharing and unsharing
- Comments and discussions

To view team activity:
1. Navigate to the Teams page
2. Select a team
3. Go to the "Activity" tab
4. View the chronological list of team activities

## API Endpoints

The collaboration features are accessible through the following API endpoints:

### Teams
- `GET /api/v1/teams` - Get all teams for the current user
- `POST /api/v1/teams` - Create a new team
- `GET /api/v1/teams/:id` - Get a team by ID
- `PUT /api/v1/teams/:id` - Update a team
- `DELETE /api/v1/teams/:id` - Delete a team

### Team Members
- `GET /api/v1/teams/:id/members` - Get all members of a team
- `POST /api/v1/teams/:id/members` - Add a member to a team
- `PUT /api/v1/teams/:id/members/:userId` - Update a team member's role
- `DELETE /api/v1/teams/:id/members/:userId` - Remove a member from a team

### Shared Resources
- `GET /api/v1/teams/:id/resources` - Get all resources shared with a team
- `POST /api/v1/teams/:id/resources` - Share a resource with a team
- `DELETE /api/v1/teams/:id/resources/:type/:resourceId` - Unshare a resource from a team

### Comments
- `GET /api/v1/comments/:resourceType/:resourceId` - Get all comments for a resource
- `POST /api/v1/comments/:resourceType/:resourceId` - Add a comment to a resource
- `PUT /api/v1/comments/:id` - Update a comment
- `DELETE /api/v1/comments/:id` - Delete a comment

### Notifications
- `GET /api/v1/notifications` - Get notifications for the current user
- `PUT /api/v1/notifications/:id` - Mark a notification as read
- `PUT /api/v1/notifications` - Mark all notifications as read
- `DELETE /api/v1/notifications/:id` - Delete a notification

### Activity
- `GET /api/v1/teams/:id/activity` - Get activity log for a team

## Database Schema

The collaboration features use the following database tables:

### teams
- `id` - UUID primary key
- `name` - Team name
- `description` - Team description
- `created_at` - Creation timestamp
- `updated_at` - Update timestamp
- `created_by` - User ID of the creator

### team_members
- `id` - UUID primary key
- `team_id` - Reference to teams.id
- `user_id` - Reference to auth.users.id
- `role` - Role of the user in the team (admin, editor, viewer)
- `created_at` - Creation timestamp
- `updated_at` - Update timestamp
- `created_by` - User ID of the creator

### shared_resources
- `id` - UUID primary key
- `team_id` - Reference to teams.id
- `resource_type` - Type of resource (project, mixture, experiment)
- `resource_id` - ID of the resource
- `created_at` - Creation timestamp
- `created_by` - User ID of the creator

### comments
- `id` - UUID primary key
- `resource_type` - Type of resource (project, mixture, experiment)
- `resource_id` - ID of the resource
- `content` - Comment content
- `parent_id` - Reference to comments.id for replies
- `created_at` - Creation timestamp
- `updated_at` - Update timestamp
- `created_by` - User ID of the creator

### notifications
- `id` - UUID primary key
- `user_id` - Reference to auth.users.id
- `team_id` - Reference to teams.id (optional)
- `title` - Notification title
- `message` - Notification message
- `resource_type` - Type of resource (team, project, mixture, experiment, comment)
- `resource_id` - ID of the resource
- `is_read` - Whether the notification has been read
- `created_at` - Creation timestamp
- `created_by` - User ID of the creator

### activity_log
- `id` - UUID primary key
- `team_id` - Reference to teams.id
- `user_id` - Reference to auth.users.id
- `action` - Action performed
- `resource_type` - Type of resource
- `resource_id` - ID of the resource
- `details` - Additional details as JSON
- `created_at` - Creation timestamp