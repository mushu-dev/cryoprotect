"""
CryoProtect Analyzer API - Team Models

This module contains data models for team collaboration features.
"""

from typing import Dict, List, Any
from flask_restful import fields
from marshmallow import Schema, fields as ma_fields, validate

from api.models import BaseModel
from api.utils import get_user_id

# Team-related fields
team_fields = {
    'id': fields.String,
    'name': fields.String,
    'description': fields.String,
    'created_at': fields.DateTime(dt_format='iso8601'),
    'updated_at': fields.DateTime(dt_format='iso8601'),
    'created_by': fields.String,
    'member_count': fields.Integer
}

team_member_fields = {
    'id': fields.String,
    'team_id': fields.String,
    'user_id': fields.String,
    'role': fields.String,
    'created_at': fields.DateTime(dt_format='iso8601'),
    'email': fields.String,
    'display_name': fields.String
}

shared_resource_fields = {
    'id': fields.String,
    'team_id': fields.String,
    'resource_type': fields.String,
    'resource_id': fields.String,
    'created_at': fields.DateTime(dt_format='iso8601'),
    'created_by': fields.String
}

comment_fields = {
    'id': fields.String,
    'resource_type': fields.String,
    'resource_id': fields.String,
    'content': fields.String,
    'parent_id': fields.String,
    'created_at': fields.DateTime(dt_format='iso8601'),
    'updated_at': fields.DateTime(dt_format='iso8601'),
    'created_by': fields.String,
    'user_email': fields.String,
    'user_display_name': fields.String
}

notification_fields = {
    'id': fields.String,
    'user_id': fields.String,
    'team_id': fields.String,
    'title': fields.String,
    'message': fields.String,
    'resource_type': fields.String,
    'resource_id': fields.String,
    'is_read': fields.Boolean,
    'created_at': fields.DateTime(dt_format='iso8601'),
    'created_by': fields.String
}

activity_log_fields = {
    'id': fields.String,
    'team_id': fields.String,
    'user_id': fields.String,
    'user_email': fields.String,
    'action': fields.String,
    'resource_type': fields.String,
    'resource_id': fields.String,
    'details': fields.Raw,
    'created_at': fields.DateTime(dt_format='iso8601')
}

# Team-related schemas
class TeamSchema(Schema):
    """Schema for creating or updating a team."""
    name = ma_fields.String(required=True)
    description = ma_fields.String()

class TeamMemberSchema(Schema):
    """Schema for adding or updating a team member."""
    user_id = ma_fields.UUID(required=True)
    role = ma_fields.String(required=True, validate=validate.OneOf(['admin', 'editor', 'viewer']))

class SharedResourceSchema(Schema):
    """Schema for sharing a resource with a team."""
    team_id = ma_fields.UUID(required=True)
    resource_type = ma_fields.String(required=True, validate=validate.OneOf(['project', 'mixture', 'experiment']))
    resource_id = ma_fields.UUID(required=True)

class CommentSchema(Schema):
    """Schema for creating or updating a comment."""
    resource_type = ma_fields.String(required=True, validate=validate.OneOf(['project', 'mixture', 'experiment']))
    resource_id = ma_fields.UUID(required=True)
    content = ma_fields.String(required=True)
    parent_id = ma_fields.UUID(allow_none=True)

class NotificationSchema(Schema):
    """Schema for creating a notification."""
    user_id = ma_fields.UUID(required=True)
    team_id = ma_fields.UUID()
    title = ma_fields.String(required=True)
    message = ma_fields.String(required=True)
    resource_type = ma_fields.String(validate=validate.OneOf(['team', 'project', 'mixture', 'experiment', 'comment']))
    resource_id = ma_fields.UUID()

# Team model
class Team(BaseModel):
    """Model for teams in the database."""
    table_name = 'teams'

    @classmethod
    def get_with_member_count(cls, id_value: str) -> Dict[str, Any]:
        """
        Get a team with member count.
        Args:
            id_value: The ID of the team
        Returns:
            Team with member count
        """
        response = cls.get_supabase().rpc(
            "get_team_members_with_details",
            {"p_team_id": id_value}
        ).execute()

        if response.error:
            raise Exception(f"Error getting team members: {response.error}")

        team = cls.get(id_value)
        if team:
            team['member_count'] = len(response.data)

        return team

    @classmethod
    def get_user_teams(cls, user_id: str) -> List[Dict[str, Any]]:
        """
        Get all teams for a user.
        Args:
            user_id: The ID of the user
        Returns:
            List of teams with role and member count
        """
        response = cls.get_supabase().rpc(
            "get_user_teams",
            {"p_user_id": user_id}
        ).execute()

        if response.error:
            raise Exception(f"Error getting user teams: {response.error}")

        return response.data

    @classmethod
    def create_with_creator_as_admin(cls, name: str, description: str = None) -> Dict[str, Any]:
        """
        Create a new team with the creator as an admin.
        Args:
            name: Name of the team
            description: Optional description of the team
        Returns:
            The created team
        """
        user_id = get_user_id()
        if not user_id:
            raise Exception("User must be authenticated to create a team")

        # Create the team
        team_data = {
            "name": name,
            "description": description,
            "created_by": user_id
        }

        team = cls.create(team_data)

        # Add the creator as an admin
        if team:
            supabase = cls.get_supabase()
            response = supabase.table("team_members").insert({
                "team_id": team["id"],
                "user_id": user_id,
                "role": "admin",
                "created_by": user_id
            }).execute()

            if response.error:
                # Rollback team creation if adding admin fails
                cls.delete(team["id"])
                raise Exception(f"Error adding creator as admin: {response.error}")

            # Log the activity
            supabase.rpc(
                "log_team_activity",
                {
                    "p_team_id": team["id"],
                    "p_user_id": user_id,
                    "p_action": "created",
                    "p_resource_type": "team",
                    "p_resource_id": team["id"],
                    "p_details": {"team_name": name}
                }
            ).execute()

        return team

# TeamMember model
class TeamMember(BaseModel):
    """Model for team members in the database."""
    table_name = 'team_members'

    @classmethod
    def get_team_members(cls, team_id: str) -> List[Dict[str, Any]]:
        """
        Get all members of a team with user details.
        Args:
            team_id: The ID of the team
        Returns:
            List of team members with user details
        """
        response = cls.get_supabase().rpc(
            "get_team_members_with_details",
            {"p_team_id": team_id}
        ).execute()

        if response.error:
            raise Exception(f"Error getting team members: {response.error}")

        return response.data

    @classmethod
    def add_member(cls, team_id: str, user_id: str, role: str) -> Dict[str, Any]:
        """
        Add a member to a team.
        Args:
            team_id: The ID of the team
            user_id: The ID of the user to add
            role: The role of the user in the team ('admin', 'editor', 'viewer')
        Returns:
            The created team member record
        """
        current_user_id = get_user_id()

        # Check if the current user is an admin of the team
        supabase = cls.get_supabase()
        is_admin_response = supabase.rpc(
            "has_team_role",
            {"p_team_id": team_id, "p_user_id": current_user_id, "p_role": "admin"}
        ).execute()

        if is_admin_response.error:
            raise Exception(f"Error checking team role: {is_admin_response.error}")

        if not is_admin_response.data:
            raise Exception("Only team admins can add members")

        # Add the member
        member_data = {
            "team_id": team_id,
            "user_id": user_id,
            "role": role,
            "created_by": current_user_id
        }

        member = cls.create(member_data)

        # Log the activity
        if member:
            supabase.rpc(
                "log_team_activity",
                {
                    "p_team_id": team_id,
                    "p_user_id": current_user_id,
                    "p_action": "added_member",
                    "p_resource_type": "team",
                    "p_resource_id": team_id,
                    "p_details": {"member_id": user_id, "role": role}
                }
            ).execute()

            # Create notification for the added user
            supabase.table("notifications").insert({
                "user_id": user_id,
                "team_id": team_id,
                "title": "Added to Team",
                "message": f"You have been added to a team with the role of {role}",
                "resource_type": "team",
                "resource_id": team_id,
                "created_by": current_user_id
            }).execute()

        return member

    @classmethod
    def update_role(cls, team_id: str, user_id: str, new_role: str) -> Dict[str, Any]:
        """
        Update a team member's role.
        Args:
            team_id: The ID of the team
            user_id: The ID of the user
            new_role: The new role for the user
        Returns:
            The updated team member record
        """
        current_user_id = get_user_id()

        # Check if the current user is an admin of the team
        supabase = cls.get_supabase()
        is_admin_response = supabase.rpc(
            "has_team_role",
            {"p_team_id": team_id, "p_user_id": current_user_id, "p_role": "admin"}
        ).execute()

        if is_admin_response.error:
            raise Exception(f"Error checking team role: {is_admin_response.error}")

        if not is_admin_response.data:
            raise Exception("Only team admins can update member roles")

        # Get the member record
        response = supabase.table("team_members").select("id").eq("team_id", team_id).eq("user_id", user_id).execute()

        if response.error:
            raise Exception(f"Error getting team member: {response.error}")

        if not response.data:
            raise Exception("User is not a member of this team")

        member_id = response.data[0]["id"]

        # Update the role
        updated_member = cls.update(member_id, {"role": new_role})

        # Log the activity
        if updated_member:
            supabase.rpc(
                "log_team_activity",
                {
                    "p_team_id": team_id,
                    "p_user_id": current_user_id,
                    "p_action": "updated_member_role",
                    "p_resource_type": "team",
                    "p_resource_id": team_id,
                    "p_details": {"member_id": user_id, "new_role": new_role}
                }
            ).execute()

            # Create notification for the user
            supabase.table("notifications").insert({
                "user_id": user_id,
                "team_id": team_id,
                "title": "Role Updated",
                "message": f"Your role in the team has been updated to {new_role}",
                "resource_type": "team",
                "resource_id": team_id,
                "created_by": current_user_id
            }).execute()

        return updated_member

    @classmethod
    def remove_member(cls, team_id: str, user_id: str) -> bool:
        """
        Remove a member from a team.
        Args:
            team_id: The ID of the team
            user_id: The ID of the user to remove
        Returns:
            True if successful, False otherwise
        """
        current_user_id = get_user_id()

        # Check if the current user is an admin of the team
        supabase = cls.get_supabase()
        is_admin_response = supabase.rpc(
            "has_team_role",
            {"p_team_id": team_id, "p_user_id": current_user_id, "p_role": "admin"}
        ).execute()

        if is_admin_response.error:
            raise Exception(f"Error checking team role: {is_admin_response.error}")

        if not is_admin_response.data and current_user_id != user_id:
            raise Exception("Only team admins can remove other members")

        # Get the member record
        response = supabase.table("team_members").select("id").eq("team_id", team_id).eq("user_id", user_id).execute()

        if response.error:
            raise Exception(f"Error getting team member: {response.error}")

        if not response.data:
            raise Exception("User is not a member of this team")

        member_id = response.data[0]["id"]

        # Remove the member
        result = cls.delete(member_id)

        # Log the activity
        if result:
            supabase.rpc(
                "log_team_activity",
                {
                    "p_team_id": team_id,
                    "p_user_id": current_user_id,
                    "p_action": "removed_member",
                    "p_resource_type": "team",
                    "p_resource_id": team_id,
                    "p_details": {"member_id": user_id}
                }
            ).execute()

            # Create notification for the removed user if they didn't remove themselves
            if current_user_id != user_id:
                supabase.table("notifications").insert({
                    "user_id": user_id,
                    "title": "Removed from Team",
                    "message": "You have been removed from a team",
                    "resource_type": "team",
                    "resource_id": team_id,
                    "created_by": current_user_id
                }).execute()

        return result

# SharedResource model
class SharedResource(BaseModel):
    """Model for shared resources in the database."""
    table_name = 'shared_resources'

    @classmethod
    def share_resource(cls, team_id: str, resource_type: str, resource_id: str) -> Dict[str, Any]:
        """
        Share a resource with a team.
        Args:
            team_id: The ID of the team
            resource_type: The type of resource ('project', 'mixture', 'experiment')
            resource_id: The ID of the resource
        Returns:
            The created shared resource record
        """
        user_id = get_user_id()

        # Check if the user has permission to share the resource
        supabase = cls.get_supabase()

        # Check if the user is a member of the team
        is_member_response = supabase.rpc(
            "is_team_member",
            {"p_team_id": team_id, "p_user_id": user_id}
        ).execute()

        if is_member_response.error:
            raise Exception(f"Error checking team membership: {is_member_response.error}")

        if not is_member_response.data:
            raise Exception("You must be a member of the team to share resources")

        # Share the resource
        shared_resource_data = {
            "team_id": team_id,
            "resource_type": resource_type,
            "resource_id": resource_id,
            "created_by": user_id
        }

        shared_resource = cls.create(shared_resource_data)

        # Log the activity
        if shared_resource:
            supabase.rpc(
                "log_team_activity",
                {
                    "p_team_id": team_id,
                    "p_user_id": user_id,
                    "p_action": "shared_resource",
                    "p_resource_type": resource_type,
                    "p_resource_id": resource_id,
                    "p_details": {"resource_type": resource_type}
                }
            ).execute()

            # Create notification for team members
            supabase.rpc(
                "create_team_notification",
                {
                    "p_team_id": team_id,
                    "p_title": "New Shared Resource",
                    "p_message": f"A new {resource_type} has been shared with your team",
                    "p_resource_type": resource_type,
                    "p_resource_id": resource_id,
                    "p_created_by": user_id
                }
            ).execute()

        return shared_resource

    @classmethod
    def get_team_shared_resources(cls, team_id: str, resource_type: str = None) -> List[Dict[str, Any]]:
        """
        Get all resources shared with a team.
        Args:
            team_id: The ID of the team
            resource_type: Optional filter by resource type
        Returns:
            List of shared resources
        """
        query = cls.get_table().select('*').eq('team_id', team_id)

        if resource_type:
            query = query.eq('resource_type', resource_type)

        response = query.execute()

        if response.error:
            raise Exception(f"Error getting shared resources: {response.error}")

        return response.data

    @classmethod
    def get_resource_teams(cls, resource_type: str, resource_id: str) -> List[Dict[str, Any]]:
        """
        Get all teams a resource is shared with.
        Args:
            resource_type: The type of resource
            resource_id: The ID of the resource
        Returns:
            List of teams
        """
        response = cls.get_table().select('team_id').eq('resource_type', resource_type).eq('resource_id', resource_id).execute()

        if response.error:
            raise Exception(f"Error getting resource teams: {response.error}")

        team_ids = [item['team_id'] for item in response.data]

        if not team_ids:
            return []

        teams_response = cls.get_supabase().table('teams').select('*').in_('id', team_ids).execute()

        if teams_response.error:
            raise Exception(f"Error getting teams: {teams_response.error}")

        return teams_response.data

    @classmethod
    def unshare_resource(cls, team_id: str, resource_type: str, resource_id: str) -> bool:
        """
        Unshare a resource from a team.
        Args:
            team_id: The ID of the team
            resource_type: The type of resource
            resource_id: The ID of the resource
        Returns:
            True if successful, False otherwise
        """
        user_id = get_user_id()

        # Check if the user has permission to unshare the resource
        supabase = cls.get_supabase()

        # Check if the user is an admin of the team
        is_admin_response = supabase.rpc(
            "has_team_role",
            {"p_team_id": team_id, "p_user_id": user_id, "p_role": "admin"}
        ).execute()

        if is_admin_response.error:
            raise Exception(f"Error checking team role: {is_admin_response.error}")

        if not is_admin_response.data:
            raise Exception("Only team admins can unshare resources")

        # Get the shared resource record
        response = supabase.table("shared_resources").select("id").eq("team_id", team_id).eq("resource_type", resource_type).eq("resource_id", resource_id).execute()

        if response.error:
            raise Exception(f"Error getting shared resource: {response.error}")

        if not response.data:
            raise Exception("Resource is not shared with this team")

        shared_resource_id = response.data[0]["id"]

        # Unshare the resource
        result = cls.delete(shared_resource_id)

        # Log the activity
        if result:
            supabase.rpc(
                "log_team_activity",
                {
                    "p_team_id": team_id,
                    "p_user_id": user_id,
                    "p_action": "unshared_resource",
                    "p_resource_type": resource_type,
                    "p_resource_id": resource_id,
                    "p_details": {"resource_type": resource_type}
                }
            ).execute()

        return result

# Comment model
class Comment(BaseModel):
    """Model for comments in the database."""
    table_name = 'comments'

    @classmethod
    def add_comment(cls, resource_type: str, resource_id: str, content: str, parent_id: str = None) -> Dict[str, Any]:
        """
        Add a comment to a resource.
        Args:
            resource_type: The type of resource ('project', 'mixture', 'experiment')
            resource_id: The ID of the resource
            content: The comment content
            parent_id: Optional parent comment ID for replies
        Returns:
            The created comment
        """
        user_id = get_user_id()

        # Check if the user has permission to comment on the resource
        # This would involve checking if the user owns the resource or is a member of a team that has access
        # For simplicity, we'll assume the RLS policies handle this

        comment_data = {
            "resource_type": resource_type,
            "resource_id": resource_id,
            "content": content,
            "parent_id": parent_id,
            "created_by": user_id
        }

        comment = cls.create(comment_data)

        # If this is a reply to another comment, create a notification for the parent comment author
        if parent_id and comment:
            supabase = cls.get_supabase()
            parent_comment = cls.get(parent_id)

            if parent_comment and parent_comment['created_by'] != user_id:
                supabase.table("notifications").insert({
                    "user_id": parent_comment['created_by'],
                    "title": "New Reply",
                    "message": "Someone replied to your comment",
                    "resource_type": "comment",
                    "resource_id": parent_id,
                    "created_by": user_id
                }).execute()

        return comment

    @classmethod
    def get_resource_comments(cls, resource_type: str, resource_id: str) -> List[Dict[str, Any]]:
        """
        Get all comments for a resource.
        Args:
            resource_type: The type of resource
            resource_id: The ID of the resource
        Returns:
            List of comments with user details
        """
        supabase = cls.get_supabase()

        # Get comments with user details
        response = supabase.from_("comments").select("""
            *,
            users:created_by (
                email,
                raw_user_meta_data
            )
        """).eq("resource_type", resource_type).eq("resource_id", resource_id).order("created_at", {"ascending": True}).execute()

        if response.error:
            raise Exception(f"Error getting comments: {response.error}")

        # Format the response to include user details
        comments = []
        for comment in response.data:
            user = comment.get('users', {})
            comments.append({
                "id": comment["id"],
                "resource_type": comment["resource_type"],
                "resource_id": comment["resource_id"],
                "content": comment["content"],
                "parent_id": comment["parent_id"],
                "created_at": comment["created_at"],
                "updated_at": comment["updated_at"],
                "created_by": comment["created_by"],
                "user_email": user.get("email"),
                "user_display_name": user.get("raw_user_meta_data", {}).get("name")
            })

        return comments

    @classmethod
    def update_comment(cls, comment_id: str, content: str) -> Dict[str, Any]:
        """
        Update a comment.
        Args:
            comment_id: The ID of the comment
            content: The new content
        Returns:
            The updated comment
        """
        user_id = get_user_id()

        # Check if the user is the author of the comment
        comment = cls.get(comment_id)

        if not comment:
            raise Exception("Comment not found")

        if comment['created_by'] != user_id:
            raise Exception("Only the author can update the comment")

        return cls.update(comment_id, {"content": content})

    @classmethod
    def delete_comment(cls, comment_id: str) -> bool:
        """
        Delete a comment.
        Args:
            comment_id: The ID of the comment
        Returns:
            True if successful, False otherwise
        """
        user_id = get_user_id()

        # Check if the user is the author of the comment
        comment = cls.get(comment_id)

        if not comment:
            raise Exception("Comment not found")

        if comment['created_by'] != user_id:
            raise Exception("Only the author can delete the comment")

        return cls.delete(comment_id)

# Notification model
class Notification(BaseModel):
    """Model for notifications in the database."""
    table_name = 'notifications'

    @classmethod
    def get_user_notifications(cls, user_id: str, limit: int = 20, offset: int = 0, unread_only: bool = False) -> List[Dict[str, Any]]:
        """
        Get notifications for a user.
        Args:
            user_id: The ID of the user
            limit: Maximum number of notifications to return
            offset: Number of notifications to skip
            unread_only: Whether to return only unread notifications
        Returns:
            List of notifications
        """
        query = cls.get_table().select('*').eq('user_id', user_id)

        if unread_only:
            query = query.eq('is_read', False)

        query = query.order('created_at', {'ascending': False}).range(offset, offset + limit - 1)

        response = query.execute()

        if response.error:
            raise Exception(f"Error getting notifications: {response.error}")

        return response.data

    @classmethod
    def mark_as_read(cls, notification_id: str) -> Dict[str, Any]:
        """
        Mark a notification as read.
        Args:
            notification_id: The ID of the notification
        Returns:
            The updated notification
        """
        return cls.update(notification_id, {"is_read": True})

    @classmethod
    def mark_all_as_read(cls, user_id: str) -> bool:
        """
        Mark all notifications for a user as read.
        Args:
            user_id: The ID of the user
        Returns:
            True if successful, False otherwise
        """
        supabase = cls.get_supabase()

        response = supabase.table("notifications").update({"is_read": True}).eq("user_id", user_id).eq("is_read", False).execute()

        if response.error:
            raise Exception(f"Error marking notifications as read: {response.error}")

        return True

    @classmethod
    def create_notification(cls, user_id: str, title: str, message: str, resource_type: str = None, resource_id: str = None, team_id: str = None) -> Dict[str, Any]:
        """
        Create a notification for a user.
        Args:
            user_id: The ID of the user
            title: The notification title
            message: The notification message
            resource_type: Optional resource type
            resource_id: Optional resource ID
            team_id: Optional team ID
        Returns:
            The created notification
        """
        current_user_id = get_user_id()

        notification_data = {
            "user_id": user_id,
            "title": title,
            "message": message,
            "resource_type": resource_type,
            "resource_id": resource_id,
            "team_id": team_id,
            "created_by": current_user_id
        }

        return cls.create(notification_data)

# ActivityLog model
class ActivityLog(BaseModel):
    """Model for activity logs in the database."""
    table_name = 'activity_log'

    @classmethod
    def get_team_activity(cls, team_id: str, limit: int = 20) -> List[Dict[str, Any]]:
        """
        Get activity log for a team.
        Args:
            team_id: The ID of the team
            limit: Maximum number of activities to return
        Returns:
            List of activities with user details
        """
        response = cls.get_supabase().rpc(
            "get_team_activity",
            {"p_team_id": team_id, "p_limit": limit}
        ).execute()

        if response.error:
            raise Exception(f"Error getting team activity: {response.error}")

        return response.data

    @classmethod
    def log_activity(cls, team_id: str, action: str, resource_type: str, resource_id: str, details: Dict[str, Any] = None) -> Dict[str, Any]:
        """
        Log an activity for a team.
        Args:
            team_id: The ID of the team
            action: The action performed
            resource_type: The type of resource
            resource_id: The ID of the resource
            details: Optional details about the activity
        Returns:
            The created activity log
        """
        user_id = get_user_id()

        response = cls.get_supabase().rpc(
            "log_team_activity",
            {
                "p_team_id": team_id,
                "p_user_id": user_id,
                "p_action": action,
                "p_resource_type": resource_type,
                "p_resource_id": resource_id,
                "p_details": details
            }
        ).execute()

        if response.error:
            raise Exception(f"Error logging activity: {response.error}")

        return {"id": response.data}
