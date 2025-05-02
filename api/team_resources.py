"""
CryoProtect Analyzer API - Team Resources

This module contains API resources for team collaboration features.
"""

from flask import request, g
from flask_restful import Resource
from flask_apispec import marshal_with, doc, use_kwargs
from flask_apispec.views import MethodResource
from marshmallow import Schema, fields

from api.utils import token_required, get_user_id, _handle_json_serialization, handle_supabase_error
from api.models import (
    BaseModel
)
from api.team_models import (
    team_fields, team_member_fields, shared_resource_fields,
    comment_fields, notification_fields, activity_log_fields,
    TeamSchema, TeamMemberSchema, SharedResourceSchema,
    CommentSchema, NotificationSchema,
    Team, TeamMember, SharedResource, Comment, Notification, ActivityLog
)

# Team resources
class TeamListResource(MethodResource, Resource):
    """Resource for listing and creating teams."""
    
    @token_required
    @doc(description='Get all teams for the current user', tags=['Teams'])
    @marshal_with(team_fields)
    def get(self):
        """Get all teams for the current user."""
        user_id = get_user_id()
        teams = Team.get_user_teams(user_id)
        return _handle_json_serialization(teams), 200
    
    @token_required
    @doc(description='Create a new team', tags=['Teams'])
    @use_kwargs(TeamSchema)
    @marshal_with(team_fields)
    def post(self, **kwargs):
        """Create a new team."""
        name = kwargs.get('name')
        description = kwargs.get('description')
        
        team = Team.create_with_creator_as_admin(name, description)
        return _handle_json_serialization(team), 201

class TeamResource(MethodResource, Resource):
    """Resource for managing a specific team."""
    
    @token_required
    @doc(description='Get a team by ID', tags=['Teams'])
    @marshal_with(team_fields)
    def get(self, team_id):
        """Get a team by ID."""
        team = Team.get_with_member_count(team_id)
        if not team:
            return _handle_json_serialization({'message': 'Team not found'}), 404
        return _handle_json_serialization(team), 200
    
    @token_required
    @doc(description='Update a team', tags=['Teams'])
    @use_kwargs(TeamSchema)
    @marshal_with(team_fields)
    def put(self, team_id, **kwargs):
        """Update a team."""
        team = Team.get(team_id)
        if not team:
            return _handle_json_serialization({'message': 'Team not found'}), 404
        
        # Check if user is admin
        user_id = get_user_id()
        supabase = BaseModel.get_supabase()
        is_admin_response = supabase.rpc(
            "has_team_role",
            {"p_team_id": team_id, "p_user_id": user_id, "p_role": "admin"}
        ).execute()
        
        if is_admin_response.error or not is_admin_response.data:
            return _handle_json_serialization({'message': 'Only team admins can update the team'}), 403
        
        updated_team = Team.update(team_id, kwargs)
        return _handle_json_serialization(updated_team), 200
    
    @token_required
    @doc(description='Delete a team', tags=['Teams'])
    def delete(self, team_id):
        """Delete a team."""
        team = Team.get(team_id)
        if not team:
            return _handle_json_serialization({'message': 'Team not found'}), 404
        
        # Check if user is admin
        user_id = get_user_id()
        supabase = BaseModel.get_supabase()
        is_admin_response = supabase.rpc(
            "has_team_role",
            {"p_team_id": team_id, "p_user_id": user_id, "p_role": "admin"}
        ).execute()
        
        if is_admin_response.error or not is_admin_response.data:
            return _handle_json_serialization({'message': 'Only team admins can delete the team'}), 403
        
        Team.delete(team_id)
        return _handle_json_serialization({'message': 'Team deleted successfully'}), 200

# Team member resources
class TeamMemberListResource(MethodResource, Resource):
    """Resource for listing and adding team members."""
    
    @token_required
    @doc(description='Get all members of a team', tags=['Team Members'])
    @marshal_with(team_member_fields)
    def get(self, team_id):
        """Get all members of a team."""
        members = TeamMember.get_team_members(team_id)
        return _handle_json_serialization(members), 200
    
    @token_required
    @doc(description='Add a member to a team', tags=['Team Members'])
    @use_kwargs(TeamMemberSchema)
    @marshal_with(team_member_fields)
    def post(self, team_id, **kwargs):
        """Add a member to a team."""
        user_id = kwargs.get('user_id')
        role = kwargs.get('role')
        
        try:
            member = TeamMember.add_member(team_id, user_id, role)
            return _handle_json_serialization(member), 201
        except Exception as e:
            return _handle_json_serialization(handle_supabase_error(e)), 400

class TeamMemberResource(MethodResource, Resource):
    """Resource for managing a specific team member."""
    
    @token_required
    @doc(description='Update a team member\'s role', tags=['Team Members'])
    @use_kwargs(TeamMemberSchema(only=['role']))
    @marshal_with(team_member_fields)
    def put(self, team_id, user_id, **kwargs):
        """Update a team member's role."""
        role = kwargs.get('role')
        
        try:
            member = TeamMember.update_role(team_id, user_id, role)
            return _handle_json_serialization(member), 200
        except Exception as e:
            return _handle_json_serialization(handle_supabase_error(e)), 400
    
    @token_required
    @doc(description='Remove a member from a team', tags=['Team Members'])
    def delete(self, team_id, user_id):
        """Remove a member from a team."""
        try:
            TeamMember.remove_member(team_id, user_id)
            return _handle_json_serialization({'message': 'Member removed successfully'}), 200
        except Exception as e:
            return _handle_json_serialization(handle_supabase_error(e)), 400

# Shared resource resources
class SharedResourceListResource(MethodResource, Resource):
    """Resource for listing and sharing resources with teams."""
    
    @token_required
    @doc(description='Get all resources shared with a team', tags=['Shared Resources'])
    @marshal_with(shared_resource_fields)
    def get(self, team_id):
        """Get all resources shared with a team."""
        resource_type = request.args.get('resource_type')
        resources = SharedResource.get_team_shared_resources(team_id, resource_type)
        return _handle_json_serialization(resources), 200
    
    @token_required
    @doc(description='Share a resource with a team', tags=['Shared Resources'])
    @use_kwargs(SharedResourceSchema(only=['resource_type', 'resource_id']))
    @marshal_with(shared_resource_fields)
    def post(self, team_id, **kwargs):
        """Share a resource with a team."""
        resource_type = kwargs.get('resource_type')
        resource_id = kwargs.get('resource_id')
        
        try:
            shared_resource = SharedResource.share_resource(team_id, resource_type, resource_id)
            return _handle_json_serialization(shared_resource), 201
        except Exception as e:
            return _handle_json_serialization(handle_supabase_error(e)), 400

class SharedResourceResource(MethodResource, Resource):
    """Resource for managing a specific shared resource."""
    
    @token_required
    @doc(description='Unshare a resource from a team', tags=['Shared Resources'])
    def delete(self, team_id, resource_type, resource_id):
        """Unshare a resource from a team."""
        try:
            SharedResource.unshare_resource(team_id, resource_type, resource_id)
            return _handle_json_serialization({'message': 'Resource unshared successfully'}), 200
        except Exception as e:
            return _handle_json_serialization(handle_supabase_error(e)), 400

# Comment resources
class CommentListResource(MethodResource, Resource):
    """Resource for listing and creating comments."""
    
    @token_required
    @doc(description='Get all comments for a resource', tags=['Comments'])
    @marshal_with(comment_fields)
    def get(self, resource_type, resource_id):
        """Get all comments for a resource."""
        comments = Comment.get_resource_comments(resource_type, resource_id)
        return _handle_json_serialization(comments), 200
    
    @token_required
    @doc(description='Add a comment to a resource', tags=['Comments'])
    @use_kwargs(CommentSchema(only=['content', 'parent_id']))
    @marshal_with(comment_fields)
    def post(self, resource_type, resource_id, **kwargs):
        """Add a comment to a resource."""
        content = kwargs.get('content')
        parent_id = kwargs.get('parent_id')
        
        try:
            comment = Comment.add_comment(resource_type, resource_id, content, parent_id)
            return _handle_json_serialization(comment), 201
        except Exception as e:
            return _handle_json_serialization(handle_supabase_error(e)), 400

class CommentResource(MethodResource, Resource):
    """Resource for managing a specific comment."""
    
    @token_required
    @doc(description='Update a comment', tags=['Comments'])
    @use_kwargs(CommentSchema(only=['content']))
    @marshal_with(comment_fields)
    def put(self, comment_id, **kwargs):
        """Update a comment."""
        content = kwargs.get('content')
        
        try:
            comment = Comment.update_comment(comment_id, content)
            return _handle_json_serialization(comment), 200
        except Exception as e:
            return _handle_json_serialization(handle_supabase_error(e)), 400
    
    @token_required
    @doc(description='Delete a comment', tags=['Comments'])
    def delete(self, comment_id):
        """Delete a comment."""
        try:
            Comment.delete_comment(comment_id)
            return _handle_json_serialization({'message': 'Comment deleted successfully'}), 200
        except Exception as e:
            return _handle_json_serialization(handle_supabase_error(e)), 400

# Notification resources
class NotificationListResource(MethodResource, Resource):
    """Resource for listing and managing notifications."""
    
    @token_required
    @doc(description='Get notifications for the current user', tags=['Notifications'])
    @marshal_with(notification_fields)
    def get(self):
        """Get notifications for the current user."""
        user_id = get_user_id()
        limit = int(request.args.get('limit', 20))
        offset = int(request.args.get('offset', 0))
        unread_only = request.args.get('unread_only', 'false').lower() == 'true'
        
        notifications = Notification.get_user_notifications(user_id, limit, offset, unread_only)
        return _handle_json_serialization(notifications), 200
    
    @token_required
    @doc(description='Mark all notifications as read', tags=['Notifications'])
    def put(self):
        """Mark all notifications as read."""
        user_id = get_user_id()
        
        try:
            Notification.mark_all_as_read(user_id)
            return _handle_json_serialization({'message': 'All notifications marked as read'}), 200
        except Exception as e:
            return _handle_json_serialization(handle_supabase_error(e)), 400

class NotificationResource(MethodResource, Resource):
    """Resource for managing a specific notification."""
    
    @token_required
    @doc(description='Mark a notification as read', tags=['Notifications'])
    def put(self, notification_id):
        """Mark a notification as read."""
        try:
            Notification.mark_as_read(notification_id)
            return _handle_json_serialization({'message': 'Notification marked as read'}), 200
        except Exception as e:
            return _handle_json_serialization(handle_supabase_error(e)), 400
    
    @token_required
    @doc(description='Delete a notification', tags=['Notifications'])
    def delete(self, notification_id):
        """Delete a notification."""
        try:
            Notification.delete(notification_id)
            return _handle_json_serialization({'message': 'Notification deleted successfully'}), 200
        except Exception as e:
            return _handle_json_serialization(handle_supabase_error(e)), 400

# Activity log resources
class ActivityLogResource(MethodResource, Resource):
    """Resource for viewing team activity logs."""
    
    @token_required
    @doc(description='Get activity log for a team', tags=['Activity Log'])
    @marshal_with(activity_log_fields)
    def get(self, team_id):
        """Get activity log for a team."""
        limit = int(request.args.get('limit', 20))
        
        try:
            activities = ActivityLog.get_team_activity(team_id, limit)
            return _handle_json_serialization(activities), 200
        except Exception as e:
            return _handle_json_serialization(handle_supabase_error(e)), 400

def register_resources(api):
    """Register team resources with the API."""
    api.add_resource(TeamListResource, '/teams')
    api.add_resource(TeamResource, '/teams/<string:team_id>')
    api.add_resource(TeamMemberListResource, '/teams/<string:team_id>/members')
    api.add_resource(TeamMemberResource, '/teams/<string:team_id>/members/<string:user_id>')
    api.add_resource(SharedResourceListResource, '/teams/<string:team_id>/resources')
    api.add_resource(SharedResourceResource, '/teams/<string:team_id>/resources/<string:resource_type>/<string:resource_id>')
    api.add_resource(CommentListResource, '/comments/<string:resource_type>/<string:resource_id>')
    api.add_resource(CommentResource, '/comments/<string:comment_id>')
    api.add_resource(NotificationListResource, '/notifications')
    api.add_resource(NotificationResource, '/notifications/<string:notification_id>')
    api.add_resource(ActivityLogResource, '/teams/<string:team_id>/activity')