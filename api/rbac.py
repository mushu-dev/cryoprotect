"""
CryoProtect v2 - Role-Based Access Control (RBAC)

This module provides functionality for role-based access control, including:
- Role and permission management
- Role and permission verification
- Decorators for protecting routes based on roles and permissions
"""

import os
import logging
import functools
from typing import Dict, List, Union, Optional, Any, Callable
from flask import request, g, jsonify, current_app
from .jwt_auth import jwt_required, get_current_user, extract_user_from_token

logger = logging.getLogger(__name__)

def get_supabase():
    """Get the Supabase client from the current app context."""
    if hasattr(g, 'supabase'):
        return g.supabase
    
    from supabase import create_client
    
    supabase_url = os.environ.get('SUPABASE_URL') or current_app.config.get('SUPABASE_URL')
    supabase_key = os.environ.get('SUPABASE_KEY') or current_app.config.get('SUPABASE_KEY')
    
    if not supabase_url or not supabase_key:
        raise ValueError("SUPABASE_URL and SUPABASE_KEY must be set")
    
    g.supabase = create_client(supabase_url, supabase_key)
    return g.supabase

class RoleManager:
    """
    Manages roles in the system.
    """
    
    @staticmethod
    def get_all_roles() -> List[Dict[str, Any]]:
        """
        Get all roles in the system.
        
        Returns:
            List of role objects
        """
        supabase = get_supabase()
        response = supabase.table("roles").select("*").execute()
        
        if response.error:
            logger.error(f"Error getting roles: {response.error}")
            raise Exception(f"Error getting roles: {response.error}")
        
        return response.data
    
    @staticmethod
    def get_role(role_id: str) -> Optional[Dict[str, Any]]:
        """
        Get a role by ID.
        
        Args:
            role_id: The ID of the role
            
        Returns:
            The role object or None if not found
        """
        supabase = get_supabase()
        response = supabase.table("roles").select("*").eq("id", role_id).execute()
        
        if response.error:
            logger.error(f"Error getting role: {response.error}")
            raise Exception(f"Error getting role: {response.error}")
        
        if not response.data:
            return None
        
        return response.data[0]
    
    @staticmethod
    def get_role_by_name(name: str) -> Optional[Dict[str, Any]]:
        """
        Get a role by name.
        
        Args:
            name: The name of the role
            
        Returns:
            The role object or None if not found
        """
        supabase = get_supabase()
        response = supabase.table("roles").select("*").eq("name", name).execute()
        
        if response.error:
            logger.error(f"Error getting role: {response.error}")
            raise Exception(f"Error getting role: {response.error}")
        
        if not response.data:
            return None
        
        return response.data[0]
    
    @staticmethod
    def create_role(name: str, description: str = None, is_system_role: bool = False) -> Dict[str, Any]:
        """
        Create a new role.
        
        Args:
            name: The name of the role
            description: The description of the role
            is_system_role: Whether this is a system role
            
        Returns:
            The created role object
        """
        user = get_current_user()
        if not user:
            raise Exception("User not authenticated")
        
        supabase = get_supabase()
        response = supabase.table("roles").insert({
            "name": name,
            "description": description,
            "is_system_role": is_system_role,
            "created_by": user.get("id")
        }).execute()
        
        if response.error:
            logger.error(f"Error creating role: {response.error}")
            raise Exception(f"Error creating role: {response.error}")
        
        return response.data[0]
    
    @staticmethod
    def update_role(role_id: str, data: Dict[str, Any]) -> Dict[str, Any]:
        """
        Update a role.
        
        Args:
            role_id: The ID of the role to update
            data: The data to update
            
        Returns:
            The updated role object
        """
        supabase = get_supabase()
        response = supabase.table("roles").update(data).eq("id", role_id).execute()
        
        if response.error:
            logger.error(f"Error updating role: {response.error}")
            raise Exception(f"Error updating role: {response.error}")
        
        if not response.data:
            raise Exception(f"Role with ID {role_id} not found")
        
        return response.data[0]
    
    @staticmethod
    def delete_role(role_id: str) -> bool:
        """
        Delete a role.
        
        Args:
            role_id: The ID of the role to delete
            
        Returns:
            True if successful, False otherwise
        """
        # Check if it's a system role
        role = RoleManager.get_role(role_id)
        if role and role.get("is_system_role"):
            raise Exception("Cannot delete system roles")
        
        supabase = get_supabase()
        response = supabase.table("roles").delete().eq("id", role_id).execute()
        
        if response.error:
            logger.error(f"Error deleting role: {response.error}")
            raise Exception(f"Error deleting role: {response.error}")
        
        return bool(response.data)
    
    @staticmethod
    def get_role_permissions(role_id: str) -> List[Dict[str, Any]]:
        """
        Get all permissions for a role.
        
        Args:
            role_id: The ID of the role
            
        Returns:
            List of permission objects
        """
        supabase = get_supabase()
        response = supabase.rpc(
            "get_role_permissions",
            {"p_role_id": role_id}
        ).execute()
        
        if response.error:
            logger.error(f"Error getting role permissions: {response.error}")
            raise Exception(f"Error getting role permissions: {response.error}")
        
        return response.data
    
    @staticmethod
    def assign_permission_to_role(role_id: str, permission_id: str) -> Dict[str, Any]:
        """
        Assign a permission to a role.
        
        Args:
            role_id: The ID of the role
            permission_id: The ID of the permission
            
        Returns:
            The created role_permission object
        """
        user = get_current_user()
        if not user:
            raise Exception("User not authenticated")
        
        supabase = get_supabase()
        response = supabase.table("role_permissions").insert({
            "role_id": role_id,
            "permission_id": permission_id,
            "created_by": user.get("id")
        }).execute()
        
        if response.error:
            logger.error(f"Error assigning permission to role: {response.error}")
            raise Exception(f"Error assigning permission to role: {response.error}")
        
        return response.data[0]
    
    @staticmethod
    def remove_permission_from_role(role_id: str, permission_id: str) -> bool:
        """
        Remove a permission from a role.
        
        Args:
            role_id: The ID of the role
            permission_id: The ID of the permission
            
        Returns:
            True if successful, False otherwise
        """
        supabase = get_supabase()
        response = supabase.table("role_permissions").delete().eq("role_id", role_id).eq("permission_id", permission_id).execute()
        
        if response.error:
            logger.error(f"Error removing permission from role: {response.error}")
            raise Exception(f"Error removing permission from role: {response.error}")
        
        return bool(response.data)


class PermissionManager:
    """
    Manages permissions in the system.
    """
    
    @staticmethod
    def get_all_permissions() -> List[Dict[str, Any]]:
        """
        Get all permissions in the system.
        
        Returns:
            List of permission objects
        """
        supabase = get_supabase()
        response = supabase.table("permissions").select("*").execute()
        
        if response.error:
            logger.error(f"Error getting permissions: {response.error}")
            raise Exception(f"Error getting permissions: {response.error}")
        
        return response.data
    
    @staticmethod
    def get_permission(permission_id: str) -> Optional[Dict[str, Any]]:
        """
        Get a permission by ID.
        
        Args:
            permission_id: The ID of the permission
            
        Returns:
            The permission object or None if not found
        """
        supabase = get_supabase()
        response = supabase.table("permissions").select("*").eq("id", permission_id).execute()
        
        if response.error:
            logger.error(f"Error getting permission: {response.error}")
            raise Exception(f"Error getting permission: {response.error}")
        
        if not response.data:
            return None
        
        return response.data[0]
    
    @staticmethod
    def get_permission_by_resource_action(resource_type: str, action: str) -> Optional[Dict[str, Any]]:
        """
        Get a permission by resource type and action.
        
        Args:
            resource_type: The resource type
            action: The action
            
        Returns:
            The permission object or None if not found
        """
        supabase = get_supabase()
        response = supabase.table("permissions").select("*").eq("resource_type", resource_type).eq("action", action).execute()
        
        if response.error:
            logger.error(f"Error getting permission: {response.error}")
            raise Exception(f"Error getting permission: {response.error}")
        
        if not response.data:
            return None
        
        return response.data[0]
    
    @staticmethod
    def create_permission(name: str, resource_type: str, action: str, description: str = None) -> Dict[str, Any]:
        """
        Create a new permission.
        
        Args:
            name: The name of the permission
            resource_type: The resource type
            action: The action
            description: The description of the permission
            
        Returns:
            The created permission object
        """
        user = get_current_user()
        if not user:
            raise Exception("User not authenticated")
        
        supabase = get_supabase()
        response = supabase.table("permissions").insert({
            "name": name,
            "resource_type": resource_type,
            "action": action,
            "description": description,
            "created_by": user.get("id")
        }).execute()
        
        if response.error:
            logger.error(f"Error creating permission: {response.error}")
            raise Exception(f"Error creating permission: {response.error}")
        
        return response.data[0]
    
    @staticmethod
    def update_permission(permission_id: str, data: Dict[str, Any]) -> Dict[str, Any]:
        """
        Update a permission.
        
        Args:
            permission_id: The ID of the permission to update
            data: The data to update
            
        Returns:
            The updated permission object
        """
        supabase = get_supabase()
        response = supabase.table("permissions").update(data).eq("id", permission_id).execute()
        
        if response.error:
            logger.error(f"Error updating permission: {response.error}")
            raise Exception(f"Error updating permission: {response.error}")
        
        if not response.data:
            raise Exception(f"Permission with ID {permission_id} not found")
        
        return response.data[0]
    
    @staticmethod
    def delete_permission(permission_id: str) -> bool:
        """
        Delete a permission.
        
        Args:
            permission_id: The ID of the permission to delete
            
        Returns:
            True if successful, False otherwise
        """
        supabase = get_supabase()
        response = supabase.table("permissions").delete().eq("id", permission_id).execute()
        
        if response.error:
            logger.error(f"Error deleting permission: {response.error}")
            raise Exception(f"Error deleting permission: {response.error}")
        
        return bool(response.data)


class UserRoleManager:
    """
    Manages user roles in the system.
    """
    
    @staticmethod
    def get_user_roles(user_id: str) -> List[Dict[str, Any]]:
        """
        Get all roles for a user.
        
        Args:
            user_id: The ID of the user
            
        Returns:
            List of role objects
        """
        supabase = get_supabase()
        response = supabase.rpc(
            "get_user_roles",
            {"p_user_id": user_id}
        ).execute()
        
        if response.error:
            logger.error(f"Error getting user roles: {response.error}")
            raise Exception(f"Error getting user roles: {response.error}")
        
        return response.data
    
    @staticmethod
    def get_user_permissions(user_id: str) -> List[Dict[str, Any]]:
        """
        Get all permissions for a user.
        
        Args:
            user_id: The ID of the user
            
        Returns:
            List of permission objects
        """
        supabase = get_supabase()
        response = supabase.rpc(
            "get_user_permissions",
            {"p_user_id": user_id}
        ).execute()
        
        if response.error:
            logger.error(f"Error getting user permissions: {response.error}")
            raise Exception(f"Error getting user permissions: {response.error}")
        
        return response.data
    
    @staticmethod
    def assign_role_to_user(user_id: str, role_id: str) -> Dict[str, Any]:
        """
        Assign a role to a user.
        
        Args:
            user_id: The ID of the user
            role_id: The ID of the role
            
        Returns:
            The created user_role object
        """
        current_user = get_current_user()
        if not current_user:
            raise Exception("User not authenticated")
        
        supabase = get_supabase()
        response = supabase.table("user_roles").insert({
            "user_id": user_id,
            "role_id": role_id,
            "created_by": current_user.get("id")
        }).execute()
        
        if response.error:
            logger.error(f"Error assigning role to user: {response.error}")
            raise Exception(f"Error assigning role to user: {response.error}")
        
        return response.data[0]
    
    @staticmethod
    def remove_role_from_user(user_id: str, role_id: str) -> bool:
        """
        Remove a role from a user.
        
        Args:
            user_id: The ID of the user
            role_id: The ID of the role
            
        Returns:
            True if successful, False otherwise
        """
        supabase = get_supabase()
        response = supabase.table("user_roles").delete().eq("user_id", user_id).eq("role_id", role_id).execute()
        
        if response.error:
            logger.error(f"Error removing role from user: {response.error}")
            raise Exception(f"Error removing role from user: {response.error}")
        
        return bool(response.data)
    
    @staticmethod
    def has_role(user_id: str, role_name: str) -> bool:
        """
        Check if a user has a specific role.
        
        Args:
            user_id: The ID of the user
            role_name: The name of the role
            
        Returns:
            True if the user has the role, False otherwise
        """
        supabase = get_supabase()
        response = supabase.rpc(
            "has_role",
            {"p_user_id": user_id, "p_role_name": role_name}
        ).execute()
        
        if response.error:
            logger.error(f"Error checking if user has role: {response.error}")
            raise Exception(f"Error checking if user has role: {response.error}")
        
        return response.data
    
    @staticmethod
    def has_permission(user_id: str, resource_type: str, action: str) -> bool:
        """
        Check if a user has a specific permission.
        
        Args:
            user_id: The ID of the user
            resource_type: The resource type
            action: The action
            
        Returns:
            True if the user has the permission, False otherwise
        """
        supabase = get_supabase()
        response = supabase.rpc(
            "has_permission",
            {"p_user_id": user_id, "p_resource_type": resource_type, "p_action": action}
        ).execute()
        
        if response.error:
            logger.error(f"Error checking if user has permission: {response.error}")
            raise Exception(f"Error checking if user has permission: {response.error}")
        
        return response.data


# Decorator functions for route protection

def role_required(required_role: Union[str, List[str]]):
    """
    Decorator to require a specific role for API endpoints.
    
    Args:
        required_role: Role or list of roles required
        
    Returns:
        Decorated function
    """
    def decorator(f):
        @functools.wraps(f)
        @jwt_required
        def decorated(*args, **kwargs):
            user = get_current_user()
            if not user:
                return jsonify({'message': 'Authentication required'}), 401
            
            user_id = user.get('id')
            
            # Convert single role to list
            roles = required_role if isinstance(required_role, list) else [required_role]
            
            # Check if user has any of the required roles
            has_required_role = False
            for role in roles:
                if UserRoleManager.has_role(user_id, role):
                    has_required_role = True
                    break
            
            if not has_required_role:
                return jsonify({'message': 'Insufficient permissions'}), 403
            
            return f(*args, **kwargs)
        
        return decorated
    
    return decorator


def permission_required(resource_type: str, action: str):
    """
    Decorator to require a specific permission for API endpoints.
    
    Args:
        resource_type: The resource type
        action: The action
        
    Returns:
        Decorated function
    """
    def decorator(f):
        @functools.wraps(f)
        @jwt_required
        def decorated(*args, **kwargs):
            user = get_current_user()
            if not user:
                return jsonify({'message': 'Authentication required'}), 401
            
            user_id = user.get('id')
            
            # Check if user has the required permission
            if not UserRoleManager.has_permission(user_id, resource_type, action):
                return jsonify({'message': 'Insufficient permissions'}), 403
            
            return f(*args, **kwargs)
        
        return decorated
    
    return decorator


def admin_required(f):
    """
    Decorator to require admin role for API endpoints.
    
    Args:
        f: Function to decorate
        
    Returns:
        Decorated function
    """
    return role_required('admin')(f)


# Helper functions

def get_current_user_roles():
    """
    Get the roles of the current user.
    
    Returns:
        List of role objects or None if not authenticated
    """
    user = get_current_user()
    if not user:
        return None
    
    return UserRoleManager.get_user_roles(user.get('id'))


def get_current_user_permissions():
    """
    Get the permissions of the current user.
    
    Returns:
        List of permission objects or None if not authenticated
    """
    user = get_current_user()
    if not user:
        return None
    
    return UserRoleManager.get_user_permissions(user.get('id'))


def current_user_has_role(role_name: str) -> bool:
    """
    Check if the current user has a specific role.
    
    Args:
        role_name: The name of the role
        
    Returns:
        True if the user has the role, False otherwise
    """
    user = get_current_user()
    if not user:
        return False
    
    return UserRoleManager.has_role(user.get('id'), role_name)


def current_user_has_permission(resource_type: str, action: str) -> bool:
    """
    Check if the current user has a specific permission.
    
    Args:
        resource_type: The resource type
        action: The action
        
    Returns:
        True if the user has the permission, False otherwise
    """
    user = get_current_user()
    if not user:
        return False
    
    return UserRoleManager.has_permission(user.get('id'), resource_type, action)