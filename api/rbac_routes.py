"""
CryoProtect v2 - Role-Based Access Control (RBAC) API Routes

This module provides API endpoints for managing roles and permissions.
"""

from flask import Blueprint, request, jsonify, g
from .rbac import (
    RoleManager, PermissionManager, UserRoleManager,
    admin_required, role_required, permission_required
)
from .jwt_auth import jwt_required, get_current_user

rbac_bp = Blueprint('rbac', __name__)

# Role management endpoints

@rbac_bp.route('/roles', methods=['GET'])
@jwt_required
@permission_required('roles', 'read')
def get_roles():
    """Get all roles."""
    try:
        roles = RoleManager.get_all_roles()
        return jsonify(roles), 200
    except Exception as e:
        return jsonify({'message': str(e)}), 500


@rbac_bp.route('/roles/<role_id>', methods=['GET'])
@jwt_required
@permission_required('roles', 'read')
def get_role(role_id):
    """Get a role by ID."""
    try:
        role = RoleManager.get_role(role_id)
        if not role:
            return jsonify({'message': 'Role not found'}), 404
        return jsonify(role), 200
    except Exception as e:
        return jsonify({'message': str(e)}), 500


@rbac_bp.route('/roles', methods=['POST'])
@jwt_required
@permission_required('roles', 'create')
def create_role():
    """Create a new role."""
    try:
        data = request.get_json()
        if not data:
            return jsonify({'message': 'No input data provided'}), 400
        
        name = data.get('name')
        description = data.get('description')
        is_system_role = data.get('is_system_role', False)
        
        if not name:
            return jsonify({'message': 'Role name is required'}), 400
        
        # Check if role already exists
        existing_role = RoleManager.get_role_by_name(name)
        if existing_role:
            return jsonify({'message': f'Role with name {name} already exists'}), 409
        
        role = RoleManager.create_role(name, description, is_system_role)
        return jsonify(role), 201
    except Exception as e:
        return jsonify({'message': str(e)}), 500


@rbac_bp.route('/roles/<role_id>', methods=['PUT'])
@jwt_required
@permission_required('roles', 'update')
def update_role(role_id):
    """Update a role."""
    try:
        data = request.get_json()
        if not data:
            return jsonify({'message': 'No input data provided'}), 400
        
        # Check if role exists
        role = RoleManager.get_role(role_id)
        if not role:
            return jsonify({'message': 'Role not found'}), 404
        
        # Don't allow changing system role status
        if 'is_system_role' in data and role.get('is_system_role') != data.get('is_system_role'):
            return jsonify({'message': 'Cannot change system role status'}), 403
        
        updated_role = RoleManager.update_role(role_id, data)
        return jsonify(updated_role), 200
    except Exception as e:
        return jsonify({'message': str(e)}), 500


@rbac_bp.route('/roles/<role_id>', methods=['DELETE'])
@jwt_required
@permission_required('roles', 'delete')
def delete_role(role_id):
    """Delete a role."""
    try:
        # Check if role exists
        role = RoleManager.get_role(role_id)
        if not role:
            return jsonify({'message': 'Role not found'}), 404
        
        # Don't allow deleting system roles
        if role.get('is_system_role'):
            return jsonify({'message': 'Cannot delete system roles'}), 403
        
        success = RoleManager.delete_role(role_id)
        if success:
            return jsonify({'message': 'Role deleted successfully'}), 200
        else:
            return jsonify({'message': 'Failed to delete role'}), 500
    except Exception as e:
        return jsonify({'message': str(e)}), 500


@rbac_bp.route('/roles/<role_id>/permissions', methods=['GET'])
@jwt_required
@permission_required('roles', 'read')
def get_role_permissions(role_id):
    """Get all permissions for a role."""
    try:
        # Check if role exists
        role = RoleManager.get_role(role_id)
        if not role:
            return jsonify({'message': 'Role not found'}), 404
        
        permissions = RoleManager.get_role_permissions(role_id)
        return jsonify(permissions), 200
    except Exception as e:
        return jsonify({'message': str(e)}), 500


@rbac_bp.route('/roles/<role_id>/permissions/<permission_id>', methods=['POST'])
@jwt_required
@permission_required('roles', 'update')
def assign_permission_to_role(role_id, permission_id):
    """Assign a permission to a role."""
    try:
        # Check if role exists
        role = RoleManager.get_role(role_id)
        if not role:
            return jsonify({'message': 'Role not found'}), 404
        
        # Check if permission exists
        permission = PermissionManager.get_permission(permission_id)
        if not permission:
            return jsonify({'message': 'Permission not found'}), 404
        
        role_permission = RoleManager.assign_permission_to_role(role_id, permission_id)
        return jsonify(role_permission), 201
    except Exception as e:
        return jsonify({'message': str(e)}), 500


@rbac_bp.route('/roles/<role_id>/permissions/<permission_id>', methods=['DELETE'])
@jwt_required
@permission_required('roles', 'update')
def remove_permission_from_role(role_id, permission_id):
    """Remove a permission from a role."""
    try:
        # Check if role exists
        role = RoleManager.get_role(role_id)
        if not role:
            return jsonify({'message': 'Role not found'}), 404
        
        # Check if permission exists
        permission = PermissionManager.get_permission(permission_id)
        if not permission:
            return jsonify({'message': 'Permission not found'}), 404
        
        success = RoleManager.remove_permission_from_role(role_id, permission_id)
        if success:
            return jsonify({'message': 'Permission removed from role successfully'}), 200
        else:
            return jsonify({'message': 'Failed to remove permission from role'}), 500
    except Exception as e:
        return jsonify({'message': str(e)}), 500


# Permission management endpoints

@rbac_bp.route('/permissions', methods=['GET'])
@jwt_required
@permission_required('permissions', 'read')
def get_permissions():
    """Get all permissions."""
    try:
        permissions = PermissionManager.get_all_permissions()
        return jsonify(permissions), 200
    except Exception as e:
        return jsonify({'message': str(e)}), 500


@rbac_bp.route('/permissions/<permission_id>', methods=['GET'])
@jwt_required
@permission_required('permissions', 'read')
def get_permission(permission_id):
    """Get a permission by ID."""
    try:
        permission = PermissionManager.get_permission(permission_id)
        if not permission:
            return jsonify({'message': 'Permission not found'}), 404
        return jsonify(permission), 200
    except Exception as e:
        return jsonify({'message': str(e)}), 500


@rbac_bp.route('/permissions', methods=['POST'])
@jwt_required
@permission_required('permissions', 'create')
def create_permission():
    """Create a new permission."""
    try:
        data = request.get_json()
        if not data:
            return jsonify({'message': 'No input data provided'}), 400
        
        name = data.get('name')
        resource_type = data.get('resource_type')
        action = data.get('action')
        description = data.get('description')
        
        if not name or not resource_type or not action:
            return jsonify({'message': 'Name, resource_type, and action are required'}), 400
        
        # Check if permission already exists
        existing_permission = PermissionManager.get_permission_by_resource_action(resource_type, action)
        if existing_permission:
            return jsonify({'message': f'Permission for {resource_type}:{action} already exists'}), 409
        
        permission = PermissionManager.create_permission(name, resource_type, action, description)
        return jsonify(permission), 201
    except Exception as e:
        return jsonify({'message': str(e)}), 500


@rbac_bp.route('/permissions/<permission_id>', methods=['PUT'])
@jwt_required
@permission_required('permissions', 'update')
def update_permission(permission_id):
    """Update a permission."""
    try:
        data = request.get_json()
        if not data:
            return jsonify({'message': 'No input data provided'}), 400
        
        # Check if permission exists
        permission = PermissionManager.get_permission(permission_id)
        if not permission:
            return jsonify({'message': 'Permission not found'}), 404
        
        # If resource_type or action is being changed, check for conflicts
        if ('resource_type' in data and data['resource_type'] != permission['resource_type']) or \
           ('action' in data and data['action'] != permission['action']):
            existing = PermissionManager.get_permission_by_resource_action(
                data.get('resource_type', permission['resource_type']),
                data.get('action', permission['action'])
            )
            if existing and existing['id'] != permission_id:
                return jsonify({'message': 'Permission with this resource_type and action already exists'}), 409
        
        updated_permission = PermissionManager.update_permission(permission_id, data)
        return jsonify(updated_permission), 200
    except Exception as e:
        return jsonify({'message': str(e)}), 500


@rbac_bp.route('/permissions/<permission_id>', methods=['DELETE'])
@jwt_required
@permission_required('permissions', 'delete')
def delete_permission(permission_id):
    """Delete a permission."""
    try:
        # Check if permission exists
        permission = PermissionManager.get_permission(permission_id)
        if not permission:
            return jsonify({'message': 'Permission not found'}), 404
        
        success = PermissionManager.delete_permission(permission_id)
        if success:
            return jsonify({'message': 'Permission deleted successfully'}), 200
        else:
            return jsonify({'message': 'Failed to delete permission'}), 500
    except Exception as e:
        return jsonify({'message': str(e)}), 500


# User role management endpoints

@rbac_bp.route('/users/<user_id>/roles', methods=['GET'])
@jwt_required
@permission_required('users', 'read')
def get_user_roles(user_id):
    """Get all roles for a user."""
    try:
        roles = UserRoleManager.get_user_roles(user_id)
        return jsonify(roles), 200
    except Exception as e:
        return jsonify({'message': str(e)}), 500


@rbac_bp.route('/users/<user_id>/permissions', methods=['GET'])
@jwt_required
@permission_required('users', 'read')
def get_user_permissions(user_id):
    """Get all permissions for a user."""
    try:
        permissions = UserRoleManager.get_user_permissions(user_id)
        return jsonify(permissions), 200
    except Exception as e:
        return jsonify({'message': str(e)}), 500


@rbac_bp.route('/users/<user_id>/roles/<role_id>', methods=['POST'])
@jwt_required
@permission_required('users', 'manage')
def assign_role_to_user(user_id, role_id):
    """Assign a role to a user."""
    try:
        # Check if role exists
        role = RoleManager.get_role(role_id)
        if not role:
            return jsonify({'message': 'Role not found'}), 404
        
        # Check if user already has this role
        if UserRoleManager.has_role(user_id, role['name']):
            return jsonify({'message': f'User already has the role {role["name"]}'}), 409
        
        user_role = UserRoleManager.assign_role_to_user(user_id, role_id)
        return jsonify(user_role), 201
    except Exception as e:
        return jsonify({'message': str(e)}), 500


@rbac_bp.route('/users/<user_id>/roles/<role_id>', methods=['DELETE'])
@jwt_required
@permission_required('users', 'manage')
def remove_role_from_user(user_id, role_id):
    """Remove a role from a user."""
    try:
        # Check if role exists
        role = RoleManager.get_role(role_id)
        if not role:
            return jsonify({'message': 'Role not found'}), 404
        
        # Don't allow removing the last role from a user
        user_roles = UserRoleManager.get_user_roles(user_id)
        if len(user_roles) <= 1:
            return jsonify({'message': 'Cannot remove the last role from a user'}), 403
        
        success = UserRoleManager.remove_role_from_user(user_id, role_id)
        if success:
            return jsonify({'message': 'Role removed from user successfully'}), 200
        else:
            return jsonify({'message': 'Failed to remove role from user'}), 500
    except Exception as e:
        return jsonify({'message': str(e)}), 500


# Current user endpoints

@rbac_bp.route('/me/roles', methods=['GET'])
@jwt_required
def get_my_roles():
    """Get all roles for the current user."""
    try:
        user = get_current_user()
        if not user:
            return jsonify({'message': 'User not authenticated'}), 401
        
        roles = UserRoleManager.get_user_roles(user['id'])
        return jsonify(roles), 200
    except Exception as e:
        return jsonify({'message': str(e)}), 500


@rbac_bp.route('/me/permissions', methods=['GET'])
@jwt_required
def get_my_permissions():
    """Get all permissions for the current user."""
    try:
        user = get_current_user()
        if not user:
            return jsonify({'message': 'User not authenticated'}), 401
        
        permissions = UserRoleManager.get_user_permissions(user['id'])
        return jsonify(permissions), 200
    except Exception as e:
        return jsonify({'message': str(e)}), 500


@rbac_bp.route('/me/has-role/<role_name>', methods=['GET'])
@jwt_required
def check_has_role(role_name):
    """Check if the current user has a specific role."""
    try:
        user = get_current_user()
        if not user:
            return jsonify({'message': 'User not authenticated'}), 401
        
        has_role = UserRoleManager.has_role(user['id'], role_name)
        return jsonify({'has_role': has_role}), 200
    except Exception as e:
        return jsonify({'message': str(e)}), 500


@rbac_bp.route('/me/has-permission/<resource_type>/<action>', methods=['GET'])
@jwt_required
def check_has_permission(resource_type, action):
    """Check if the current user has a specific permission."""
    try:
        user = get_current_user()
        if not user:
            return jsonify({'message': 'User not authenticated'}), 401
        
        has_permission = UserRoleManager.has_permission(user['id'], resource_type, action)
        return jsonify({'has_permission': has_permission}), 200
    except Exception as e:
        return jsonify({'message': str(e)}), 500