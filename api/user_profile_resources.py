"""
CryoProtect Analyzer API - User Profile Resource

Provides endpoints for managing the authenticated user's profile.
"""

from flask import request, jsonify
from api.utils import _handle_json_serialization
from flask_restful import Resource
from api.models import UserProfile
from api.utils import token_required, get_user_id

from api.schemas import UserProfileSchema

class UserProfileResource(Resource):
    method_decorators = [token_required]

    def get(self):
        """
        Get the current user's profile.
        """
        user_id = get_user_id()
        if not user_id:
            return jsonify(_handle_json_serialization({"error": "Not authenticated"})), 401
        profile = UserProfile.get_by_user_id(user_id)
        if not profile:
            return jsonify(_handle_json_serialization({"error": "Profile not found"})), 404
        return jsonify(_handle_json_serialization(profile)), 200

    def post(self):
        """
        Create or upsert the current user's profile.
        Expects JSON: { "email": ..., "name": ... }
        """
        user_id = get_user_id()
        if not user_id:
            return jsonify(_handle_json_serialization({"error": "Not authenticated"})), 401
        data = request.get_json(force=True)
        schema = UserProfileSchema()
        try:
            validated = schema.load(data)
        except Exception as err:
            return jsonify(_handle_json_serialization({"error": "Validation failed", "messages": err.messages})), 400
        email = validated.get("email")
        name = validated.get("name")
        profile = UserProfile.create_or_update(user_id, email, name)
        return jsonify(_handle_json_serialization(profile)), 201

    def put(self):
        """
        Update the current user's profile.
        Expects JSON: { "name": ..., ... }
        """
        user_id = get_user_id()
        if not user_id:
            return jsonify(_handle_json_serialization({"error": "Not authenticated"})), 401
        updates = request.get_json(force=True)
        # Prevent changing user_id/email via update
        updates.pop("user_id", None)
        updates.pop("email", None)
        schema = UserProfileSchema()
        try:
            validated = schema.load(updates, partial=True)
        except Exception as err:
            return jsonify(_handle_json_serialization({"error": "Validation failed", "messages": err.messages})), 400
        profile = UserProfile.create_or_update(user_id, None, validated.get("name"))
        return jsonify(_handle_json_serialization(profile)), 200
        profile = UserProfile.update_profile(user_id, updates)
        return jsonify(_handle_json_serialization(profile)), 200