"""
CryoProtect Analyzer - Team Models & Resources Tests

This module contains unit tests for:
- Marshmallow schemas in api/team_models.py (TeamSchema, TeamMemberSchema, SharedResourceSchema, CommentSchema, NotificationSchema)
- (To be expanded) API resource endpoints in api/team_resources.py

Covers: serialization, required/optional fields, allowed values, and error handling.
"""

import uuid
from marshmallow import ValidationError
import sys
import os
import unittest

# Add parent directory to path to import from api
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from tests.base_test_case import BaseTestCase

from api.team_models import (
    TeamSchema,
    TeamMemberSchema,
    SharedResourceSchema,
    CommentSchema,
    NotificationSchema
)

class TestTeamSchema(BaseTestCase):
    def setUp(self):
        self.schema = TeamSchema()

    def test_valid_team(self):
        data = {"name": "Cryo Team", "description": "Research group"}
        result = self.schema.load(data)
        self.assertEqual(result["name"], "Cryo Team")
        self.assertEqual(result["description"], "Research group")

    def test_missing_name(self):
        data = {"description": "No name"}
        with self.assertRaises(ValidationError) as ctx:
            self.schema.load(data)
        self.assertIn("name", ctx.exception.messages)

class TestTeamMemberSchema(BaseTestCase):
    def setUp(self):
        self.schema = TeamMemberSchema()

    def test_valid_member(self):
        data = {
            "user_id": str(uuid.uuid4()),
            "role": "admin"
        }
        result = self.schema.load(data)
        self.assertEqual(result["role"], "admin")

    def test_invalid_role(self):
        data = {
            "user_id": str(uuid.uuid4()),
            "role": "invalid_role"
        }
        with self.assertRaises(ValidationError) as ctx:
            self.schema.load(data)
        self.assertIn("role", ctx.exception.messages)

    def test_missing_user_id(self):
        data = {"role": "viewer"}
        with self.assertRaises(ValidationError) as ctx:
            self.schema.load(data)
        self.assertIn("user_id", ctx.exception.messages)

class TestSharedResourceSchema(BaseTestCase):
    def setUp(self):
        self.schema = SharedResourceSchema()

    def test_valid_resource(self):
        data = {
            "team_id": str(uuid.uuid4()),
            "resource_type": "project",
            "resource_id": str(uuid.uuid4())
        }
        result = self.schema.load(data)
        self.assertEqual(result["resource_type"], "project")

    def test_invalid_resource_type(self):
        data = {
            "team_id": str(uuid.uuid4()),
            "resource_type": "invalid",
            "resource_id": str(uuid.uuid4())
        }
        with self.assertRaises(ValidationError) as ctx:
            self.schema.load(data)
        self.assertIn("resource_type", ctx.exception.messages)

    def test_missing_team_id(self):
        data = {
            "resource_type": "project",
            "resource_id": str(uuid.uuid4())
        }
        with self.assertRaises(ValidationError) as ctx:
            self.schema.load(data)
        self.assertIn("team_id", ctx.exception.messages)

class TestCommentSchema(BaseTestCase):
    def setUp(self):
        self.schema = CommentSchema()

    def test_valid_comment(self):
        data = {
            "resource_type": "mixture",
            "resource_id": str(uuid.uuid4()),
            "content": "Looks good!",
            "parent_id": None
        }
        result = self.schema.load(data)
        self.assertEqual(result["resource_type"], "mixture")
        self.assertEqual(result["content"], "Looks good!")

    def test_missing_content(self):
        data = {
            "resource_type": "experiment",
            "resource_id": str(uuid.uuid4())
        }
        with self.assertRaises(ValidationError) as ctx:
            self.schema.load(data)
        self.assertIn("content", ctx.exception.messages)

    def test_invalid_resource_type(self):
        data = {
            "resource_type": "invalid",
            "resource_id": str(uuid.uuid4()),
            "content": "Test"
        }
        with self.assertRaises(ValidationError) as ctx:
            self.schema.load(data)
        self.assertIn("resource_type", ctx.exception.messages)

class TestNotificationSchema(BaseTestCase):
    def setUp(self):
        self.schema = NotificationSchema()

    def test_valid_notification(self):
        data = {
            "user_id": str(uuid.uuid4()),
            "team_id": str(uuid.uuid4()),
            "title": "New Comment",
            "message": "A new comment was added.",
            "resource_type": "comment",
            "resource_id": str(uuid.uuid4())
        }
        result = self.schema.load(data)
        self.assertEqual(result["title"], "New Comment")
        self.assertEqual(result["resource_type"], "comment")

    def test_missing_user_id(self):
        data = {
            "title": "Missing User",
            "message": "No user_id"
        }
        with self.assertRaises(ValidationError) as ctx:
            self.schema.load(data)
        self.assertIn("user_id", ctx.exception.messages)

    def test_invalid_resource_type(self):
        data = {
            "user_id": str(uuid.uuid4()),
            "title": "Invalid Resource",
            "message": "Bad type",
            "resource_type": "badtype"
        }
        with self.assertRaises(ValidationError) as ctx:
            self.schema.load(data)
        self.assertIn("resource_type", ctx.exception.messages)

if __name__ == "__main__":
    unittest.main()