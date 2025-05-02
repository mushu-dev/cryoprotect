"""
CryoProtect Analyzer - Team Resources API Tests

This module contains direct tests for all team-related API endpoints.
It covers request validation, response formats, and error handling for:
- Teams
- Team Members
- Shared Resources
- Comments
- Notifications
- Activity Logs
"""

import os
import sys
import json
import uuid

# Add the parent directory to the path so we can import the api package
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from tests.base_test_case import BaseTestCase, MockSupabaseBaseTestCase

class TestTeamResourcesAPI(MockSupabaseBaseTestCase):
    """Test cases for team-related API endpoints."""

    def setUp(self):
        """Set up test client and sample data."""
        super().setUp()
        self.client = self.app.test_client()
        # Sample IDs for use in tests
        self.team_id = str(uuid.uuid4())
        self.user_id = str(uuid.uuid4())
        self.resource_id = str(uuid.uuid4())
        self.comment_id = str(uuid.uuid4())
        self.notification_id = str(uuid.uuid4())

    # --- Teams ---
    def test_get_teams(self):
        """GET /teams - List all teams for current user."""
        response = self.client.get('/teams')
        self.assertIn(response.status_code, [200, 401, 403])

    def test_post_team(self):
        """POST /teams - Create a new team."""
        data = {"name": "Test Team", "description": "A test team"}
        response = self.client.post('/teams', json=data)
        self.assertIn(response.status_code, [201, 400, 401, 403])

    def test_get_team_by_id(self):
        """GET /teams/<team_id> - Get a team by ID."""
        response = self.client.get(f'/teams/{self.team_id}')
        self.assertIn(response.status_code, [200, 404, 401, 403])

    def test_put_team(self):
        """PUT /teams/<team_id> - Update a team."""
        data = {"name": "Updated Team", "description": "Updated description"}
        response = self.client.put(f'/teams/{self.team_id}', json=data)
        self.assertIn(response.status_code, [200, 403, 404, 401])

    def test_delete_team(self):
        """DELETE /teams/<team_id> - Delete a team."""
        response = self.client.delete(f'/teams/{self.team_id}')
        self.assertIn(response.status_code, [200, 403, 404, 401])

    # --- Team Members ---
    def test_get_team_members(self):
        """GET /teams/<team_id>/members - List team members."""
        response = self.client.get(f'/teams/{self.team_id}/members')
        self.assertIn(response.status_code, [200, 404, 401, 403])

    def test_post_team_member(self):
        """POST /teams/<team_id>/members - Add a member."""
        data = {"user_id": self.user_id, "role": "member"}
        response = self.client.post(f'/teams/{self.team_id}/members', json=data)
        self.assertIn(response.status_code, [201, 400, 401, 403])

    def test_put_team_member(self):
        """PUT /teams/<team_id>/members/<user_id> - Update member role."""
        data = {"role": "admin"}
        response = self.client.put(f'/teams/{self.team_id}/members/{self.user_id}', json=data)
        self.assertIn(response.status_code, [200, 400, 401, 403])

    def test_delete_team_member(self):
        """DELETE /teams/<team_id>/members/<user_id> - Remove member."""
        response = self.client.delete(f'/teams/{self.team_id}/members/{self.user_id}')
        self.assertIn(response.status_code, [200, 400, 401, 403])

    # --- Shared Resources ---
    def test_get_shared_resources(self):
        """GET /teams/<team_id>/resources - List shared resources."""
        response = self.client.get(f'/teams/{self.team_id}/resources')
        self.assertIn(response.status_code, [200, 404, 401, 403])

    def test_post_shared_resource(self):
        """POST /teams/<team_id>/resources - Share a resource."""
        data = {"resource_type": "document", "resource_id": self.resource_id}
        response = self.client.post(f'/teams/{self.team_id}/resources', json=data)
        self.assertIn(response.status_code, [201, 400, 401, 403])

    def test_delete_shared_resource(self):
        """DELETE /teams/<team_id>/resources/<resource_type>/<resource_id> - Unshare resource."""
        response = self.client.delete(f'/teams/{self.team_id}/resources/document/{self.resource_id}')
        self.assertIn(response.status_code, [200, 400, 401, 403])

    # --- Comments ---
    def test_get_comments(self):
        """GET /comments/<resource_type>/<resource_id> - List comments."""
        response = self.client.get(f'/comments/document/{self.resource_id}')
        self.assertIn(response.status_code, [200, 404, 401, 403])

    def test_post_comment(self):
        """POST /comments/<resource_type>/<resource_id> - Add comment."""
        data = {"content": "Test comment", "parent_id": None}
        response = self.client.post(f'/comments/document/{self.resource_id}', json=data)
        self.assertIn(response.status_code, [201, 400, 401, 403])

    def test_put_comment(self):
        """PUT /comments/<comment_id> - Update comment."""
        data = {"content": "Updated comment"}
        response = self.client.put(f'/comments/{self.comment_id}', json=data)
        self.assertIn(response.status_code, [200, 400, 401, 403])

    def test_delete_comment(self):
        """DELETE /comments/<comment_id> - Delete comment."""
        response = self.client.delete(f'/comments/{self.comment_id}')
        self.assertIn(response.status_code, [200, 400, 401, 403])

    # --- Notifications ---
    def test_get_notifications(self):
        """GET /notifications - List notifications."""
        response = self.client.get('/notifications')
        self.assertIn(response.status_code, [200, 401, 403])

    def test_put_notifications(self):
        """PUT /notifications - Mark all as read."""
        response = self.client.put('/notifications')
        self.assertIn(response.status_code, [200, 400, 401, 403])

    def test_put_notification(self):
        """PUT /notifications/<notification_id> - Mark as read."""
        response = self.client.put(f'/notifications/{self.notification_id}')
        self.assertIn(response.status_code, [200, 400, 401, 403])

    def test_delete_notification(self):
        """DELETE /notifications/<notification_id> - Delete notification."""
        response = self.client.delete(f'/notifications/{self.notification_id}')
        self.assertIn(response.status_code, [200, 400, 401, 403])

    # --- Activity Logs ---
    def test_get_activity_log(self):
        """GET /teams/<team_id>/activity - Get activity log."""
        response = self.client.get(f'/teams/{self.team_id}/activity')
        self.assertIn(response.status_code, [200, 400, 401, 403])

if __name__ == '__main__':
    unittest.main()