import pytest
import uuid
from datetime import datetime, date
from unittest.mock import patch, MagicMock
from flask import Flask, request

from api.models import LabVerification
from api.lab_verification_resources import LabVerificationResource, VerificationStatsResource

# Test fixtures
@pytest.fixture
def verification_data():
    """Create test verification data."""
    experiment_id = str(uuid.uuid4())
    verification_id = str(uuid.uuid4())
    verifier = str(uuid.uuid4())
    equipment = "Test Equipment XYZ"
    mock_verification = {
        "id": verification_id,
        "experiment_id": experiment_id,
        "verification_status": "verified",
        "verifier": verifier,
        "equipment_used": equipment,
        "comments": "Test verification",
        "created_at": datetime.now().isoformat(),
        "updated_at": datetime.now().isoformat()
    }
    
    return {
        "experiment_id": experiment_id,
        "verification_id": verification_id,
        "verifier": verifier,
        "equipment": equipment,
        "mock_verification": mock_verification
    }

@pytest.mark.usefixtures("app")
def test_record_verification(app, verification_data):
    """Test recording a new verification."""
    with app.app_context():
        # Setup mock
        with patch('api.models.get_supabase_client') as mock_get_supabase:
            mock_supabase = MagicMock()
            mock_table = MagicMock()
            mock_execute = MagicMock()
            mock_execute.data = [verification_data["mock_verification"]]
            # Set up the error attribute to be None to avoid the error handling path
            mock_execute.error = None
            mock_table.insert.return_value.execute.return_value = mock_execute
            mock_supabase.table.return_value = mock_table
            mock_get_supabase.return_value = mock_supabase
        
            # Test
            result = LabVerification.record_verification(
                experiment_id=verification_data["experiment_id"],
                verification_status="verified",
                verifier=verification_data["verifier"],
                equipment_used=verification_data["equipment"],
                comments="Test verification"
            )
            
            # Assert
            assert result == verification_data["mock_verification"]
            mock_supabase.table.assert_called_once_with('lab_verifications')
            mock_table.insert.assert_called_once()

@pytest.mark.usefixtures("app")
def test_get_verification(app, verification_data):
    """Test getting a verification by experiment ID."""
    with app.app_context():
        # Setup mock
        with patch('api.models.get_supabase_client') as mock_get_supabase:
            mock_supabase = MagicMock()
            mock_table = MagicMock()
            mock_execute = MagicMock()
            mock_execute.data = [verification_data["mock_verification"]]
            mock_table.select.return_value.eq.return_value.execute.return_value = mock_execute
            mock_supabase.table.return_value = mock_table
            mock_get_supabase.return_value = mock_supabase
        
        # Test
            # Set up the error attribute to be None to avoid the error handling path
            mock_execute.error = None
            
            # Test
            result = LabVerification.get_verification(verification_data["experiment_id"])
            
            # Assert
            assert result == verification_data["mock_verification"]
            mock_supabase.table.assert_called_once_with('lab_verifications')
            mock_table.select.assert_called_once_with('*')
            mock_table.select.return_value.eq.assert_called_once_with('experiment_id', verification_data["experiment_id"])

@pytest.mark.usefixtures("app")
def test_update_verification_status(app, verification_data):
    """Test updating a verification status."""
    with app.app_context():
        # Setup mock
        with patch('api.models.get_supabase_client') as mock_get_supabase:
            mock_supabase = MagicMock()
            mock_table = MagicMock()
            mock_execute = MagicMock()
            updated_verification = dict(verification_data["mock_verification"])
            updated_verification["verification_status"] = "rejected"
            updated_verification["comments"] = "Updated comments"
            mock_execute.data = [updated_verification]
            # Set up the error attribute to be None to avoid the error handling path
            mock_execute.error = None
            # Make sure the error attribute is properly mocked for the handle_supabase_error function
            mock_execute.error = None
            mock_table.update.return_value.eq.return_value.execute.return_value = mock_execute
            mock_supabase.table.return_value = mock_table
            mock_get_supabase.return_value = mock_supabase
            # Test
            result = LabVerification.update_verification_status(
                verification_id=verification_data["verification_id"],
                new_status="rejected",
                comments="Updated comments"
            )
            
            # Assert
            assert result == updated_verification
            mock_supabase.table.assert_called_once_with('lab_verifications')
            mock_table.update.assert_called_once()
            mock_table.update.return_value.eq.assert_called_once_with('id', verification_data["verification_id"])

@pytest.mark.usefixtures("app")
def test_lab_verification_resource_get(app, verification_data):
    """Test LabVerificationResource GET method."""
    with app.app_context():
        with app.test_request_context():
            # Mock the jwt_auth functions to bypass authentication
            with patch('api.jwt_auth.get_token_from_request', return_value='fake-token'):
                with patch('api.jwt_auth.extract_user_from_token', return_value=({'id': '123'}, '123')):
                    # Mock the marshal_with decorator to return the data directly
                    with patch('api.lab_verification_resources.marshal_with', lambda x: lambda f: f):
                        with patch('api.lab_verification_resources.LabVerification') as mock_lab_verification:
                            # Setup mock
                            mock_lab_verification.get_verification.return_value = verification_data["mock_verification"]
                            
                            # Create resource
                            resource = LabVerificationResource()
                            
                            # Test
                            result = resource.get(verification_data["experiment_id"])
                            
                            # Assert
                            assert result == verification_data["mock_verification"]
                            mock_lab_verification.get_verification.assert_called_once_with(verification_data["experiment_id"])
    
@pytest.mark.usefixtures("app")
def test_verification_stats_resource_get(app, verification_data):
    """Test VerificationStatsResource GET method."""
    with app.app_context():
        with app.test_request_context():
            # Mock the jwt_auth functions to bypass authentication
            with patch('api.jwt_auth.get_token_from_request', return_value='fake-token'):
                with patch('api.jwt_auth.extract_user_from_token', return_value=({'id': '123'}, '123')):
                    # Mock the marshal_with decorator to return the data directly
                    with patch('api.lab_verification_resources.marshal_with', lambda x: lambda f: f):
                        # Setup mock with multiple verifications
                        with patch('api.models.get_supabase_client') as mock_get_supabase:
                            mock_supabase = MagicMock()
                            mock_table = MagicMock()
                            mock_execute = MagicMock()
                            verifications = [
                                {
                                    "id": str(uuid.uuid4()),
                                    "experiment_id": str(uuid.uuid4()),
                                    "verification_status": "verified",
                                    "verifier": verification_data["verifier"],
                                    "equipment_used": verification_data["equipment"],
                                    "comments": "Verification 1",
                                    "created_at": datetime.now().isoformat(),
                                    "updated_at": datetime.now().isoformat()
                                },
                                {
                                    "id": str(uuid.uuid4()),
                                    "experiment_id": str(uuid.uuid4()),
                                    "verification_status": "pending",
                                    "verifier": str(uuid.uuid4()),
                                    "equipment_used": "Other Equipment",
                                    "comments": "Verification 2",
                                    "created_at": datetime.now().isoformat(),
                                    "updated_at": datetime.now().isoformat()
                                },
                                {
                                    "id": str(uuid.uuid4()),
                                    "experiment_id": str(uuid.uuid4()),
                                    "verification_status": "rejected",
                                    "verifier": verification_data["verifier"],
                                    "equipment_used": verification_data["equipment"],
                                    "comments": "Verification 3",
                                    "created_at": datetime.now().isoformat(),
                                    "updated_at": datetime.now().isoformat()
                                }
                            ]
                            mock_execute.data = verifications
                            # Set up the error attribute to be None to avoid the error handling path
                            mock_execute.error = None
                            mock_table.select.return_value.execute.return_value = mock_execute
                            mock_supabase.table.return_value = mock_table
                            mock_get_supabase.return_value = mock_supabase
                            
                            # Create resource
                            resource = VerificationStatsResource()
                            
                            # Test
                            with patch('api.lab_verification_resources.get_supabase_client', return_value=mock_supabase):
                                result = resource.get()
            
                                # Assert
                                assert result['total_count'] == 3
                                assert result['verified_count'] == 1
                                assert result['pending_count'] == 1
                                assert result['rejected_count'] == 1
                                assert result['verification_rate'] == pytest.approx(1/3)
                                assert len(result['by_equipment']) == 2
                                assert len(result['by_verifier']) == 2
                                mock_supabase.table.assert_called_with('lab_verifications')
                                mock_table.select.assert_called_once_with('*')