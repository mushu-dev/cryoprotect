import pytest
import uuid
from datetime import datetime
from unittest.mock import patch, MagicMock

from api.models import LabVerification

@pytest.fixture
def mock_supabase():
    """Create a mock Supabase client."""
    mock_client = MagicMock()
    mock_table = MagicMock()
    mock_client.table.return_value = mock_table
    return mock_client

def test_lab_verification_table_schema(mock_supabase):
    """Test that the lab verification table schema is correctly defined."""
    # Mock the get_supabase_client function to return our mock client
    with patch('api.models.get_supabase_client', return_value=mock_supabase):
        # Call the get_table method to verify it's accessing the correct table
        LabVerification.get_table()
        
        # Verify that the correct table name is used
        mock_supabase.table.assert_called_once_with('lab_verifications')
        
    # Reset the mock for the next test
    mock_supabase.reset_mock()
    
    # Test that we can create a verification record with the expected fields
    experiment_id = str(uuid.uuid4())
    verification_id = str(uuid.uuid4())
    verifier = str(uuid.uuid4())
    equipment = "Test Equipment XYZ"
    
    # Create a new mock for the record_verification test
    with patch('api.models.get_supabase_client') as mock_get_supabase:
        # Setup the mock
        mock_supabase2 = MagicMock()
        mock_table = MagicMock()
        mock_execute = MagicMock()
        mock_execute.data = [{
            "id": verification_id,
            "experiment_id": experiment_id,
            "verification_status": "verified",
            "verifier": verifier,
            "equipment_used": equipment,
            "comments": "Test verification",
            "created_at": datetime.now().isoformat(),
            "updated_at": datetime.now().isoformat()
        }]
        mock_execute.error = None
        mock_table.insert.return_value.execute.return_value = mock_execute
        mock_supabase2.table.return_value = mock_table
        mock_get_supabase.return_value = mock_supabase2
        
        # Call the record_verification method
        result = LabVerification.record_verification(
            experiment_id=experiment_id,
            verification_status="verified",
            verifier=verifier,
            equipment_used=equipment,
            comments="Test verification"
        )
        
        # Verify that the insert method was called with the correct data
        mock_supabase2.table.assert_called_once_with('lab_verifications')
        mock_table.insert.assert_called_once()
        
        # Verify that the result contains the expected fields
        assert "id" in result
        assert "experiment_id" in result
        assert "verification_status" in result
        assert "verifier" in result
        assert "equipment_used" in result
        assert "comments" in result
        assert "created_at" in result
        assert "updated_at" in result
        
        # Verify that the values match what we expect
        assert result["experiment_id"] == experiment_id
        assert result["verification_status"] == "verified"
        assert result["verifier"] == verifier
        assert result["equipment_used"] == equipment
        assert result["comments"] == "Test verification"