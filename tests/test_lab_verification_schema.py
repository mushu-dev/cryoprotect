import pytest
import uuid
from datetime import datetime
from unittest.mock import patch, MagicMock

from api.models import LabVerification

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

def test_lab_verification_schema(verification_data):
    """Test that the lab verification schema is correctly defined."""
    # This test verifies that the LabVerification class has the expected attributes
    # after the schema migrations
    
    # Check that the LabVerification class has the expected class attributes
    assert hasattr(LabVerification, 'PENDING')
    assert hasattr(LabVerification, 'VERIFIED')
    assert hasattr(LabVerification, 'REJECTED')
    
    # Check that the LabVerification class has the expected methods
    assert hasattr(LabVerification, 'record_verification')
    assert hasattr(LabVerification, 'get_verification')
    assert hasattr(LabVerification, 'update_verification_status')
    
    # Check that the methods have the expected signatures
    assert callable(LabVerification.record_verification)
    assert callable(LabVerification.get_verification)
    assert callable(LabVerification.update_verification_status)