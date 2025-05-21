# Roo PM Lab Verification Completion Prompt

## Current Implementation Status

The lab verification workflow implementation is partially complete. Here's what has been done and what still needs to be completed:

### Completed Components:
1. ✅ **Database Schema**: The `lab_verifications` table schema exists in `migrations/013_lab_verification_schema.sql`
2. ✅ **Data Model**: The `LabVerification` class is implemented in `api/models.py` (lines 2604-2740)
3. ✅ **API Resources**: The `LabVerificationResource` and `VerificationStatsResource` classes exist in `api/lab_verification_resources.py`

### Missing Components:
1. ❌ **Schema Definitions**: The `lab_verification_fields` and `verification_stats_fields` are missing in `api/schemas.py`
2. ❌ **API Route Registration**: The lab verification endpoints need to be registered in `app.py`
3. ❌ **Frontend Implementation**: The `verification.js` file needs to be created
4. ❌ **Tests**: No tests have been created for lab verification functionality
5. ❌ **VerificationStatsResource Implementation**: The implementation detail for getting statistics in `VerificationStatsResource.get()` is missing

## Your Task

Your role is to delegate and coordinate the completion of the lab verification workflow implementation. Focus on the missing components listed above to ensure full functionality of the verification system.

## Detailed Implementation Guidelines

### 1. Schema Definitions (api/schemas.py)

Add the following field definitions to `api/schemas.py`:

```python
# Lab verification fields
lab_verification_fields = {
    'id': fields.String,
    'experiment_id': fields.String,
    'verification_status': fields.String,
    'verifier': fields.String,
    'equipment_used': fields.String,
    'comments': fields.String,
    'created_at': FlexibleDateTime(dt_format='iso8601'),
    'updated_at': FlexibleDateTime(dt_format='iso8601')
}

verification_stats_fields = {
    'total_count': fields.Integer,
    'verified_count': fields.Integer,
    'pending_count': fields.Integer,
    'rejected_count': fields.Integer,
    'verification_rate': fields.Float,
    'by_equipment': fields.Raw,
    'by_verifier': fields.Raw
}

# Lab verification schemas
class LabVerificationSchema(Schema):
    """
    Schema for lab verification data.
    
    This schema defines the structure for validation data associated with 
    experiments, tracking verification status and metadata.
    """
    id = fields.UUID(
        description="Unique identifier for the verification record",
        example="623e4567-e89b-12d3-a456-426614174005"
    )
    experiment_id = fields.UUID(
        required=True,
        description="ID of the experiment being verified",
        example="523e4567-e89b-12d3-a456-426614174004"
    )
    verification_status = fields.String(
        required=True,
        description="Status of verification",
        example="verified",
        validate=validate.OneOf(["pending", "verified", "rejected"])
    )
    verifier = fields.String(
        required=True,
        description="User ID or name of the person performing verification",
        example="a23e4567-e89b-12d3-a456-426614174000"
    )
    equipment_used = fields.String(
        required=True,
        description="Equipment used for verification",
        example="Differential Scanning Calorimeter Model XYZ-123"
    )
    comments = fields.String(
        description="Additional notes or comments about verification",
        example="Results match prediction within 2% margin of error"
    )
    created_at = fields.DateTime(
        description="Timestamp of when the verification was recorded",
        example="2023-04-05T14:20:00Z"
    )
    updated_at = fields.DateTime(
        description="Timestamp of when the verification was last updated",
        example="2023-04-06T09:15:00Z"
    )
```

### 2. API Route Registration (app.py)

Update `app.py` to register the lab verification endpoints:

```python
# Add to imports
from api.lab_verification_resources import LabVerificationResource, VerificationStatsResource

# Add to API resource registration (in the appropriate section)
api.add_resource(LabVerificationResource, '/api/experiments/<string:experiment_id>/verification', endpoint='experiment_verification')
api.add_resource(LabVerificationResource, '/api/verifications/<string:verification_id>', endpoint='verification_update')
api.add_resource(VerificationStatsResource, '/api/verification/stats', endpoint='verification_stats')
```

### 3. Frontend Implementation (static/js/verification.js)

Create a new file `static/js/verification.js` with the following content:

```javascript
/**
 * Lab Verification Module
 * 
 * This module handles the user interface for lab verification workflow,
 * including recording verification results, updating verification status,
 * and visualizing verification statistics.
 */

const VerificationModule = (function() {
    // Cache DOM elements
    let $verificationForm;
    let $statusSelect;
    let $equipmentInput;
    let $commentsInput;
    let $verificationList;
    let $statsContainer;
    
    // API endpoints
    const API_ENDPOINTS = {
        getVerification: (experimentId) => `/api/experiments/${experimentId}/verification`,
        updateVerification: (verificationId) => `/api/verifications/${verificationId}`,
        getStats: () => '/api/verification/stats'
    };
    
    // Verification status constants
    const STATUS = {
        PENDING: 'pending',
        VERIFIED: 'verified',
        REJECTED: 'rejected'
    };
    
    // Initialize the module
    function init() {
        // Initialize DOM elements
        cacheDomElements();
        bindEvents();
        
        // Load verification stats on initial load
        loadVerificationStats();
    }
    
    // Cache DOM elements for later use
    function cacheDomElements() {
        $verificationForm = $('#verification-form');
        $statusSelect = $('#verification-status');
        $equipmentInput = $('#verification-equipment');
        $commentsInput = $('#verification-comments');
        $verificationList = $('#verification-list');
        $statsContainer = $('#verification-stats');
    }
    
    // Bind event handlers
    function bindEvents() {
        // Form submission for recording verification
        $verificationForm.on('submit', handleVerificationSubmit);
        
        // Status update in verification list
        $verificationList.on('change', '.status-select', handleStatusUpdate);
        
        // Refresh stats button
        $('#refresh-stats').on('click', loadVerificationStats);
    }
    
    // Handle verification form submission
    function handleVerificationSubmit(e) {
        e.preventDefault();
        
        const experimentId = $verificationForm.data('experiment-id');
        const data = {
            verification_status: $statusSelect.val(),
            verifier: getCurrentUserId(),
            equipment_used: $equipmentInput.val(),
            comments: $commentsInput.val()
        };
        
        // Submit verification via API
        $.ajax({
            url: API_ENDPOINTS.getVerification(experimentId),
            method: 'POST',
            data: JSON.stringify(data),
            contentType: 'application/json',
            success: function(response) {
                showNotification('Verification recorded successfully', 'success');
                $verificationForm[0].reset();
                
                // Refresh verification list and stats
                loadExperimentVerification(experimentId);
                loadVerificationStats();
            },
            error: function(err) {
                showNotification('Error recording verification: ' + (err.responseJSON?.message || 'Unknown error'), 'error');
            }
        });
    }
    
    // Handle status update from verification list
    function handleStatusUpdate(e) {
        const $select = $(e.target);
        const verificationId = $select.data('verification-id');
        const newStatus = $select.val();
        const comments = prompt('Add comments for this status update (optional):');
        
        // Update verification status via API
        $.ajax({
            url: API_ENDPOINTS.updateVerification(verificationId),
            method: 'PUT',
            data: JSON.stringify({
                verification_status: newStatus,
                comments: comments
            }),
            contentType: 'application/json',
            success: function(response) {
                showNotification('Verification status updated successfully', 'success');
                
                // Refresh stats
                loadVerificationStats();
            },
            error: function(err) {
                showNotification('Error updating verification status: ' + (err.responseJSON?.message || 'Unknown error'), 'error');
                // Reset select to previous value
                $select.val($select.data('original-value'));
            }
        });
    }
    
    // Load verification for a specific experiment
    function loadExperimentVerification(experimentId) {
        $.ajax({
            url: API_ENDPOINTS.getVerification(experimentId),
            method: 'GET',
            success: function(verification) {
                if (verification) {
                    renderVerification(verification, experimentId);
                } else {
                    $verificationList.html('<p>No verification recorded for this experiment.</p>');
                }
            },
            error: function(err) {
                showNotification('Error loading verification: ' + (err.responseJSON?.message || 'Unknown error'), 'error');
            }
        });
    }
    
    // Load verification statistics
    function loadVerificationStats() {
        $.ajax({
            url: API_ENDPOINTS.getStats(),
            method: 'GET',
            success: function(stats) {
                renderStats(stats);
            },
            error: function(err) {
                showNotification('Error loading verification statistics: ' + (err.responseJSON?.message || 'Unknown error'), 'error');
            }
        });
    }
    
    // Render verification details
    function renderVerification(verification, experimentId) {
        const statusClass = {
            [STATUS.PENDING]: 'bg-warning',
            [STATUS.VERIFIED]: 'bg-success',
            [STATUS.REJECTED]: 'bg-danger'
        };
        
        const html = `
            <div class="card mt-3">
                <div class="card-header ${statusClass[verification.verification_status] || ''} text-white">
                    Verification Status: ${verification.verification_status.toUpperCase()}
                </div>
                <div class="card-body">
                    <h5 class="card-title">Verified by: ${verification.verifier}</h5>
                    <p><strong>Equipment Used:</strong> ${verification.equipment_used}</p>
                    <p><strong>Date:</strong> ${new Date(verification.created_at).toLocaleString()}</p>
                    ${verification.comments ? `<p><strong>Comments:</strong> ${verification.comments}</p>` : ''}
                    
                    <div class="form-group mt-3">
                        <label for="status-update-${verification.id}">Update Status:</label>
                        <select id="status-update-${verification.id}" 
                                class="form-control status-select" 
                                data-verification-id="${verification.id}"
                                data-original-value="${verification.verification_status}">
                            <option value="${STATUS.PENDING}" ${verification.verification_status === STATUS.PENDING ? 'selected' : ''}>Pending</option>
                            <option value="${STATUS.VERIFIED}" ${verification.verification_status === STATUS.VERIFIED ? 'selected' : ''}>Verified</option>
                            <option value="${STATUS.REJECTED}" ${verification.verification_status === STATUS.REJECTED ? 'selected' : ''}>Rejected</option>
                        </select>
                    </div>
                </div>
            </div>
        `;
        
        $verificationList.html(html);
    }
    
    // Render verification statistics
    function renderStats(stats) {
        if (!stats) {
            $statsContainer.html('<p>No verification statistics available.</p>');
            return;
        }
        
        const html = `
            <div class="card mt-3">
                <div class="card-header bg-primary text-white">
                    Verification Statistics
                </div>
                <div class="card-body">
                    <div class="row">
                        <div class="col-md-3">
                            <div class="stats-box">
                                <h5>Total Verifications</h5>
                                <div class="stats-value">${stats.total_count}</div>
                            </div>
                        </div>
                        <div class="col-md-3">
                            <div class="stats-box bg-success text-white">
                                <h5>Verified</h5>
                                <div class="stats-value">${stats.verified_count}</div>
                            </div>
                        </div>
                        <div class="col-md-3">
                            <div class="stats-box bg-warning text-white">
                                <h5>Pending</h5>
                                <div class="stats-value">${stats.pending_count}</div>
                            </div>
                        </div>
                        <div class="col-md-3">
                            <div class="stats-box bg-danger text-white">
                                <h5>Rejected</h5>
                                <div class="stats-value">${stats.rejected_count}</div>
                            </div>
                        </div>
                    </div>
                    
                    <div class="mt-4">
                        <h5>Verification Rate</h5>
                        <div class="progress">
                            <div class="progress-bar bg-success" role="progressbar" 
                                style="width: ${stats.verification_rate * 100}%" 
                                aria-valuenow="${stats.verification_rate * 100}" 
                                aria-valuemin="0" 
                                aria-valuemax="100">
                                ${Math.round(stats.verification_rate * 100)}%
                            </div>
                        </div>
                    </div>
                </div>
            </div>
        `;
        
        $statsContainer.html(html);
    }
    
    // Helper function to get current user ID
    function getCurrentUserId() {
        return window.currentUser?.id || 'unknown_user';
    }
    
    // Helper function to show notifications
    function showNotification(message, type) {
        // Use application's notification system
        if (window.app && window.app.showNotification) {
            window.app.showNotification(message, type);
        } else {
            alert(message);
        }
    }
    
    // Public API
    return {
        init: init,
        loadExperimentVerification: loadExperimentVerification,
        loadVerificationStats: loadVerificationStats,
        // Export constants
        STATUS: STATUS
    };
})();

// Initialize on document ready
$(document).ready(function() {
    VerificationModule.init();
});
```

Also add to `/mnt/c/Users/1edwa/Documents/CryoProtect v2/templates/experiments.html` to include the verification UI components.

### 4. VerificationStatsResource Implementation

Complete the implementation of `VerificationStatsResource.get()` method in `api/lab_verification_resources.py`:

```python
@token_required
@marshal_with(verification_stats_fields)
def get(self):
    """Get verification statistics."""
    try:
        supabase = get_supabase_client()
        
        # Get total counts by status
        stats = {
            'total_count': 0,
            'verified_count': 0,
            'pending_count': 0,
            'rejected_count': 0,
            'verification_rate': 0.0,
            'by_equipment': {},
            'by_verifier': {}
        }
        
        # Get all verifications
        response = supabase.table('lab_verifications').select('*').execute()
        if response.error:
            raise Exception(f"Error fetching verifications: {response.error}")
            
        verifications = response.data or []
        
        # Calculate counts
        stats['total_count'] = len(verifications)
        
        if stats['total_count'] > 0:
            # Count by status
            status_counts = {}
            equipment_counts = {}
            verifier_counts = {}
            
            for v in verifications:
                status = v.get('verification_status')
                equipment = v.get('equipment_used')
                verifier = v.get('verifier')
                
                # Count by status
                status_counts[status] = status_counts.get(status, 0) + 1
                
                # Count by equipment
                equipment_counts[equipment] = equipment_counts.get(equipment, 0) + 1
                
                # Count by verifier
                verifier_counts[verifier] = verifier_counts.get(verifier, 0) + 1
            
            # Set status counts
            stats['verified_count'] = status_counts.get('verified', 0)
            stats['pending_count'] = status_counts.get('pending', 0)
            stats['rejected_count'] = status_counts.get('rejected', 0)
            
            # Calculate verification rate (verified / total)
            stats['verification_rate'] = stats['verified_count'] / stats['total_count']
            
            # Set equipment and verifier stats
            stats['by_equipment'] = equipment_counts
            stats['by_verifier'] = verifier_counts
        
        return stats
    except Exception as e:
        return handle_error(e)
```

### 5. Tests for Lab Verification (tests/test_lab_verification.py)

Create a new file `tests/test_lab_verification.py` with the following content:

```python
import unittest
import uuid
from datetime import datetime, date
from unittest.mock import patch, MagicMock

from api.models import LabVerification
from api.lab_verification_resources import LabVerificationResource, VerificationStatsResource

class TestLabVerification(unittest.TestCase):
    """Test cases for lab verification model and API resources."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.experiment_id = str(uuid.uuid4())
        self.verification_id = str(uuid.uuid4())
        self.verifier = str(uuid.uuid4())
        self.equipment = "Test Equipment XYZ"
        self.mock_verification = {
            "id": self.verification_id,
            "experiment_id": self.experiment_id,
            "verification_status": "verified",
            "verifier": self.verifier,
            "equipment_used": self.equipment,
            "comments": "Test verification",
            "created_at": datetime.now().isoformat(),
            "updated_at": datetime.now().isoformat()
        }
    
    @patch('api.models.get_supabase_client')
    def test_record_verification(self, mock_get_supabase):
        """Test recording a new verification."""
        # Setup mock
        mock_supabase = MagicMock()
        mock_table = MagicMock()
        mock_execute = MagicMock()
        mock_execute.data = [self.mock_verification]
        mock_table.insert.return_value.execute.return_value = mock_execute
        mock_supabase.table.return_value = mock_table
        mock_get_supabase.return_value = mock_supabase
        
        # Test
        result = LabVerification.record_verification(
            experiment_id=self.experiment_id,
            verification_status="verified",
            verifier=self.verifier,
            equipment_used=self.equipment,
            comments="Test verification"
        )
        
        # Assert
        self.assertEqual(result, self.mock_verification)
        mock_supabase.table.assert_called_once_with('lab_verifications')
        mock_table.insert.assert_called_once()
    
    @patch('api.models.get_supabase_client')
    def test_get_verification(self, mock_get_supabase):
        """Test getting a verification by experiment ID."""
        # Setup mock
        mock_supabase = MagicMock()
        mock_table = MagicMock()
        mock_execute = MagicMock()
        mock_execute.data = [self.mock_verification]
        mock_table.select.return_value.eq.return_value.execute.return_value = mock_execute
        mock_supabase.table.return_value = mock_table
        mock_get_supabase.return_value = mock_supabase
        
        # Test
        result = LabVerification.get_verification(self.experiment_id)
        
        # Assert
        self.assertEqual(result, self.mock_verification)
        mock_supabase.table.assert_called_once_with('lab_verifications')
        mock_table.select.assert_called_once_with('*')
        mock_table.select.return_value.eq.assert_called_once_with('experiment_id', self.experiment_id)
    
    @patch('api.models.get_supabase_client')
    def test_update_verification_status(self, mock_get_supabase):
        """Test updating a verification status."""
        # Setup mock
        mock_supabase = MagicMock()
        mock_table = MagicMock()
        mock_execute = MagicMock()
        updated_verification = dict(self.mock_verification)
        updated_verification["verification_status"] = "rejected"
        updated_verification["comments"] = "Updated comments"
        mock_execute.data = [updated_verification]
        mock_table.update.return_value.eq.return_value.execute.return_value = mock_execute
        mock_supabase.table.return_value = mock_table
        mock_get_supabase.return_value = mock_supabase
        
        # Test
        result = LabVerification.update_verification_status(
            verification_id=self.verification_id,
            new_status="rejected",
            comments="Updated comments"
        )
        
        # Assert
        self.assertEqual(result, updated_verification)
        mock_supabase.table.assert_called_once_with('lab_verifications')
        mock_table.update.assert_called_once()
        mock_table.update.return_value.eq.assert_called_once_with('id', self.verification_id)
    
    @patch('api.lab_verification_resources.LabVerification')
    def test_lab_verification_resource_get(self, mock_lab_verification):
        """Test LabVerificationResource GET method."""
        # Setup mock
        mock_lab_verification.get_verification.return_value = self.mock_verification
        
        # Create resource
        resource = LabVerificationResource()
        
        # Test
        result = resource.get(self.experiment_id)
        
        # Assert
        self.assertEqual(result, self.mock_verification)
        mock_lab_verification.get_verification.assert_called_once_with(self.experiment_id)
    
    @patch('api.models.get_supabase_client')
    def test_verification_stats_resource_get(self, mock_get_supabase):
        """Test VerificationStatsResource GET method."""
        # Setup mock with multiple verifications
        mock_supabase = MagicMock()
        mock_table = MagicMock()
        mock_execute = MagicMock()
        verifications = [
            {
                "id": str(uuid.uuid4()),
                "experiment_id": str(uuid.uuid4()),
                "verification_status": "verified",
                "verifier": self.verifier,
                "equipment_used": self.equipment,
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
                "verifier": self.verifier,
                "equipment_used": self.equipment,
                "comments": "Verification 3",
                "created_at": datetime.now().isoformat(),
                "updated_at": datetime.now().isoformat()
            }
        ]
        mock_execute.data = verifications
        mock_table.select.return_value.execute.return_value = mock_execute
        mock_supabase.table.return_value = mock_table
        mock_get_supabase.return_value = mock_supabase
        
        # Create resource
        resource = VerificationStatsResource()
        
        # Test
        with patch('api.lab_verification_resources.get_supabase_client', return_value=mock_supabase):
            result = resource.get()
        
        # Assert
        self.assertEqual(result['total_count'], 3)
        self.assertEqual(result['verified_count'], 1)
        self.assertEqual(result['pending_count'], 1)
        self.assertEqual(result['rejected_count'], 1)
        self.assertAlmostEqual(result['verification_rate'], 1/3)
        self.assertEqual(len(result['by_equipment']), 2)
        self.assertEqual(len(result['by_verifier']), 2)
        mock_supabase.table.assert_called_with('lab_verifications')
        mock_table.select.assert_called_once_with('*')

if __name__ == '__main__':
    unittest.main()
```

## Task Prioritization

1. **Schema Definitions**: Implement first as other components depend on these field definitions
2. **API Route Registration**: Register the endpoints in app.py
3. **VerificationStatsResource Implementation**: Complete the implementation for statistics
4. **Frontend Implementation**: Create the verification.js file
5. **Tests**: Implement tests to validate all functionality

## Coordination and Delegation Strategy

1. **Backend Engineer**: Assign the schema definitions, API route registration, and VerificationStatsResource implementation
2. **Frontend Developer**: Assign the verification.js implementation and template updates
3. **QA Engineer**: Assign the test implementation

Have team members work in parallel on these tasks with coordination points to ensure all components integrate correctly. Schedule a brief meeting after schema definitions are completed to ensure everyone has a shared understanding of the data structure.

## Testing and Verification Plan

1. **Unit Tests**: Validate all model methods and API resources
2. **Integration Test**: Verify full verification workflow from API to database
3. **Frontend Test**: Manually test the UI implementation
4. **End-to-end Test**: Verify complete workflow from UI to database and back

## Documentation Updates

After implementation, ensure the following documentation updates:

1. Update API documentation to include verification endpoints
2. Add user guide for the verification workflow
3. Update developer documentation with the new components

## Next Steps After Completion

Once the lab verification workflow is fully implemented, proceed to:

1. The API Architecture Standardization task
2. Testing Framework Completion
3. Any critical bug fixes identified during implementation