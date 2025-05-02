# Phase 2.1: API Layer Completion - Initial Task Delegation

## Task Overview: API Response Format Standardization

This document provides a detailed implementation plan for standardizing the API response format across all endpoints in the CryoProtect v2 project. This is the first subtask in Phase 2.1 (API Layer Completion).

## Files to Examine

```
api/
├── __init__.py
├── models.py
├── resources.py
├── schemas.py
├── utils.py
├── batch_resources.py
├── comparisons.py
├── export_resources.py
├── mixture_analysis_resources.py
├── predictive_models_resources.py
├── protocol_designer_resources.py
├── rdkit_resources.py
├── scoring_resources.py
└── team_resources.py
```

## Implementation Steps

### Step 1: Create Standard Response Envelope

Create a unified response format in `api/utils.py` that all endpoints will use:

```python
def create_standard_response(data=None, message=None, status="success", meta=None):
    """
    Create a standardized API response.
    
    Args:
        data: The main response data
        message: Optional message about the response
        status: Status of the response ("success", "error", "warning")
        meta: Optional metadata (pagination, etc.)
        
    Returns:
        dict: Standardized response dictionary
    """
    response = {
        "status": status,
        "data": data or {},
    }
    
    if message:
        response["message"] = message
        
    if meta:
        response["meta"] = meta
        
    return response


def create_error_response(message, error_code=None, details=None):
    """
    Create a standardized error response.
    
    Args:
        message: Error message
        error_code: Optional error code identifier
        details: Optional detailed error information
        
    Returns:
        dict: Standardized error response
    """
    error_data = {
        "message": message
    }
    
    if error_code:
        error_data["code"] = error_code
        
    if details:
        error_data["details"] = details
    
    return create_standard_response(
        data=error_data,
        status="error",
        message=message
    )
```

### Step 2: Implement Schema Validation

Enhance `api/schemas.py` with Marshmallow schemas for response validation:

```python
from marshmallow import Schema, fields

class MetaSchema(Schema):
    """Schema for response metadata"""
    page = fields.Int(required=False)
    per_page = fields.Int(required=False) 
    total = fields.Int(required=False)
    total_pages = fields.Int(required=False)

class StandardResponseSchema(Schema):
    """Schema for standardized API responses"""
    status = fields.Str(required=True)
    message = fields.Str(required=False)
    data = fields.Dict(required=True)
    meta = fields.Nested(MetaSchema, required=False)
```

### Step 3: Update Base Resource Class

Modify the base resource class in `api/resources.py` to use the standard response format:

```python
class BaseResource(Resource):
    """Base resource class with standardized response handling"""
    
    def respond(self, data=None, message=None, status_code=200, **kwargs):
        """
        Create a standardized Flask response
        
        Args:
            data: The response data
            message: Response message
            status_code: HTTP status code
            **kwargs: Additional parameters for create_standard_response
            
        Returns:
            Flask response with standardized format
        """
        response = create_standard_response(data=data, message=message, **kwargs)
        return response, status_code
    
    def respond_error(self, message, status_code=400, error_code=None, details=None):
        """
        Create a standardized error response
        
        Args:
            message: Error message
            status_code: HTTP status code
            error_code: Error code identifier
            details: Detailed error information
            
        Returns:
            Flask response with standardized error format
        """
        response = create_error_response(
            message=message, 
            error_code=error_code,
            details=details
        )
        return response, status_code
```

### Step 4: Update HTTP Status Codes Map

Add a status code mapping in `api/utils.py`:

```python
# HTTP Status Code Mapping
HTTP_STATUS_CODES = {
    # Success
    "ok": 200,
    "created": 201,
    "accepted": 202,
    "no_content": 204,
    
    # Client Errors
    "bad_request": 400,
    "unauthorized": 401,
    "forbidden": 403,
    "not_found": 404,
    "method_not_allowed": 405,
    "conflict": 409,
    "gone": 410,
    "precondition_failed": 412,
    "unsupported_media_type": 415,
    "too_many_requests": 429,
    
    # Server Errors
    "server_error": 500,
    "not_implemented": 501,
    "bad_gateway": 502,
    "service_unavailable": 503
}

def get_status_code(status_key):
    """Get HTTP status code by key"""
    return HTTP_STATUS_CODES.get(status_key, 200)
```

### Step 5: Update API Endpoints

Modify priority endpoints in `api/resources.py` to use the new standard format.

Here's a sample implementation for one endpoint:

```python
class MoleculeResource(BaseResource):
    def get(self, molecule_id):
        try:
            molecule = Molecule.query.get(molecule_id)
            
            if not molecule:
                return self.respond_error(
                    message=f"Molecule with ID {molecule_id} not found",
                    status_code=404,
                    error_code="molecule_not_found"
                )
            
            # Get the molecule data
            molecule_data = molecule.to_dict()
            
            return self.respond(
                data=molecule_data,
                message=f"Successfully retrieved molecule {molecule_id}",
                meta={"last_updated": molecule.updated_at.isoformat() if hasattr(molecule, 'updated_at') else None}
            )
            
        except Exception as e:
            return self.respond_error(
                message="Failed to retrieve molecule",
                status_code=500,
                error_code="server_error",
                details=str(e)
            )
```

### Step 6: Write Tests

Create tests for the standardized response format in `tests/test_api_standardization.py`:

```python
import unittest
import json
from app import create_app
from api.utils import create_standard_response, create_error_response

class TestApiStandardization(unittest.TestCase):
    def setUp(self):
        self.app = create_app('testing')
        self.client = self.app.test_client()
        self.app_context = self.app.app_context()
        self.app_context.push()
        
    def tearDown(self):
        self.app_context.pop()
        
    def test_standard_response_format(self):
        # Test standard response creation
        response = create_standard_response(data={"id": 1, "name": "Test"}, message="Test message")
        self.assertEqual(response["status"], "success")
        self.assertEqual(response["message"], "Test message")
        self.assertEqual(response["data"]["id"], 1)
        
    def test_error_response_format(self):
        # Test error response creation
        response = create_error_response(message="Error message", error_code="test_error")
        self.assertEqual(response["status"], "error")
        self.assertEqual(response["message"], "Error message")
        self.assertEqual(response["data"]["code"], "test_error")
        
    def test_api_endpoint_response(self):
        # Test an actual API endpoint response
        response = self.client.get('/api/molecules/1')
        data = json.loads(response.data)
        
        # Verify response structure
        self.assertIn("status", data)
        self.assertIn("data", data)
```

## Acceptance Criteria

1. All API utility functions for standardized responses are implemented
2. Response schemas are defined for validation
3. The base resource class uses the standardized format
4. Status code mapping is implemented
5. At least one priority endpoint is updated to use the new format
6. Tests verify the response format implementation
7. No regressions in existing functionality

## Next Steps After Completion

1. Extend implementation to all remaining endpoints
2. Create documentation for the standard response format
3. Add API response examples to the OpenAPI documentation
4. Implement proper error handling middleware

## Reporting Requirements

After implementation, provide:
1. Summary of changes made
2. Any challenges encountered
3. Test results
4. Recommendations for the next subtask