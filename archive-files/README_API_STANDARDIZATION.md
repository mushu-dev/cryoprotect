# API Endpoint Standardization for CryoProtect v2

This document provides an overview of the API endpoint standardization implemented in CryoProtect v2. It explains the standardization approach, the changes made, and how to use the standardized API.

## Overview

The API endpoint standardization project aimed to:

1. Ensure consistent response formats across all endpoints
2. Implement proper HTTP status codes for all API responses
3. Standardize error handling across the API
4. Document all endpoints in OpenAPI format

## Standardized Response Format

All API endpoints now return responses in a consistent format:

```json
{
  "status": "success",
  "timestamp": "2025-04-21T17:45:00.000Z",
  "code": 200,
  "message": "Request succeeded",
  "data": {
    // Response data here
  }
}
```

For error responses:

```json
{
  "status": "error",
  "timestamp": "2025-04-21T17:45:00.000Z",
  "code": 400,
  "message": "Bad request",
  "errors": [
    {
      "type": "ValidationError",
      "message": "Invalid input data",
      "context": "Validating request data"
    }
  ]
}
```

## HTTP Status Codes

The API now uses standard HTTP status codes consistently across all endpoints:

| Status Code | Description | Usage |
|-------------|-------------|-------|
| 200 | OK | Successful request |
| 201 | Created | Resource created successfully |
| 204 | No Content | Request succeeded with no content to return |
| 400 | Bad Request | Invalid request parameters or data |
| 401 | Unauthorized | Authentication required |
| 403 | Forbidden | Permission denied |
| 404 | Not Found | Resource not found |
| 409 | Conflict | Request conflicts with current state of the resource |
| 429 | Too Many Requests | Rate limit exceeded |
| 500 | Internal Server Error | Server error |
| 503 | Service Unavailable | Service temporarily unavailable |

## Error Handling

All API endpoints now use standardized error handling:

1. Consistent error response format
2. Proper HTTP status codes for different error types
3. Detailed error messages with context
4. Proper logging of errors

## API Documentation

All API endpoints are now documented in OpenAPI 3.0.2 format. The documentation is available at:

- `/api/v1/docs/` - ReDoc UI
- `/api/v1/docs/swagger` - Swagger UI
- `/api/v1/docs/openapi.json` - OpenAPI specification in JSON format
- `/api/v1/docs/openapi.yaml` - OpenAPI specification in YAML format

## Implementation Details

### Standardization Utilities

The standardization is implemented using the following utilities:

1. `api/api_standards.py` - Core utilities for standardized responses and error handling
2. `api/api_decorators.py` - Decorators for API endpoints
3. `api/api_docs.py` - Utilities for API documentation
4. `api/openapi.py` - OpenAPI documentation generator

### Decorators

The following decorators can be used to standardize API endpoints:

#### `@standardize_response`

This decorator ensures that all responses from an endpoint follow the standardized format:

```python
from api.api_decorators import standardize_response

@standardize_response
def get(self):
    # Your code here
    return data, 200
```

#### `@validate_request_schema`

This decorator validates request data against a schema:

```python
from api.api_decorators import validate_request_schema
from marshmallow import Schema, fields

class MySchema(Schema):
    name = fields.String(required=True)
    age = fields.Integer(required=True)

@validate_request_schema(MySchema)
def post(self, name, age):
    # Your code here
    return {"name": name, "age": age}, 201
```

#### `@document_endpoint`

This decorator documents an endpoint in OpenAPI format:

```python
from api.api_docs import document_endpoint

@document_endpoint(
    summary="Get a resource",
    description="Get a resource by ID",
    tags=["Resources"]
)
def get(self, resource_id):
    # Your code here
    return data, 200
```

### Response Utilities

The following utilities can be used to create standardized responses:

#### `create_success_response`

Creates a standardized success response:

```python
from api.api_standards import create_success_response, jsonify_standard_response

def get(self):
    # Your code here
    return jsonify_standard_response(
        *create_success_response(
            data={"id": 1, "name": "Example"},
            message="Resource retrieved successfully",
            status_code=200
        )
    )
```

#### `create_error_response`

Creates a standardized error response:

```python
from api.api_standards import create_error_response, jsonify_standard_response

def get(self):
    try:
        # Your code here
    except Exception as e:
        return jsonify_standard_response(
            *create_error_response(
                error=e,
                status_code=500,
                context="Error retrieving resource"
            )
        )
```

## Audit and Standardization Tools

Two tools were created to help with the standardization process:

1. `api_audit.py` - Audits all API endpoints for consistency and generates a report
2. `apply_api_standardization.py` - Applies standardization to API endpoints based on the audit results

### Running the Audit

```bash
python api_audit.py --output-file api_audit_report.md --verbose
```

### Applying Standardization

```bash
python apply_api_standardization.py --audit-file api_audit_report.md --verbose
```

## Best Practices

When creating new API endpoints, follow these best practices:

1. Use the `@standardize_response` decorator on all endpoint methods
2. Use the `HTTPStatus` enum for all status codes
3. Use try-except blocks with standardized error handling
4. Document all endpoints using the `@document_endpoint` decorator
5. Validate request data using the `@validate_request_schema` decorator

## Example

Here's an example of a fully standardized API endpoint:

```python
from flask_restful import Resource
from http import HTTPStatus
from marshmallow import Schema, fields

from api.api_decorators import standardize_response, validate_request_schema
from api.api_docs import document_endpoint
from api.api_standards import create_success_response, create_error_response, jsonify_standard_response

class ExampleSchema(Schema):
    name = fields.String(required=True)
    description = fields.String(required=False)

class ExampleResource(Resource):
    """Resource for managing examples."""
    
    @document_endpoint(
        summary="Get an example",
        description="Get an example by ID",
        tags=["Examples"]
    )
    @standardize_response
    def get(self, example_id):
        """
        Get an example by ID.
        
        Args:
            example_id: ID of the example
            
        Returns:
            Standardized API response
        """
        try:
            # Your code here
            example = {"id": example_id, "name": "Example", "description": "This is an example"}
            
            return example, HTTPStatus.OK
        except Exception as e:
            return jsonify_standard_response(
                *create_error_response(
                    error=e,
                    context=f"Error retrieving example {example_id}"
                )
            )
    
    @document_endpoint(
        summary="Create an example",
        description="Create a new example",
        tags=["Examples"]
    )
    @validate_request_schema(ExampleSchema)
    @standardize_response
    def post(self, name, description=None):
        """
        Create a new example.
        
        Args:
            name: Name of the example
            description: Description of the example (optional)
            
        Returns:
            Standardized API response
        """
        try:
            # Your code here
            example = {"id": 1, "name": name, "description": description}
            
            return example, HTTPStatus.CREATED
        except Exception as e:
            return jsonify_standard_response(
                *create_error_response(
                    error=e,
                    context="Error creating example"
                )
            )