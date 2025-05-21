# API Module Reference

## Overview

The API module is the core of the CryoProtect application, providing RESTful endpoints for molecular analysis operations. It follows a standardized pattern for responses, error handling, and documentation.

## Key Components

### Resource Structure
- **Resource Classes**: Implemented using Flask-RESTful
- **Endpoint Registration**: Centralized in `__init__.py`
- **Authentication**: JWT-based with service role fallback
- **Documentation**: Auto-generated OpenAPI specs

### Core Resource Groups
- **Molecule Resources**: `/api/v1/molecules` - Molecular data management
- **Mixture Resources**: `/api/v1/mixtures` - Mixture creation and analysis
- **RDKit Resources**: `/api/v1/rdkit/*` - Molecular property calculation
- **Scoring Resources**: `/api/v1/scoring/*` - Molecular scoring algorithms
- **Protocol Designer**: `/api/v1/protocols/*` - Protocol creation and management
- **Lab Verification**: `/api/v1/verifications/*` - Experimental verification
- **Team Resources**: `/api/v1/teams/*` - Team and collaboration management
- **Export Resources**: `/api/v1/export/*` - Data export and sharing

### Standardization Patterns
- **Response Format**: All responses follow a standard envelope pattern
  ```json
  {
    "success": true,
    "data": {...},
    "errors": null,
    "metadata": {...}
  }
  ```
- **Error Handling**: Standardized error responses with context
- **Pagination**: Consistent pagination across list endpoints
- **Filtering**: Query parameter-based filtering

### Key Files
- `api/__init__.py`: Main API initialization and endpoint registration
- `api/resources.py`: Core resource implementations
- `api/schemas.py`: Request/response validation schemas
- `api/models.py`: Data transfer objects
- `api/api_standards.py`: Response standardization utilities
- `api/api_decorators.py`: Auth, validation, and other decorators
- `api/jwt_auth.py`: JWT authentication implementation
- `api/rbac.py`: Role-based access control

## Extension Points

When adding new API functionality:

1. **Create Resource Class**: Extend `Resource` and implement methods
2. **Define Schema**: Add validation schema in `schemas.py`
3. **Register Endpoints**: Add to API in `__init__.py`
4. **Add Documentation**: Document with apispec decorators
5. **Implement Tests**: Add tests in `tests/api/`

## Best Practices

1. **Follow Patterns**: Maintain consistent response formats
2. **Use Decorators**: Use the standard decorators for validation and auth
3. **Validate Input**: Always validate request data
4. **Transaction Safety**: Use connection context managers
5. **Proper Error Handling**: Use try/except with appropriate status codes
6. **Complete Documentation**: Document all endpoints

## Common Pitfalls

1. **Connection Management**: Always close connections properly
2. **Rate Limiting**: Be aware of rate limits for external services
3. **Auth Bypass**: Avoid bypassing authentication accidentally
4. **Response Size**: Large responses can cause performance issues
5. **N+1 Queries**: Optimize database queries to avoid multiple round trips